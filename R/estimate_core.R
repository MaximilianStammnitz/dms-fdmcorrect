# Ported unchanged from nf-core/deepmutscan (modules/local/dmsanalysis/fdm_estimate/templates/fdm_estimate.R);
# only the function names of the entry points differ. Do not edit without re-running the identity tests.

## nf-core/deepmutscan - false-double-mutant (FDM) error-rate estimation for one pool of libraries.
##
## Input: the per-library evidence tables of fdm_evidence() (F_j, M_j per SNV) of every
## library in the pool . F and M are summed per SNV across the pool, the
## error rate e_j is estimated once, and fdm_correct() applies the same
## e_j to every library of the pool (expected error reads = e_j * coverage of that library).
##
## Estimators :
## Also writes the same estimator applied to each library of the pool alone (fdm_rates_per_library.csv),
## which the error-correction report compares with the pooled rates.
##
##   mle  e_j = F_j / M_j. Undefined when M_j = 0 (no anchor in the window): reported as NA and the
##        variant is left uncorrected.
##   eb   Gamma-Poisson empirical Bayes. F_j ~ Poisson(M_j e_j), e_j ~ Gamma(nu mu, nu) with prior mean
##        mu and variance mu / nu, posterior mean e_j = (F_j + nu mu) / (M_j + nu). Priors form a
##        hierarchy :
##          substitution_type      global -> 6 pyrimidine-folded substitution types -> SNV
##          trinucleotide_context  global -> 6 types -> 96 trinucleotide contexts -> SNV
##        Each level's mean is its pooled rate shrunk towards the level above with that level's nu
##        as pseudo-exposure; each level's nu is fitted by maximising the negative-binomial marginal
##        likelihood of its SNVs (falling back to the parent's nu when a class has fewer than
##        `min_class_snvs` SNVs with exposure, or no false doubles at all). SNVs with M_j = 0 get the
##        prior mean of their class. SNVs at the first/last reference base have no trinucleotide
##        context and use their substitution-type prior.

NU_INTERVAL        <- c(1e-6, 1e12)  # search interval for nu (in exposure units, i.e. reads)
MIN_CLASS_SNVS     <- 10L            # fewer exposed SNVs in a class -> the class uses its parent's nu
MIN_POOL_FD_READS  <- 100            # fewer false-double reads in the whole pool -> no correction

## negative-binomial marginal of Gamma-Poisson with fixed prior mean mu, as a function of log(nu)
gp_negloglik <- function(log_nu, y, m, mu) {
    nu <- exp(log_nu)
    a  <- max(nu * mu, 1e-300)
    -sum(lgamma(y + a) - lgamma(a) - lgamma(y + 1) + a * log(nu / (nu + m)) + y * log(m / (nu + m)))
}

fit_nu <- function(y, m, mu, fallback) {
    ok <- m > 0
    y <- y[ok]; m <- m[ok]
    if (length(y) < MIN_CLASS_SNVS || sum(y) == 0 || !is.finite(mu) || mu <= 0)
        return(list(nu = fallback, fitted = FALSE, at_bound = NA))
    o <- try(optimize(gp_negloglik, interval = log(NU_INTERVAL), y = y, m = m, mu = mu), silent = TRUE)
    if (inherits(o, "try-error") || !is.finite(o[["objective"]]))
        return(list(nu = fallback, fitted = FALSE, at_bound = NA))
    at_bound <- min(abs(o[["minimum"]] - log(NU_INTERVAL))) < 1e-3 * diff(log(NU_INTERVAL))
    list(nu = exp(o[["minimum"]]), fitted = TRUE, at_bound = at_bound)
}

## one level of the hierarchy: classes `cls` shrunk to parent means `mu_parent`/`nu_parent` (per SNV)
fit_level <- function(level, cls, y, m, mu_parent, nu_parent) {
    ks <- sort(unique(cls[!is.na(cls)]))
    nk <- length(ks)
    mu_k <- nu_k <- Fk <- Mk <- numeric(nk); n_k <- n_exp <- integer(nk); fitted <- at_bound <- logical(nk)
    grp <- split(seq_along(cls), factor(cls, levels = ks))
    for (q in seq_len(nk)) {
        idx <- grp[[q]]
        mp <- mu_parent[idx[1]]; np <- nu_parent[idx[1]]
        Fk[q] <- sum(y[idx]); Mk[q] <- sum(m[idx])
        mu_k[q] <- max((Fk[q] + np * mp) / (Mk[q] + np), 1e-300)
        f <- fit_nu(y[idx], m[idx], mu_k[q], fallback = np)
        nu_k[q] <- f[["nu"]]; fitted[q] <- f[["fitted"]]; at_bound[q] <- f[["at_bound"]]
        n_k[q] <- length(idx); n_exp[q] <- sum(m[idx] > 0)
    }
    mu <- mu_parent; nu <- nu_parent
    has <- !is.na(cls)
    ki <- match(cls[has], ks)
    mu[has] <- mu_k[ki]
    nu[has] <- nu_k[ki]
    tab <- data.frame(level = rep(level, nk), class = ks, n_snv = n_k, n_snv_exposed = n_exp,
                      F = Fk, M = Mk, raw_rate = ifelse(Mk > 0, Fk / Mk, NA_real_),
                      mu = mu_k, nu = nu_k, nu_fitted = fitted, nu_at_bound = at_bound,
                      stringsAsFactors = FALSE)
    list(mu = mu, nu = nu, table = tab)
}

## ---- genuine double mutants ----------------------------------------------------------------------
## Two mutagenesis primers can integrate into the same molecule, so a read can carry a programmed anchor
## next to variant j without any error having occurred. Those reads are in F_j and are not sequencing
## error. kappa(d), the rate of genuine co-integration relative to independence at codon separation d, is
## measured on reads carrying TWO programmed multi-nucleotide variants (which no sequencing error can
## create), pooled over the libraries of this pool. The expected genuine-double reads for variant j are
##   G_j = f_j^true * sum_p kappa(d_jp) m_jp
## with f_j^true the library frequency of j minus its own error rate. kappa rises from ~0 while two primer
## footprints overlap to a plateau once they are clear, so each anchor must be weighted at its own
## distance; that is why the evidence step reports m_jp resolved by distance.

## exposure-weighted local mean of r over a window of +/- SMOOTH_HALFWIDTH codon distances. No shape is
## imposed: in nicking libraries kappa rises with distance until the two primer footprints clear each
## other, but in a tile-synthesised library (OCT1) it peaks around 20-25 codons and falls again, which a
## monotone fit would misrepresent.
SMOOTH_HALFWIDTH <- 2L

smooth_weighted <- function(r, w) {
    n <- length(r); out <- numeric(n)
    for (i in seq_len(n)) {
        j <- max(1L, i - SMOOTH_HALFWIDTH):min(n, i + SMOOTH_HALFWIDTH)
        tot <- sum(w[j])
        out[i] <- if (tot > 0) sum(r[j] * w[j]) / tot else r[i]
    }
    out
}

## pooled distance profiles -> kappa(d), one value per codon distance, clipped to [0, 1]
kappa_profile <- function(paths, max_dist) {
    d <- seq_len(max_dist)
    out <- data.frame(dist = d, gd_reads = 0, gd_exposure = 0, kappa_raw = NA_real_, kappa = 0)
    if (length(paths) == 0) return(out)
    p <- do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
    p <- p[is.finite(p[["gd_reads"]]) & is.finite(p[["gd_exposure"]]) & p[["dist"]] %in% d, , drop = FALSE]
    if (nrow(p) == 0) return(out)
    out[["gd_reads"]] <- as.numeric(tapply(p[["gd_reads"]], factor(p[["dist"]], levels = d), sum, default = 0))
    out[["gd_exposure"]] <- as.numeric(tapply(p[["gd_exposure"]], factor(p[["dist"]], levels = d), sum, default = 0))
    ok <- out[["gd_exposure"]] > 0
    if (!any(ok)) return(out)
    out[["kappa_raw"]][ok] <- out[["gd_reads"]][ok] / out[["gd_exposure"]][ok]
    r <- out[["kappa_raw"]]; r[!ok] <- 0
    k <- smooth_weighted(r, ifelse(ok, out[["gd_exposure"]], 0))
    ## a distance with no genuine-double exposure of its own still takes the window mean; only a distance
    ## whose whole window is empty falls back to the nearest value from below
    empty <- vapply(seq_along(k), function(i) {
        j <- max(1L, i - SMOOTH_HALFWIDTH):min(length(k), i + SMOOTH_HALFWIDTH); !any(ok[j]) }, TRUE)
    k[empty] <- NA_real_
    for (i in seq_along(k)) if (is.na(k[i])) k[i] <- if (i == 1L) 0 else k[i - 1L]
    out[["kappa"]] <- pmin(pmax(k, 0), 1)
    out
}

## pooled per-SNV kappa-weighted exposure KW_j = sum_p kappa(d_jp) m_jp
kappa_exposure <- function(paths, kappa, base_mut) {
    kw <- setNames(rep(0, length(base_mut)), base_mut)
    if (length(paths) == 0) return(unname(kw))
    for (f in paths) {
        e <- read.csv(f, stringsAsFactors = FALSE)
        e <- e[is.finite(e[["m"]]) & e[["dist"]] >= 1 & e[["dist"]] <= length(kappa), , drop = FALSE]
        if (nrow(e) == 0) next
        v <- tapply(kappa[e[["dist"]]] * e[["m"]], e[["base_mut"]], sum)
        j <- match(names(v), base_mut)
        ok <- !is.na(j)
        kw[j[ok]] <- kw[j[ok]] + as.numeric(v)[ok]
    }
    unname(kw)
}

## pooled evidence table (columns base_mut, type, context, F, M) -> per-SNV rates + prior table
## kw: kappa-weighted exposure per SNV; f_obs: the SNV's pooled library frequency. With kw = 0 (or the
## genuine-double subtraction switched off) the function is exactly the plain Gamma-Poisson fit.
estimate_core <- function(ev, method = "eb", prior = "substitution_type", kw = NULL, f_obs = NULL,
                         n_iter = 5L) {
    y <- ev[["F"]]; m <- ev[["M"]]
    if (is.null(kw)) kw <- rep(0, length(y))
    if (is.null(f_obs)) f_obs <- rep(0, length(y))
    kw <- ifelse(is.finite(kw), kw, 0); f_obs <- ifelse(is.finite(f_obs), f_obs, 0)
    est <- data.frame(base_mut = ev[["base_mut"]], type = ev[["type"]], context = ev[["context"]],
                      F = y, M = m, kappa_exposure = kw, G = 0, F_error = y,
                      prior_mean = NA_real_, prior_nu = NA_real_, shrinkage_weight = NA_real_,
                      e = NA_real_, stringsAsFactors = FALSE)
    usable <- is.finite(m) & is.finite(y)
    if (!any(usable & m > 0) || sum(y[usable]) < MIN_POOL_FD_READS) {
        ## too few false doubles in the whole pool to estimate an error rate: leave counts as they are
        return(list(rates = est, prior = NULL,
                    note = paste0("only ", sum(y[usable]), " false-double read(s) in the pool (minimum ",
                                  MIN_POOL_FD_READS, ") - no correction applied")))
    }
    mm <- ifelse(usable, m, 0)

    ## one pass of the hierarchy on the error reads yy
    fit_once <- function(yy) {
        mu_g <- max(sum(yy) / sum(mm), 1e-300)
        g <- fit_nu(yy, mm, mu_g, fallback = NA_real_)
        if (!g[["fitted"]]) g[["nu"]] <- 1 / mu_g   # degenerate pool: weak prior worth ~1 false double
        tab <- data.frame(level = "global", class = "all", n_snv = length(yy), n_snv_exposed = sum(mm > 0),
                          F = sum(yy), M = sum(mm), raw_rate = mu_g, mu = mu_g, nu = g[["nu"]],
                          nu_fitted = g[["fitted"]], nu_at_bound = g[["at_bound"]], stringsAsFactors = FALSE)
        mu <- rep(mu_g, length(yy)); nu <- rep(g[["nu"]], length(yy))
        lv <- fit_level("substitution_type", ev[["type"]], yy, mm, mu, nu)
        mu <- lv[["mu"]]; nu <- lv[["nu"]]; tab <- rbind(tab, lv[["table"]])
        if (prior == "trinucleotide_context") {
            lv <- fit_level("trinucleotide_context", ev[["context"]], yy, mm, mu, nu)
            mu <- lv[["mu"]]; nu <- lv[["nu"]]; tab <- rbind(tab, lv[["table"]])
        }
        list(mu = mu, nu = nu, table = tab, e = (yy + nu * mu) / (mm + nu))
    }

    ## Expectation-maximisation over the two sources of a false-double read. A read at variant j is a
    ## sequencing error with probability M_j e_j / (M_j e_j + G_j) and a genuine double integration
    ## otherwise, so the error reads are F_error_j = F_j * M_j e_j / (M_j e_j + G_j) (E step) and the
    ## hierarchy is refitted on those (M step). G_j itself depends on e_j through f^true = f_obs - e,
    ## so it is updated each round. This keeps the counts non-negative, which a plain F - G subtraction
    ## does not, and leaves the estimator untouched when G = 0.
    yy <- ifelse(usable, y, 0)
    use_gd <- any(kw > 0)
    G <- rep(0, length(yy))
    if (method == "mle") {
        e <- ifelse(mm > 0, yy / mm, 0)
        for (it in seq_len(if (use_gd) n_iter else 1L)) {
            if (!use_gd) break
            G <- pmax(f_obs - e, 0) * kw
            yy <- ifelse(mm * e + G > 0, y * (mm * e) / (mm * e + G), y)
            e <- ifelse(mm > 0, yy / mm, 0)
        }
        est[["G"]] <- G; est[["F_error"]] <- yy
        est[["e"]] <- ifelse(usable & m > 0, yy / mm, NA_real_)
        return(list(rates = est, prior = NULL, note = if (use_gd) "mle, genuine doubles removed" else "mle"))
    }
    fit <- fit_once(yy)
    if (use_gd) for (it in seq_len(n_iter)) {
        G <- pmax(f_obs - fit[["e"]], 0) * kw
        yy <- ifelse(mm * fit[["e"]] + G > 0, y * (mm * fit[["e"]]) / (mm * fit[["e"]] + G), y)
        yy <- ifelse(usable, yy, 0)
        fit <- fit_once(yy)
    }
    est[["G"]] <- G
    est[["F_error"]] <- ifelse(usable, yy, NA_real_)
    est[["prior_mean"]] <- fit[["mu"]]
    est[["prior_nu"]] <- fit[["nu"]]
    est[["shrinkage_weight"]] <- ifelse(usable, mm / (mm + fit[["nu"]]), NA_real_)  # weight on the SNV's own data
    est[["e"]] <- ifelse(usable, fit[["e"]], NA_real_)
    list(rates = est, prior = fit[["table"]],
         note = paste("eb", prior, if (use_gd) "- genuine doubles removed" else ""))
}

## sum per-library evidence tables per SNV
pool_evidence <- function(evs) {
    all <- do.call(rbind, lapply(evs, function(e) e[, c("base_mut", "type", "context", "F", "M", "counts", "cov")]))
    ## a library whose decay was not estimable contributes no evidence (NA), not zeros
    all <- all[is.finite(all[["F"]]) & is.finite(all[["M"]]), , drop = FALSE]
    keys <- unique(do.call(rbind, lapply(evs, function(e) e[, c("base_mut", "type", "context")])))
    if (nrow(all) == 0) {
        keys[["F"]] <- NA_real_; keys[["M"]] <- NA_real_
        keys[["counts"]] <- NA_real_; keys[["cov"]] <- NA_real_
        return(list(pooled = keys, n_libraries = length(evs), n_informative = 0L))
    }
    s <- aggregate(cbind(F, M, counts, cov) ~ base_mut, data = all, FUN = sum)
    pooled <- merge(keys, s, by = "base_mut", all.x = TRUE, sort = FALSE)
    n_inf <- sum(vapply(evs, function(e) any(is.finite(e[["M"]])), TRUE))
    list(pooled = pooled, n_libraries = length(evs), n_informative = n_inf)
}
