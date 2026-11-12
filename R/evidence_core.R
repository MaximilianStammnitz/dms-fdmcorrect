# Ported unchanged from nf-core/deepmutscan (modules/local/dmsanalysis/fdm_evidence/templates/fdm_evidence.R);
# only the function names of the entry points differ. Do not edit without re-running the identity tests.

## nf-core/deepmutscan - false-double-mutant (FDM) error evidence, one sequencing library.
##
## For every observed single-nucleotide library variant (SNV) j this collects the two numbers the
## FDM error model needs:
##   F_j  false-double reads: reads carrying SNV j together with a programmed multi-nucleotide codon
##        variant (an "anchor") within +/- max_FDdist codons. A library variant changes only
##        one codon, so such reads arise when a sequencing/PCR error creates j on an anchor molecule.
##   M_j  exposure: the number of anchor-bearing reads that span codon j, i.e. the number of
##        opportunities to observe that error. Per anchor codon p: A_p * FDcov_p / Acov_p, where A_p
##        is the anchor read count, Acov_p the anchor coverage and FDcov_p the joint coverage of
##        codons j and p - measured where a false double was observed, otherwise imputed from the
##        library-wide decay of joint coverage with codon distance.
##
## Estimation (pooling across libraries, empirical-Bayes shrinkage) happens downstream in
## fdm_estimate(), so that libraries can share evidence. Logic follows
## seq_error_correct_by_false_doubles_v5 (M. Stammnitz), vectorised, with three deliberate changes:
##   1. SNVs whose candidate false doubles all fail to match an anchor are KEPT with F = 0. The
##      previous template dropped them from the fit and left their counts uncorrected, which made the
##      fitted rate depend on sequencing depth (the F = 0 fraction grows as depth falls).
##   2. The maximum possible joint coverage of a false double is the NUMERIC minimum of the two
##      codons' median coverages (the previous code took the minimum of their character strings).
##   3. Every SNV is reported, including those with no anchor in the window (M = 0).
##   4. Anchors closer than `min_fd_dist` codons are not used. Within one mutagenesis primer's footprint the
##      false doubles are dominated by the primer's OWN recurrent synthesis/sequencing errors, which travel
##      with the anchor by construction and are not the per-site error the model estimates; a genuine second
##      primer integration is near-impossible there, since the two footprints would clash. Dropping those
##      anchors replaces the earlier 99.5 % count quantile, which removed the same reads but depended on
##      sequencing depth. The excluded region is summarised (near_*) and profiled by distance for the report.
##
## Implementation notes (same results as the first vectorised version, less time):
##   - the raw count table is read with data.table::fread, only the five columns this step uses;
##   - the reference is read with readLines (no Biostrings; loading it cost more than the whole analysis);
##   - the per-codon mismatch count, the (SNV, anchor) grid and its keys are vectorised instead of
##     built per SNV, and all group keys are integers rather than pasted strings.
##   Row order inside every group is unchanged, so sums are accumulated in the same order and the
##   outputs are byte-identical.

## first record of a FASTA file as a character vector of bases (as written, no case change)
read_first_fasta <- function(path) {
    l <- readLines(path, warn = FALSE)
    h <- which(startsWith(l, ">"))
    if (length(h) == 0) stop("no FASTA header line in ", path)
    end <- if (length(h) > 1) h[2] - 1 else length(l)
    s <- paste(trimws(l[seq.int(h[1] + 1, end)]), collapse = "")
    strsplit(s, "")[[1]]
}

## text before the first ":" as an integer (NA if there is no ":" or it is not a number)
codon_of <- function(cm) {
    p <- regexpr(":", cm, fixed = TRUE)
    suppressWarnings(as.integer(where(p > 0L, substr(cm, 1L, p - 1L), cm)))
}

## number of nucleotide changes in one "c:XXX>YYY" codon change (NA unless exactly one ">" separates
## two non-empty strings of equal length)
n_codon_changes <- function(cm) {
    p <- regexpr(":", cm, fixed = TRUE)
    x <- where(p > 0L, substr(cm, p + 1L, nchar(cm)), cm)
    n_gt <- nchar(x) - nchar(gsub(">", "", x, fixed = TRUE))
    q <- regexpr(">", x, fixed = TRUE)
    a <- where(q > 0L, substr(x, 1L, q - 1L), x)
    b <- where(q > 0L, substr(x, q + 1L, nchar(x)), x)
    ok <- n_gt == 1L & nchar(a) == nchar(b) & nchar(a) > 0L
    n <- integer(length(x))
    for (k in seq_len(max(0L, nchar(a[ok])))) n <- n + (substr(a, k, k) != substr(b, k, k))
    n[!ok] <- NA_integer_
    n
}

## pyrimidine-folded substitution type and trinucleotide context (96 classes) of SNV "pos:R>A"
snv_classes <- function(base_mut, wt_seq) {
    comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
    pos  <- as.integer(sub(":.*", "", base_mut))
    ref  <- substr(sub("^[^:]*:", "", base_mut), 1, 1)
    alt  <- substr(sub("^[^:]*:", "", base_mut), 3, 3)
    L    <- length(wt_seq)
    up   <- ifelse(pos > 1, wt_seq[pmax(pos - 1, 1)], NA_character_)
    dn   <- ifelse(pos < L, wt_seq[pmin(pos + 1, L)], NA_character_)
    flip <- ref %in% c("G", "A")
    ref2 <- ifelse(flip, comp[ref], ref)
    alt2 <- ifelse(flip, comp[alt], alt)
    up2  <- ifelse(flip, comp[dn], up)
    dn2  <- ifelse(flip, comp[up], dn)
    type <- paste0(ref2, ">", alt2)
    ctx  <- ifelse(is.na(up2) | is.na(dn2), NA_character_, paste0(up2, "[", type, "]", dn2))
    data.frame(type = unname(type), context = unname(ctx), stringsAsFactors = FALSE)
}

## numeric pairwise minimum (kept as a named helper so the change from the old code is explicit)
pair_min <- function(a, b) pmin(a, b)

## (row, anchor codon) pairs for every row i with codon j_i: the anchor codons p (in the order of
## `anc_codons`) with max(1, j - W) <= p <= min(max_codon, j + W), p != j. Row-major, so the pairs come
## out in the same order as a per-row loop would produce them.
window_grid <- function(j, anc_codons, W, max_codon) {
    uj <- unique(j)
    pl <- lapply(uj, function(c0) anc_codons[anc_codons >= max(1, c0 - W) & anc_codons <= min(max_codon, c0 + W) & anc_codons != c0])
    ui <- match(j, uj)
    len <- lengths(pl)[ui]
    list(i = rep(seq_along(j), len), p = unlist(pl[ui], use.names = FALSE))
}

## b, with a wherever cond is TRUE (ifelse without its copying and attribute handling)
where <- function(cond, a, b) { b[cond] <- if (length(a) == 1L) a else a[cond]; b }

## first and second codon change of a "c1, c2" codon_mut (NA second part if there is none)
split_pair <- function(cm) {
    p <- regexpr(", ", cm, fixed = TRUE)
    has2 <- p > 0L
    rest <- substr(cm, p + 2L, nchar(cm))
    q <- regexpr(",", rest, fixed = TRUE)
    list(cm1 = where(has2, substr(cm, 1L, p - 1L), cm),
         cm2 = where(has2, where(q > 0L, substr(rest, 1L, q - 1L), rest), rep(NA_character_, length(cm))))
}

## genuine double mutants: reads carrying two PROGRAMMED multi-nucleotide codon variants (2+2, 2+3, 3+3 nt
## changes). Sequencing/PCR error cannot plausibly create such a read, so these count real double-mutant
## molecules (synthesis, PCR chimeras). Per programmed multi-nt variant v (codon j, library frequency f_v):
##   D_v = reads carrying v and a programmed multi-nt variant at codon p, min_FDdist < |p - j| <= W
##   m_v = sum_p A_p * FDcov_vp / Acov_p   (the same anchor exposure as M for SNVs)
## and E[D_v] = kappa * f_v * m_v. Returns sum D and sum f m, so kappa can be pooled across libraries.
## The same count threshold as the false doubles applies, so kappa describes the genuine doubles that survive
## into F; E[genuine-double reads in F_j] = kappa * f_j * M_j (used by the error-signature fit of nf-core/deepmutscan).
genuine_doubles <- function(raw, anc, A, Acov, decay_pct, W, max_codon) {
    empty <- data.frame(dist = integer(0), gd_reads = numeric(0), gd_exposure = numeric(0),
                        gd_anchor_exposure = numeric(0))
    if (nrow(anc) == 0) return(empty)
    two <- raw[which(raw[["varying_codons"]] == 2 & raw[["varying_bases"]] >= 4), c("counts", "cov", "codon_mut"), drop = FALSE]
    sp <- split_pair(two[["codon_mut"]])
    cm1 <- sp[["cm1"]]; cm2 <- sp[["cm2"]]
    prog <- anc[["codon_mut"]]
    i1 <- match(cm1, prog); i2 <- match(cm2, prog)
    ok <- which(!is.na(i1) & !is.na(i2))
    nt <- nrow(anc)
    stride <- max_codon + 1
    t_i <- c(i1[ok], i2[ok])
    p_c <- c(anc[["codon"]][i2[ok]], anc[["codon"]][i1[ok]])
    cnt <- rep(two[["counts"]][ok], 2)
    cvr <- rep(two[["cov"]][ok], 2)
    keep <- abs(anc[["codon"]][t_i] - p_c) <= W
    key_pr <- factor((t_i[keep] - 1) * stride + p_c[keep])
    Dp <- tapply(cnt[keep], key_pr, sum)
    FDcovp <- tapply(cvr[keep], key_pr, max)
    Dk <- as.numeric(levels(key_pr))
    ac <- as.integer(names(A))
    g <- window_grid(anc[["codon"]], ac, W, max_codon)
    gi <- g[["i"]]; gp <- g[["p"]]
    if (length(gi) == 0) return(empty)
    j  <- anc[["codon"]][gi]
    hit <- match((gi - 1) * stride + gp, Dk)
    miss <- is.na(hit)
    pi <- match(gp, ac)
    Ap <- as.numeric(A)[pi]; Acp <- as.numeric(Acov)[pi]
    d <- abs(gp - j)
    fdc <- where(miss, round(decay_pct(d) / 100 * Acp), as.numeric(FDcovp)[hit])
    f <- anc[["counts"]][gi] / anc[["cov"]][gi]
    ## every read is counted for both of its variants; halve so that gd_reads counts reads
    D <- where(miss, 0, as.numeric(Dp)[hit]) / 2
    AE <- Ap * fdc / Acp / 2          # anchor-read exposure, as for the false doubles
    E <- f * AE                       # ... times the variant's library frequency: the independence expectation
    df <- factor(d)
    out <- data.frame(dist = as.integer(levels(df)),
                      gd_reads = as.numeric(tapply(D, df, sum)),
                      gd_exposure = as.numeric(tapply(E, df, sum)),
                      gd_anchor_exposure = as.numeric(tapply(AE, df, sum)))
    out[order(out[["dist"]]), , drop = FALSE]
}

evidence_core <- function(raw_path, filtered_path, wt_path, max_fd_dist, min_fd_dist) {

    W <- as.integer(max_fd_dist)
    wt_seq <- read_first_fasta(wt_path)

    ## the raw table has no header; only counts, cov, varying_bases, varying_codons and codon_mut are used
    raw <- data.table::fread(raw_path, sep = "\t", header = FALSE, quote = "", select = c(1L, 2L, 4L, 6L, 7L),
                 col.names = c("counts", "cov", "varying_bases", "varying_codons", "codon_mut"),
                 colClasses = list(integer = c(1L, 2L, 4L, 6L), character = 7L),
                 fill = TRUE, data.table = FALSE, showProgress = FALSE, nThread = 1L)
    proc <- read.csv(filtered_path, stringsAsFactors = FALSE)
    proc[["codon"]] <- codon_of(proc[["codon_mut"]])
    max_codon <- max(codon_of(raw[["codon_mut"]]), na.rm = TRUE)

    ## ---- observed SNVs ---------------------------------------------------------------------------
    snv <- proc[!grepl(",", proc[["base_mut"]], fixed = TRUE), c("base_mut", "codon_mut", "codon", "counts", "cov")]
    cls <- snv_classes(snv[["base_mut"]], wt_seq)
    snv[["type"]] <- cls[["type"]]
    snv[["context"]] <- cls[["context"]]
    snv[["F"]] <- 0
    snv[["M"]] <- 0
    snv[["n_anchor_codons"]] <- 0L
    snv[["n_fd_codons"]] <- 0L
    snv[["n_fd_candidates"]] <- 0L
    n_snv <- nrow(snv)

    info <- list(n_snv = n_snv, n_fd_raw = 0L, n_fd_used = 0L, min_fd_dist = as.integer(min_fd_dist),
                 decay_estimable = FALSE, decay_b0 = NA_real_, decay_b1 = NA_real_, decay_b2 = NA_real_,
                 near_reads = NA_real_, near_exposure = NA_real_, far_reads = NA_real_, far_exposure = NA_real_,
                 gd_reads = NA_real_, gd_exposure = NA_real_)
    profile <- data.frame(dist = integer(0), fd_reads = numeric(0), fd_exposure = numeric(0),
                          gd_reads = numeric(0), gd_exposure = numeric(0), gd_anchor_exposure = numeric(0))
    exposure0 <- data.frame(base_mut = character(0), dist = integer(0), m = numeric(0),
                            stringsAsFactors = FALSE)

    ## ---- false doubles: two changed codons, 3-4 changed bases, one side single-nt ------------------
    fd <- raw[which(raw[["varying_bases"]] >= 3 & raw[["varying_bases"]] < 5 & raw[["varying_codons"]] == 2), , drop = FALSE]
    info[["n_fd_raw"]] <- nrow(fd)
    sp <- split_pair(fd[["codon_mut"]])
    fd[["cm1"]] <- sp[["cm1"]]
    fd[["cm2"]] <- sp[["cm2"]]
    fd[["c1"]]  <- codon_of(fd[["cm1"]])
    fd[["c2"]]  <- codon_of(fd[["cm2"]])
    fd[["codon_dist"]] <- abs(fd[["c1"]] - fd[["c2"]])
    n1 <- n_codon_changes(fd[["cm1"]])
    n2 <- n_codon_changes(fd[["cm2"]])
    keep <- !is.na(n1) & !is.na(n2) & !(n1 == 1 & n2 == 1) & !(n1 == 2 & n2 == 2) &
            !(n1 == 3 & n2 == 2) & !(n1 == 2 & n2 == 3)
    fd <- fd[keep, , drop = FALSE]
    info[["n_fd_used"]] <- nrow(fd)

    ## ---- decay of joint coverage with codon distance ------------------------------------------------
    med_cov <- tapply(proc[["cov"]], proc[["codon"]], median)
    m1 <- unname(med_cov[as.character(fd[["c1"]])])
    m2 <- unname(med_cov[as.character(fd[["c2"]])])
    perc <- 100 * fd[["cov"]] / pair_min(m1, m2)
    usable <- is.finite(log(perc))
    if (nrow(fd) == 0 || length(unique(fd[["codon_dist"]][usable])) < 3) {
        snv[["M"]] <- NA_real_
        snv[["F"]] <- NA_real_
        return(list(evidence = snv, info = info, profile = profile, exposure = exposure0))
    }
    d  <- fd[["codon_dist"]][usable]
    fit <- lm(log(perc[usable]) ~ d + I(d^2))
    b <- unname(coef(fit))
    info[["decay_estimable"]] <- TRUE
    info[["decay_b0"]] <- b[1]; info[["decay_b1"]] <- b[2]; info[["decay_b2"]] <- b[3]
    decay_pct <- function(dist) exp(b[1] + b[2] * dist + b[3] * dist^2)

    ## ---- anchors: programmed multi-nt codon variants, aggregated per codon ----------------------------
    anc <- proc[proc[["varying_bases"]] > 1, , drop = FALSE]
    if (nrow(anc) == 0 || n_snv == 0) return(list(evidence = snv, info = info, profile = profile, exposure = exposure0))
    A    <- tapply(anc[["counts"]], anc[["codon"]], sum)
    Acov <- tapply(anc[["cov"]], anc[["codon"]], max)
    gd <- genuine_doubles(raw, anc, A, Acov, decay_pct, W, max_codon)
    info[["gd_reads"]] <- sum(gd[["gd_reads"]]); info[["gd_exposure"]] <- sum(gd[["gd_exposure"]])

    ## ---- false doubles attached to (SNV, anchor) pairs, both orientations ------------------------------
    ## (first every row with cm1 as the SNV, then every row with cm2 as the SNV, as before)
    pr_snv     <- c(fd[["cm1"]], fd[["cm2"]])
    pr_partner <- c(fd[["cm2"]], fd[["cm1"]])
    pr_snvcod  <- c(fd[["c1"]], fd[["c2"]])
    pr_counts  <- rep(fd[["counts"]], 2)
    pr_cov     <- rep(fd[["cov"]], 2)
    pr_dist    <- rep(fd[["codon_dist"]], 2)
    si <- match(pr_snv, snv[["codon_mut"]])          # SNV row of the single-nt side
    sel <- !is.na(si) & pr_dist <= W
    si <- si[sel]; pr_partner <- pr_partner[sel]; pr_snvcod <- pr_snvcod[sel]
    pr_counts <- pr_counts[sel]; pr_cov <- pr_cov[sel]; pr_dist <- pr_dist[sel]
    ## candidates: false doubles carrying the SNV within the window, whether or not the partner codon
    ## change is a programmed anchor (diagnostic only - the previous template skipped SNVs whose
    ## candidates all failed to match an anchor)
    snv[["n_fd_candidates"]] <- as.integer(tabulate(si[pr_dist > min_fd_dist], nbins = n_snv))
    ai <- match(pr_partner, anc[["codon_mut"]])      # anchor row of the partner side
    sel <- !is.na(ai)
    si <- si[sel]; ai <- ai[sel]; pr_snvcod <- pr_snvcod[sel]; pr_counts <- pr_counts[sel]; pr_cov <- pr_cov[sel]
    pr_p <- anc[["codon"]][ai]
    sel <- pr_p != pr_snvcod & pr_p >= pmax(1, pr_snvcod - W) & pr_p <= pmin(max_codon, pr_snvcod + W)
    si <- si[sel]; pr_p <- pr_p[sel]; pr_counts <- pr_counts[sel]; pr_cov <- pr_cov[sel]
    stride <- max_codon + 1
    key_pr <- factor((si - 1) * stride + pr_p)        # one key per (SNV, anchor codon)
    Fp     <- tapply(pr_counts, key_pr, sum)
    FDcovp <- tapply(pr_cov, key_pr, max)
    Fk     <- as.numeric(levels(key_pr))

    ## ---- every (SNV, anchor codon) pair in the window ----------------------------------------------------
    ## (SNV-major, anchor codons in the order of A, as before)
    anc_codons <- as.integer(names(A))
    g <- window_grid(snv[["codon"]], anc_codons, W, max_codon)
    gi <- g[["i"]]; gp <- g[["p"]]
    if (length(gi) == 0) return(list(evidence = snv, info = info, profile = profile, exposure = exposure0))
    c0 <- snv[["codon"]][gi]
    hit <- match((gi - 1) * stride + gp, Fk)
    miss <- is.na(hit)
    pi <- match(gp, anc_codons)
    g_A    <- as.numeric(A)[pi]
    g_Acov <- as.numeric(Acov)[pi]
    g_dist <- abs(gp - c0)
    g_F    <- where(miss, 0, as.numeric(Fp)[hit])
    g_FDcov <- where(miss, round(decay_pct(g_dist) / 100 * g_Acov), as.numeric(FDcovp)[hit])
    g_m    <- g_A * g_FDcov / g_Acov

    ## distance profile of the false-double evidence and of the genuine doubles, for the report
    df <- factor(g_dist)
    fdp <- data.frame(dist = as.integer(levels(df)),
                      fd_reads = as.numeric(tapply(g_F, df, sum)),
                      fd_exposure = as.numeric(tapply(g_m, df, sum)))
    profile <- merge(fdp, gd, by = "dist", all = TRUE)
    profile[is.na(profile)] <- 0
    profile <- profile[order(profile[["dist"]]), , drop = FALSE]

    ## anchors closer than min_fd_dist carry the primer's own recurrent errors: summarise, then drop
    near <- g_dist <= min_fd_dist
    info[["near_reads"]] <- sum(g_F[near]); info[["near_exposure"]] <- sum(g_m[near])
    info[["far_reads"]] <- sum(g_F[!near]); info[["far_exposure"]] <- sum(g_m[!near])
    far <- !near
    gi <- gi[far]; g_F <- g_F[far]; g_m <- g_m[far]; g_dist <- g_dist[far]; hit <- hit[far]
    if (length(gi) == 0) {
        snv[["M"]] <- 0; snv[["F"]] <- 0
        return(list(evidence = snv, info = info, profile = profile, exposure = exposure0))
    }

    gf <- factor(gi, levels = seq_len(n_snv))
    snv[["F"]] <- as.numeric(tapply(g_F, gf, sum, default = 0))
    snv[["M"]] <- as.numeric(tapply(g_m, gf, sum, default = 0))
    snv[["n_anchor_codons"]] <- as.integer(tabulate(gi, nbins = n_snv))
    snv[["n_fd_codons"]] <- as.integer(tabulate(gi[!is.na(hit)], nbins = n_snv))

    ## exposure resolved by anchor distance. fdm_estimate() needs it to weight each anchor by
    ## kappa at its own distance when it subtracts the genuine double mutants: M_j alone cannot, because
    ## kappa varies by an order of magnitude across the window.
    key_id <- gi * 1000L + g_dist
    md <- tapply(g_m, key_id, sum)
    ids <- as.numeric(names(md))
    exposure <- data.frame(base_mut = snv[["base_mut"]][as.integer(ids %/% 1000)], dist = as.integer(ids %% 1000),
                           m = as.numeric(md), stringsAsFactors = FALSE)
    exposure <- exposure[order(exposure[["base_mut"]], exposure[["dist"]]), , drop = FALSE]
    list(evidence = snv, info = info, profile = profile, exposure = exposure)
}
