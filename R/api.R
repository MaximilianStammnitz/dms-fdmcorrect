#' fdmcorrect: sequencing-error correction of deep mutational scanning counts by false double mutants
#'
#' In a library of programmed single-codon variants, a read that carries a programmed multi-nucleotide
#' codon variant (the "anchor") together with a single-nucleotide change elsewhere is, with rare
#' exceptions, the product of a sequencing or PCR error on the anchor molecule. Counting these false
#' double mutants (F) against the number of anchor reads that could have shown the error (the exposure
#' M) gives the error rate of every single-nucleotide variant (SNV), which is then subtracted from its
#' count. Rates are estimated per SNV with an empirical-Bayes Gamma-Poisson model whose priors follow
#' the trinucleotide context and substitution type, pooled across libraries, with genuine double
#' mutants (two programmed mutations in one molecule) removed from the evidence.
#'
#' The workflow has three steps, each also available from the command line (`dms-fdmcorrect`):
#' \enumerate{
#'   \item [fdm_evidence()] per library: F and M for every observed SNV;
#'   \item [fdm_estimate()] per pool of libraries: the error rate of every SNV;
#'   \item [fdm_correct()] per library: the corrected counts.
#' }
#' [fdm_library_variants()] prepares the per-library input from a variant count table, and [fdm_run()]
#' runs all steps for a set of libraries.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats aggregate coef lm median optimize setNames
#' @importFrom utils packageVersion read.csv read.table write.csv
NULL

DMS_COUNT_COLUMNS <- c("counts", "cov", "mean_length_variant_reads", "varying_bases",
                       "base_mut", "varying_codons", "codon_mut", "aa_mut", "pos_mut")

write_table <- function(x, outdir, name, row.names = FALSE) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(outdir, name)
    write.csv(x, file = path, row.names = row.names)
    invisible(path)
}

#' Programmed library variants of one sequencing library
#'
#' Reads a variant count table in the format of `dms_count` (and of GATK AnalyzeSaturationMutagenesis:
#' tab-separated, no header, columns counts, cov, mean_length_variant_reads, varying_bases, base_mut,
#' varying_codons, codon_mut, aa_mut, pos_mut) and keeps the single-codon variants that the mutagenesis
#' library programs, as nf-core/deepmutscan does before error correction.
#'
#' Kept are rows that change exactly one codon, to a codon listed for that position in `codon_library`
#' (or to the wild-type codon), with all changed bases within three consecutive positions. Rows without
#' an amino-acid annotation are dropped, synonymous rows get a `pos_mut` such as `A12A`, and the column
#' `counts_per_cov` (= counts / cov) is added.
#'
#' @param raw_counts Path of the variant count table (may be gzip-compressed).
#' @param codon_library Path of a CSV with columns `Codon_Number` (1 = first codon of the reading
#'   frame), `wt_codon` and `Variant`: one row per programmed codon at each position (e.g. the 32 NNK
#'   codons other than the wild type).
#' @param out Optional path of a CSV to write the result to.
#' @return A data frame, one row per programmed variant observed.
#' @export
fdm_library_variants <- function(raw_counts, codon_library, out = NULL) {
    counts <- read.table(raw_counts, sep = "\t", header = FALSE, fill = TRUE, col.names = DMS_COUNT_COLUMNS)
    counts <- counts[!(counts$aa_mut == "" | is.na(counts$aa_mut)), ]
    syn <- (is.na(counts$pos_mut) | counts$pos_mut == "") & grepl("^S:", counts$aa_mut)
    counts$pos_mut <- ifelse(syn,
                             paste0(sub("S:([A-Z])>[A-Z]", "\\1", counts$aa_mut),
                                    sub("^(\\d+):.*", "\\1", counts$codon_mut),
                                    sub("S:[A-Z]>([A-Z])", "\\1", counts$aa_mut)),
                             counts$pos_mut)
    counts$counts_per_cov <- counts$counts / counts$cov
    rownames(counts) <- NULL

    lib <- read.csv(codon_library)
    if (!all(c("Codon_Number", "wt_codon", "Variant") %in% colnames(lib)))
        stop("Codon library must contain columns 'Codon_Number', 'wt_codon', and 'Variant'.")
    counts$codon_mut <- toupper(counts$codon_mut)
    keys <- unique(c(paste(lib$Codon_Number, lib$Variant), paste(lib$Codon_Number, lib$wt_codon)))
    span <- vapply(strsplit(counts$base_mut, ",\\s*"), function(m) {
        p <- suppressWarnings(as.numeric(regmatches(m, regexpr("^[0-9]+", m))))
        if (length(p) == 0) p <- NA_real_
        suppressWarnings(max(p, na.rm = TRUE) - min(p, na.rm = TRUE))
    }, numeric(1))
    key <- paste(as.numeric(sub(":.*", "", counts$codon_mut)), sub(".*>", "", counts$codon_mut))
    keep <- counts$varying_codons == 1 & key %in% keys & span <= 2
    res <- counts[which(keep), , drop = FALSE]
    rownames(res) <- NULL
    if (!is.null(out)) write_table(res, dirname(out), basename(out))
    res
}

#' False-double evidence of one library
#'
#' For every observed SNV j: the false-double reads F_j (reads carrying j together with a programmed
#' multi-nucleotide codon variant between `min_fd_dist` and `max_fd_dist` codons away) and the exposure
#' M_j (the anchor reads spanning codon j). Also measures the genuine double mutants (reads with two
#' programmed multi-nucleotide variants) that [fdm_estimate()] removes from the evidence.
#'
#' @param raw_counts Path of the library's variant count table (format as in [fdm_library_variants()],
#'   may be gzip-compressed).
#'   It must contain the double-codon reads, i.e. the unfiltered table.
#' @param library_variants Path of the library's programmed variants (CSV written by
#'   [fdm_library_variants()]).
#' @param reference Path of the reference FASTA the reads were aligned to (first record used).
#' @param id Library name, written into the output tables.
#' @param max_fd_dist,min_fd_dist Codon distance window of the anchors (pipeline defaults 40 and 13).
#' @param outdir Optional directory to write `fdm_evidence.csv`, `fdm_evidence_summary.csv`,
#'   `fdm_distance_profile.csv` and `fdm_exposure_by_distance.csv` to.
#' @return A list of data frames: `evidence`, `summary`, `profile`, `exposure`.
#' @export
fdm_evidence <- function(raw_counts, library_variants, reference, id = "library",
                         max_fd_dist = 40, min_fd_dist = 13, outdir = NULL) {
    res <- evidence_core(raw_path = raw_counts, filtered_path = library_variants, wt_path = reference,
                         max_fd_dist = max_fd_dist, min_fd_dist = min_fd_dist)
    ev <- res[["evidence"]]
    ev[["library"]] <- id
    info <- as.data.frame(res[["info"]])
    info[["sample"]] <- id
    info[["n_snv_F0"]] <- sum(ev[["F"]] == 0, na.rm = TRUE)
    info[["n_snv_M0"]] <- sum(ev[["M"]] == 0, na.rm = TRUE)
    info[["F_total"]] <- sum(ev[["F"]], na.rm = TRUE)
    info[["M_total"]] <- sum(ev[["M"]], na.rm = TRUE)
    info[["gd_kappa"]] <- if (isTRUE(info[["gd_exposure"]] > 0)) info[["gd_reads"]] / info[["gd_exposure"]] else NA_real_
    info[["near_rate"]] <- if (isTRUE(info[["near_exposure"]] > 0)) info[["near_reads"]] / info[["near_exposure"]] else NA_real_
    info[["far_rate"]] <- if (isTRUE(info[["far_exposure"]] > 0)) info[["far_reads"]] / info[["far_exposure"]] else NA_real_
    prof <- res[["profile"]]
    prof[["library"]] <- id
    exp_d <- res[["exposure"]]
    exp_d[["library"]] <- id
    if (!res[["info"]][["decay_estimable"]])
        message("WARNING: fewer than 3 codon distances with usable false double mutants - the joint-coverage\n",
                "decay is not estimable, so this library contributes no error evidence. If it is pooled with\n",
                "other libraries it is corrected with their rates; on its own it is left uncorrected.")
    if (!is.null(outdir)) {
        write_table(ev, outdir, "fdm_evidence.csv")
        write_table(info, outdir, "fdm_evidence_summary.csv")
        write_table(prof, outdir, "fdm_distance_profile.csv")
        write_table(exp_d, outdir, "fdm_exposure_by_distance.csv")
    }
    invisible(list(evidence = ev, summary = info, profile = prof, exposure = exp_d))
}

#' Error rates of one pool of libraries
#'
#' Sums F and M per SNV over the pooled libraries and estimates every SNV's error rate e_j once:
#' `eb` (default) is a Gamma-Poisson empirical-Bayes estimate with hierarchical priors
#' (global, substitution type, and with `prior = "trinucleotide_context"` the 96 trinucleotide
#' contexts), `mle` is F_j / M_j. With `subtract_genuine_doubles`, the reads expected from genuine
#' double integrations, measured on the pool's reads with two programmed multi-nucleotide variants,
#' are removed from F by expectation-maximisation. Pools with fewer than 100 false-double reads are
#' not corrected.
#'
#' @param evidence Paths of the `fdm_evidence.csv` files of the libraries in the pool.
#' @param profiles,exposures Paths of their `fdm_distance_profile.csv` and
#'   `fdm_exposure_by_distance.csv`; by default the files of those names next to each evidence file.
#' @param id Pool name, written into the output tables.
#' @param method `"eb"` or `"mle"`.
#' @param prior `"trinucleotide_context"` or `"substitution_type"` (for `method = "eb"`).
#' @param max_fd_dist The `max_fd_dist` used by [fdm_evidence()].
#' @param subtract_genuine_doubles Remove genuine double mutants from the evidence.
#' @param wildtype `TRUE` for a library of the unmutated template, which has no anchors: the evidence
#'   is summed and no rates are estimated.
#' @param outdir Optional directory to write `fdm_rates.csv`, `fdm_prior_fit.csv`,
#'   `fdm_rates_per_library.csv` and `fdm_kappa_profile.csv` to.
#' @return A list of data frames: `rates` (per SNV, the rate is column `e`), `prior_fit`,
#'   `rates_per_library` (each library estimated alone, for comparison), `kappa_profile`.
#' @export
fdm_estimate <- function(evidence, profiles = NULL, exposures = NULL, id = "pool", method = c("eb", "mle"),
                         prior = c("trinucleotide_context", "substitution_type"), max_fd_dist = 40,
                         subtract_genuine_doubles = TRUE, wildtype = FALSE, outdir = NULL) {
    method <- match.arg(method)
    prior <- match.arg(prior)
    sibling <- function(name) {
        f <- file.path(dirname(evidence), name)
        f[file.exists(f)]
    }
    if (is.null(profiles)) profiles <- sibling("fdm_distance_profile.csv")
    if (is.null(exposures)) exposures <- sibling("fdm_exposure_by_distance.csv")
    profiles <- Filter(nzchar, profiles)
    exposures <- Filter(nzchar, exposures)

    evidence_tables <- lapply(evidence, read.csv, stringsAsFactors = FALSE)
    pool <- pool_evidence(evidence_tables)
    message("Pool '", id, "': ", pool[["n_libraries"]], " library/libraries, ", pool[["n_informative"]],
            " with estimable false-double evidence.")

    subtract_gd <- isTRUE(subtract_genuine_doubles)
    kp <- kappa_profile(if (subtract_gd) profiles else character(0), max_dist = max_fd_dist)
    kw <- if (subtract_gd && !isTRUE(wildtype))
              kappa_exposure(exposures, kp[["kappa"]], pool[["pooled"]][["base_mut"]]) else
              rep(0, nrow(pool[["pooled"]]))
    f_obs <- with(pool[["pooled"]], ifelse(is.finite(cov) & cov > 0, counts / cov, 0))
    if (any(kw > 0))
        message("Genuine doubles: kappa ", sprintf("%.3f", min(kp[["kappa"]])), " to ",
                sprintf("%.3f", max(kp[["kappa"]])), " over ", nrow(kp), " codon distances; ",
                sprintf("%.1f", 100 * sum(kp[["gd_reads"]]) / max(sum(kp[["gd_exposure"]]), 1e-300)),
                " % of independent co-occurrence overall.")

    if (isTRUE(wildtype)) {
        wt_rates <- pool[["pooled"]]
        wt_rates[["kappa_exposure"]] <- 0
        wt_rates[["G"]] <- 0
        wt_rates[["F_error"]] <- wt_rates[["F"]]
        wt_rates[["prior_mean"]] <- NA_real_; wt_rates[["prior_nu"]] <- NA_real_
        wt_rates[["shrinkage_weight"]] <- NA_real_; wt_rates[["e"]] <- NA_real_
        res <- list(rates = wt_rates,
                    prior = NULL, note = "wildtype-template library - no false-double correction applied")
    } else {
        res <- estimate_core(pool[["pooled"]], method = method, prior = prior, kw = kw, f_obs = f_obs)
    }
    message("Estimator: ", res[["note"]])
    kp[["pool"]] <- id
    rates <- res[["rates"]]
    rates[["pool"]] <- id

    per_lib <- lapply(seq_along(evidence), function(q) {
        f <- evidence[q]; e <- evidence_tables[[q]]
        lib <- if ("library" %in% names(e)) e[["library"]][1] else basename(dirname(f))
        one <- e[is.finite(e[["F"]]) & is.finite(e[["M"]]), c("base_mut", "type", "context", "F", "M"), drop = FALSE]
        r <- if (isTRUE(wildtype) || nrow(one) == 0) {
            data.frame(base_mut = e[["base_mut"]], F = e[["F"]], M = e[["M"]], e = NA_real_, shrinkage_weight = NA_real_)
        } else {
            estimate_core(one, method = method, prior = prior)[["rates"]][, c("base_mut", "F", "M", "e", "shrinkage_weight")]
        }
        r[["library"]] <- lib
        r
    })
    per_lib <- do.call(rbind, per_lib)
    per_lib[["pool"]] <- id
    pt <- res[["prior"]]
    if (is.null(pt)) {
        pt <- data.frame(level = character(0), class = character(0), n_snv = integer(0), n_snv_exposed = integer(0),
                         F = numeric(0), M = numeric(0), raw_rate = numeric(0), mu = numeric(0), nu = numeric(0),
                         nu_fitted = logical(0), nu_at_bound = logical(0))
    }
    if (nrow(pt) > 0) {
        pt[["pool"]] <- id
        pt[["n_libraries"]] <- pool[["n_libraries"]]
        pt[["median_shrinkage_weight"]] <- NA_real_
        w <- res[["rates"]][["shrinkage_weight"]]
        pt[["median_shrinkage_weight"]][1] <- median(w[res[["rates"]][["M"]] > 0], na.rm = TRUE)
    }
    if (!is.null(outdir)) {
        write_table(kp, outdir, "fdm_kappa_profile.csv")
        write_table(rates, outdir, "fdm_rates.csv")
        write_table(per_lib, outdir, "fdm_rates_per_library.csv")
        write_table(pt, outdir, "fdm_prior_fit.csv")
    }
    invisible(list(rates = rates, prior_fit = pt, rates_per_library = per_lib, kappa_profile = kp))
}

#' Corrected counts of one library
#'
#' Subtracts the expected error reads e_j * coverage_j from every SNV that has a rate (rounded, floored
#' at zero). Multi-nucleotide codon variants are not corrected, as they need two or three independent
#' errors.
#'
#' @param library_variants Path of the library's programmed variants (CSV from
#'   [fdm_library_variants()]).
#' @param evidence Path of the library's `fdm_evidence.csv`.
#' @param rates Path of the `fdm_rates.csv` of the library's pool.
#' @param reference Path of the reference FASTA. Positions and wild-type bases of the SNVs are checked
#'   against it.
#' @param outdir Optional directory to write the two tables to.
#' @param out,seq_error_rate File names of the two tables in `outdir`.
#' @return A list: `counts`, the library variants with `counts` and `counts_per_cov` corrected and the
#'   originals kept as `counts_raw` and `counts_per_cov_raw`; and `seq_error_rate`, one row per possible
#'   SNV of the reference with its evidence, rate and prior.
#' @export
fdm_correct <- function(library_variants, evidence, rates, reference, outdir = NULL,
                        out = "fdm_corrected_counts.csv", seq_error_rate = "fdm_seq_error_rate.csv") {
    check_reference_frame(reference, library_variants)
    r <- correct_core(wt_path = reference, filtered_path = library_variants, evidence_path = evidence,
                      rates_path = rates)
    message("FDM correction: ", r[["n_corrected"]], " of ", r[["n_snv"]], " single-nucleotide variants corrected; ",
            r[["n_snv"]] - r[["n_corrected"]], " without an estimable rate left unchanged.")
    ic <- r[["counts"]]
    ic$counts_raw <- ic$counts
    ic$counts_per_cov_raw <- ic$counts_per_cov
    ic$counts <- ic$counts_corrected
    ic$counts_per_cov <- ic$counts_per_cov_corrected
    if (!is.null(outdir)) {
        write_table(ic, outdir, out)
        write_table(r[["seq_error_rate"]], outdir, seq_error_rate, row.names = TRUE)
    }
    invisible(list(counts = ic, seq_error_rate = r[["seq_error_rate"]]))
}

#' All steps for a set of libraries
#'
#' Runs [fdm_library_variants()] and [fdm_evidence()] for every library, [fdm_estimate()] once for all
#' of them as one pool, and [fdm_correct()] for every library. Pooling assumes one error process, i.e.
#' libraries prepared, processed and sequenced together; for anything else, run the steps separately.
#'
#' @param raw_counts Paths of the libraries' variant count tables.
#' @param codon_library,reference As in [fdm_library_variants()] and [fdm_evidence()].
#' @param outdir Output directory: one subdirectory per library (named by `ids`) and one for the pool.
#' @param ids Library names; by default the file names of `raw_counts` up to the first dot.
#' @param pool Name of the pool and its subdirectory.
#' @param ... Passed to [fdm_estimate()] (`method`, `prior`, `subtract_genuine_doubles`).
#' @inheritParams fdm_evidence
#' @return Invisibly, the paths of the corrected count tables.
#' @export
fdm_run <- function(raw_counts, codon_library, reference, outdir = ".", ids = NULL, pool = "pool",
                    max_fd_dist = 40, min_fd_dist = 13, ...) {
    if (is.null(ids)) ids <- sub("\\..*$", "", basename(raw_counts))
    if (length(ids) != length(raw_counts) || anyDuplicated(ids) || pool %in% ids)
        stop("library names must be unique, one per count table, and differ from the pool name")
    lv <- file.path(outdir, ids, "fdm_library_variants.csv")
    for (q in seq_along(raw_counts)) {
        fdm_library_variants(raw_counts[q], codon_library, out = lv[q])
        fdm_evidence(raw_counts[q], lv[q], reference, id = ids[q], max_fd_dist = max_fd_dist,
                     min_fd_dist = min_fd_dist, outdir = file.path(outdir, ids[q]))
    }
    fdm_estimate(file.path(outdir, ids, "fdm_evidence.csv"), id = pool, max_fd_dist = max_fd_dist,
                 outdir = file.path(outdir, pool), ...)
    for (q in seq_along(raw_counts))
        fdm_correct(lv[q], file.path(outdir, ids[q], "fdm_evidence.csv"), file.path(outdir, pool, "fdm_rates.csv"),
                    reference, outdir = file.path(outdir, ids[q]))
    invisible(file.path(outdir, ids, "fdm_corrected_counts.csv"))
}
