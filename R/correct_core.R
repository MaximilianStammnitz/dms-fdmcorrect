# Ported unchanged from nf-core/deepmutscan (modules/local/dmsanalysis/error_correction_false_doubles/templates/error_correction_false_doubles.R);
# changes: entry-point names, and apply_fdm_rates returns its two tables instead of writing them. Do not edit without re-running the identity tests.

## nf-core/deepmutscan - false-double-mutant (FDM) sequencing-error correction, one library.
##
## Applies the per-SNV error rates e_j estimated by fdm_estimate() for this library's pool:
## expected error reads = e_j * coverage_j of this library, subtracted from the SNV's count (rounded,
## floored at zero). Multi-nucleotide codon variants are not corrected - a two- or three-nucleotide
## change needs two or three independent errors and is negligibly affected.
##
## SNVs without a rate (MLE with no anchor in the window, or no usable evidence) are left as they
## are. All SNVs with a rate are corrected, including those with no observed false double (F = 0).

correct_core <- function(wt_path, filtered_path, evidence_path, rates_path) {

  wt <- read_first_fasta(wt_path)
  nts <- c("A", "T", "G", "C")
  seq.error.names <- unlist(lapply(seq_along(wt), function(i) paste0(i, ":", wt[i], ">", setdiff(nts, wt[i]))))

  input.counts.processed <- read.csv(filtered_path)
  input.counts.processed$counts_corrected <- input.counts.processed$counts
  input.counts.processed$counts_per_cov_corrected <- input.counts.processed$counts_per_cov

  ev <- read.csv(evidence_path, stringsAsFactors = FALSE)
  rates <- read.csv(rates_path, stringsAsFactors = FALSE)

  ## per-SNV table in reference order (same first three columns as before for compatibility)
  n <- length(seq.error.names)
  ser <- data.frame(check.names = FALSE,
                    "1nt counts/coverage" = rep(0, n), "1nt false counts/coverage" = rep(0, n),
                    "false double reads (F)" = rep(NA_real_, n), "exposure (M)" = rep(NA_real_, n),
                    "pooled F" = rep(NA_real_, n), "pooled M" = rep(NA_real_, n),
                    "prior mean" = rep(NA_real_, n), "prior nu" = rep(NA_real_, n),
                    "shrinkage weight" = rep(NA_real_, n))
  rownames(ser) <- seq.error.names
  snv.rows <- grep("[,]", input.counts.processed$base_mut, invert = TRUE)
  snv <- input.counts.processed$base_mut[snv.rows]
  ser[snv, "1nt counts/coverage"] <- input.counts.processed$counts_per_cov[snv.rows]
  ie <- match(snv, ev$base_mut)
  ser[snv, "false double reads (F)"] <- ev$F[ie]
  ser[snv, "exposure (M)"] <- ev$M[ie]
  ir <- match(snv, rates$base_mut)
  e <- rates$e[ir]
  ser[snv, "1nt false counts/coverage"] <- ifelse(is.na(e), 0, e)
  ser[snv, "pooled F"] <- rates$F[ir]
  ser[snv, "pooled M"] <- rates$M[ir]
  ser[snv, "prior mean"] <- rates$prior_mean[ir]
  ser[snv, "prior nu"] <- rates$prior_nu[ir]
  ser[snv, "shrinkage weight"] <- rates$shrinkage_weight[ir]

  ## correct every SNV that has a rate
  has <- !is.na(e)
  r <- snv.rows[has]
  input.counts.processed$counts_per_cov_corrected[r] <- pmax(0, input.counts.processed$counts_per_cov[r] - e[has])
  cc <- round(input.counts.processed$counts[r] - e[has] * input.counts.processed$cov[r])
  neg <- cc < 0
  cc[neg] <- 0
  input.counts.processed$counts_corrected[r] <- cc
  input.counts.processed$counts_per_cov_corrected[r[neg]] <- 0

  list(counts = input.counts.processed, seq_error_rate = ser, n_snv = length(snv.rows), n_corrected = sum(has))
}

# Defensive coordinate-contract check (general, not GID1A-specific): base_mut positions must be
# 1-based within the supplied reference, and the stated WT base must match the reference base there.
# A mismatch means the reference and the counts are in different frames -> stop rather than
# silently mis-index seq.error.rate.
check_reference_frame <- function(wt_path, processed_path){
  ref <- read_first_fasta(wt_path)
  proc <- read.csv(processed_path)
  snv <- proc[grep("[,]", proc[,"base_mut"], invert = TRUE), "base_mut"]
  if(length(snv) == 0) return(invisible(NULL))
  pos <- suppressWarnings(as.integer(sub(":.*", "", snv)))
  wtbase <- substr(sub("^[^:]*:", "", snv), 1, 1)
  ok <- !is.na(pos) & pos >= 1 & pos <= length(ref)
  if(any(!ok)){
    stop("Reference/counts coordinate mismatch: ", sum(!ok), " single-nucleotide base_mut positions ",
         "fall outside the reference (length ", length(ref), "). The reference fasta must be the same ",
         "one the reads were aligned to, in the same coordinate frame as base_mut.")
  }
  mism <- ref[pos] != wtbase
  if(any(mism, na.rm = TRUE)){
    stop("Reference/counts coordinate mismatch: ", sum(mism, na.rm = TRUE), " single-nucleotide ",
         "base_mut entries state a wild-type base that disagrees with the reference at that position ",
         "(e.g. ", snv[which(mism)[1]], "). Check that the reference matches the counts.")
  }
  invisible(NULL)
}
