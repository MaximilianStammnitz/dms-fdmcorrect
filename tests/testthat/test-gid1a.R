## The bundled GID1A example (two libraries, thinned to 5 % of reads) against the outputs of
## nf-core/deepmutscan's error-correction modules on the same input. On Linux the files are
## byte-identical; the comparison here allows for last-digit differences of other maths libraries.

ext <- function(f) system.file("extdata", f, package = "fdmcorrect", mustWork = TRUE)

out <- file.path(tempdir(), "gid1a")
suppressMessages(fdm_run(c(ext("GID1A_input.variant_counts.tsv.gz"), ext("GID1A_output.variant_counts.tsv.gz")),
                         codon_library = ext("GID1A_codon_library.csv.gz"), reference = ext("GID1A.fasta"),
                         outdir = out, pool = "GID1A"))

expected <- test_path("expected")
files <- list.files(expected, recursive = TRUE)

test_that("every expected table is produced", {
    expect_length(files, 18)
    for (f in files) expect_true(file.exists(file.path(out, sub("\\.gz$", "", f))), info = f)
})

for (f in files) {
    test_that(paste("same as the pipeline:", f), {
        rn <- if (basename(f) == "fdm_seq_error_rate.csv.gz") 1 else NULL
        exp <- read.csv(file.path(expected, f), row.names = rn, check.names = FALSE, stringsAsFactors = FALSE)
        got <- read.csv(file.path(out, sub("\\.gz$", "", f)), row.names = rn, check.names = FALSE,
                        stringsAsFactors = FALSE)
        expect_equal(got, exp, tolerance = 1e-8)
    })
}

test_that("counts are corrected, never negative, and the originals kept", {
    x <- read.csv(file.path(out, "GID1A_input", "fdm_corrected_counts.csv"))
    snv <- !grepl(",", x$base_mut)
    expect_true(all(x$counts >= 0))
    expect_true(all(x$counts <= x$counts_raw))
    expect_true(any(x$counts[snv] < x$counts_raw[snv]))
    expect_identical(x$counts[!snv], x$counts_raw[!snv])
})

test_that("a pool with too little evidence is left uncorrected", {
    ev <- file.path(out, "GID1A_input", "fdm_evidence.csv")
    small <- read.csv(ev)
    small$F <- pmin(small$F, 0)
    d <- file.path(tempdir(), "small")
    dir.create(d, showWarnings = FALSE)
    write.csv(small, file.path(d, "fdm_evidence.csv"), row.names = FALSE)
    r <- suppressMessages(fdm_estimate(file.path(d, "fdm_evidence.csv"), id = "small"))
    expect_true(all(is.na(r$rates$e)))
})
