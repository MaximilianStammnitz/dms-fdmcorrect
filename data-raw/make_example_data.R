## Makes inst/extdata from two full-depth GID1A libraries (input and output, replicate 1), thinned to
## 5 % of their reads so that the example stays small: every variant's count is drawn from
## Binomial(count, 0.05) and its coverage scaled by 0.05 (at least the new count); variants left with
## no read are dropped.
## usage: Rscript data-raw/make_example_data.R <input.variant_counts.tsv> <output.variant_counts.tsv> <outdir>
args <- commandArgs(trailingOnly = TRUE)
thin <- function(path, out, p = 0.05) {
    x <- read.table(path, sep = "\t", header = FALSE, quote = "", comment.char = "", fill = TRUE,
                    colClasses = "character", col.names = paste0("V", 1:9))
    n <- as.integer(x$V1)
    k <- rbinom(length(n), n, p)
    cv <- pmax(round(as.integer(x$V2) * p), k)
    x$V1 <- k; x$V2 <- cv
    x <- x[k > 0, , drop = FALSE]
    con <- gzfile(out, "w")
    write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(con)
}
set.seed(20261112)
thin(args[1], file.path(args[3], "GID1A_input.variant_counts.tsv.gz"))
thin(args[2], file.path(args[3], "GID1A_output.variant_counts.tsv.gz"))
