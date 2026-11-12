CLI_USAGE <- "usage: dms-fdmcorrect <command> [options]

commands (options in brackets have defaults):
  run       all steps for a set of libraries that share one error process
            --counts FILE... --codon_library FILE --reference FILE
            [--ids NAME...] [--pool pool] [--outdir .] [estimate and evidence options]
  library   programmed library variants of one library
            --counts FILE --codon_library FILE --out FILE
  evidence  false-double evidence of one library
            --counts FILE --library_variants FILE --reference FILE --id NAME
            [--max_fd_dist 40] [--min_fd_dist 13] [--outdir .]
  estimate  error rates of one pool of libraries
            --evidence FILE... [--profiles FILE...] [--exposures FILE...] --id NAME
            [--method eb|mle] [--prior trinucleotide_context|substitution_type]
            [--max_fd_dist 40] [--no_genuine_doubles] [--wildtype] [--outdir .]
  correct   corrected counts of one library
            --library_variants FILE --evidence FILE --rates FILE --reference FILE
            [--outdir .] [--out fdm_corrected_counts.csv] [--seq_error_rate fdm_seq_error_rate.csv]

  --version   print the version
  --help      print this message

Variant count tables are in the format of dms_count (https://github.com/MaximilianStammnitz/dms-count).
"

## "--key v1 v2 --flag" -> list(key = c("v1", "v2"), flag = TRUE)
parse_cli_options <- function(args) {
    opts <- list()
    key <- NULL
    for (a in args) {
        if (startsWith(a, "--")) {
            key <- substring(a, 3)
            if (!nzchar(key)) stop("empty option name")
            opts[[key]] <- TRUE
        } else {
            if (is.null(key)) stop("unexpected argument '", a, "'")
            opts[[key]] <- if (isTRUE(opts[[key]])) a else c(opts[[key]], a)
        }
    }
    opts
}

#' Command-line interface
#'
#' The function behind the `dms-fdmcorrect` script. `dms-fdmcorrect --help` lists the commands.
#'
#' @param args Command-line arguments, e.g. `commandArgs(trailingOnly = TRUE)`.
#' @return Invisibly, the exit status (0 on success).
#' @export
fdm_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
    if (length(args) == 0 || args[1] %in% c("--help", "-h", "help")) {
        cat(CLI_USAGE)
        return(invisible(0L))
    }
    if (args[1] == "--version") {
        cat("dms-fdmcorrect", as.character(packageVersion("fdmcorrect")), "\n")
        return(invisible(0L))
    }
    cmd <- args[1]
    if (length(args) == 1L && cmd %in% c("run", "library", "evidence", "estimate", "correct")) {
        cat(CLI_USAGE)
        return(invisible(0L))
    }
    o <- parse_cli_options(args[-1])
    allowed <- list(
        run = c("counts", "codon_library", "reference", "ids", "pool", "outdir", "max_fd_dist", "min_fd_dist",
                "method", "prior", "no_genuine_doubles"),
        library = c("counts", "codon_library", "out"),
        evidence = c("counts", "library_variants", "reference", "id", "max_fd_dist", "min_fd_dist", "outdir"),
        estimate = c("evidence", "profiles", "exposures", "id", "method", "prior", "max_fd_dist",
                     "no_genuine_doubles", "wildtype", "outdir"),
        correct = c("library_variants", "evidence", "rates", "reference", "outdir", "out", "seq_error_rate"))
    if (!cmd %in% names(allowed)) stop("unknown command '", cmd, "' (see dms-fdmcorrect --help)")
    bad <- setdiff(names(o), allowed[[cmd]])
    if (length(bad)) stop("unknown option(s) for '", cmd, "': ", paste0("--", bad, collapse = " "))
    req <- function(k, n = 1L) {
        v <- o[[k]]
        if (is.null(v) || isTRUE(v)) stop("'", cmd, "' needs --", k)
        if (n == 1L && length(v) != 1L) stop("--", k, " takes one value")
        v
    }
    opt <- function(k, default) {
        v <- o[[k]]
        if (is.null(v)) default else if (isTRUE(v)) stop("--", k, " needs a value") else v
    }
    num <- function(k, default) as.numeric(opt(k, default))
    flag <- function(k) {
        v <- o[[k]]
        if (!is.null(v) && !isTRUE(v)) stop("--", k, " takes no value")
        isTRUE(v)
    }
    ## arguments are collected with list(), which evaluates them in order, so a missing option is
    ## reported in the order of the usage line
    call <- switch(cmd,
        run = list(fdm_run, list(req("counts", n = Inf), req("codon_library"), req("reference"),
                   outdir = opt("outdir", "."), ids = o[["ids"]], pool = opt("pool", "pool"),
                   max_fd_dist = num("max_fd_dist", 40), min_fd_dist = num("min_fd_dist", 13),
                   method = opt("method", "eb"), prior = opt("prior", "trinucleotide_context"),
                   subtract_genuine_doubles = !flag("no_genuine_doubles"))),
        library = list(fdm_library_variants, list(req("counts"), req("codon_library"), out = req("out"))),
        evidence = list(fdm_evidence, list(req("counts"), req("library_variants"), req("reference"), id = req("id"),
                        max_fd_dist = num("max_fd_dist", 40), min_fd_dist = num("min_fd_dist", 13),
                        outdir = opt("outdir", "."))),
        estimate = list(fdm_estimate, list(req("evidence", n = Inf), profiles = o[["profiles"]],
                        exposures = o[["exposures"]], id = req("id"), method = opt("method", "eb"),
                        prior = opt("prior", "trinucleotide_context"), max_fd_dist = num("max_fd_dist", 40),
                        subtract_genuine_doubles = !flag("no_genuine_doubles"), wildtype = flag("wildtype"),
                        outdir = opt("outdir", "."))),
        correct = list(fdm_correct, list(req("library_variants"), req("evidence"), req("rates"), req("reference"),
                       outdir = opt("outdir", "."), out = opt("out", "fdm_corrected_counts.csv"),
                       seq_error_rate = opt("seq_error_rate", "fdm_seq_error_rate.csv"))))
    do.call(call[[1]], call[[2]])
    invisible(0L)
}
