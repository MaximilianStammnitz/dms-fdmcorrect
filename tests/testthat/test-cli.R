test_that("options are parsed", {
    o <- fdmcorrect:::parse_cli_options(c("--evidence", "a.csv", "b.csv", "--wildtype", "--id", "p"))
    expect_identical(o$evidence, c("a.csv", "b.csv"))
    expect_true(o$wildtype)
    expect_identical(o$id, "p")
    expect_error(fdmcorrect:::parse_cli_options("stray"))
})

test_that("unknown commands and options are refused", {
    expect_error(fdm_cli(c("nonsense")), "unknown command")
    expect_error(fdm_cli(c("correct", "--colour", "red")), "unknown option")
    expect_error(fdm_cli(c("evidence", "--id", "x")), "needs --counts")
})

test_that("--version prints the package version", {
    expect_output(fdm_cli("--version"), as.character(packageVersion("fdmcorrect")), fixed = TRUE)
})

test_that("a command without options prints the usage", {
    expect_output(fdm_cli("estimate"), "usage: dms-fdmcorrect")
})
