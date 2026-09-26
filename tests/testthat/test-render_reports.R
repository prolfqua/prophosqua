# ptm.sh starts an R process to find its own installation. R CMD check runs
# tests with the user library disabled and R_LIBS pointing at a library of its
# own, so that process cannot see prophosqua unless it is told where the copy
# under test lives -- outside check this adds the library the test already loaded
# the package from and changes nothing. R_TESTS is cleared for a similar reason:
# check sets it to a relative "startup.Rs" that a process started elsewhere
# cannot source. Both are artefacts of the harness; a caller in a work directory
# has neither. (The stub Rscript check puts on PATH is the wrapper's problem, not
# this one's: it uses $R_HOME/bin/Rscript when R_HOME is set.)
run_ptm_sh <- function(...) {
  lib <- dirname(system.file(package = "prophosqua"))
  suppressWarnings(system2(
    "bash",
    c(shQuote(system.file("application", "bin", "ptm.sh", package = "prophosqua")), ...),
    stdout = TRUE,
    stderr = TRUE,
    env = c(paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = ":")), "R_TESTS=")
  ))
}

test_that("report_file resolves every analysis report from the installed doc/", {
  # The analysis reports are the package's vignettes, so they are only installed
  # when the vignettes were built; an install that skipped them has no doc/ and
  # nothing to resolve. R CMD check --no-build-vignettes is exactly that case.
  skip_if(
    !nzchar(system.file("doc", "ptm_statistics.qmd", package = "prophosqua")),
    "package installed without vignettes built"
  )

  reports <- c("ptm_statistics.qmd", "ptm_enrichment.qmd")
  for (report in reports) {
    path <- report_file(report)
    expect_true(file.exists(path), info = report)
    expect_equal(basename(dirname(path)), "doc", info = report)
  }
})

test_that("report_file says what an install without vignettes is missing", {
  expect_error(
    report_file("Analysis_NoSuchReport.Rmd"),
    "installed with its vignettes built"
  )
})

test_that("copy_ptm_shell_script places one executable wrapper", {
  workdir <- tempfile("workdir")
  dir.create(workdir)
  copied <- suppressMessages(copy_ptm_shell_script(workdir))

  expect_equal(basename(copied), "ptm.sh")
  expect_true(file.access(copied, mode = 1) == 0)
})

test_that("ptm.sh help names every command script the package installs", {
  # The wrapper reads its command list from the installation, so this also says
  # that a newly added CMD_*.R needs no edit to the wrapper to be reachable.
  help <- run_ptm_sh("help")

  installed <- tolower(sub(
    "^CMD_(.*)\\.R$",
    "\\1",
    basename(list.files(system.file("application", package = "prophosqua"), pattern = "^CMD_.*\\.R$"))
  ))
  listed <- sub("^ +([a-z0-9_]+) +.*$", "\\1", grep("^ +[a-z0-9_]+ +\\S", help, value = TRUE))
  expect_true(length(listed) > 0)
  expect_setequal(listed, installed)

  # Each line carries the command's own first comment line as its summary.
  summaries <- sub("^ +[a-z0-9_]+ +", "", grep("^ +[a-z0-9_]+ +\\S", help, value = TRUE))
  expect_true(all(nchar(summaries) > 20))
})

test_that("ptm.sh refuses a command it does not have", {
  status <- run_ptm_sh("no_such_step")
  expect_equal(attr(status, "status"), 2L)
  expect_true(any(grepl("no such command", status)))
})
