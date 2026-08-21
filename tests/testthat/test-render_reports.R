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
    c(shQuote(application_file("bin/ptm.sh")), ...),
    stdout = TRUE,
    stderr = TRUE,
    env = c(paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = ":")), "R_TESTS=")
  ))
}

test_that("application_file resolves every report the pipeline renders", {
  reports <- c(
    "Analysis_DPA_DPU.Rmd",
    "Analysis_CorrectFirst_DEA.Rmd",
    "Analysis_n_to_c.Rmd",
    "Analysis_seqlogo.Rmd",
    "Analysis_PTMSEA.Rmd",
    "Analysis_KinaseLibrary.Rmd",
    "Analysis_MEA.Rmd",
    "_Overview_PhosphoAndIntegration_site.Rmd",
    "create_top_index.Rmd",
    "bibliography2025.bib"
  )
  for (report in reports) {
    expect_true(file.exists(application_file(report)), info = report)
  }
})

test_that("application_file resolves every command script and its wrapper", {
  scripts <- c(
    "CMD_DPA_DPU.R",
    "CMD_CF_DEA.R",
    "CMD_COMBINE_RESULTS.R",
    "CMD_PREP_PTMSIGDB.R",
    "CMD_PREP_KINASELIB.R",
    "CMD_DPU_OVERVIEW.R",
    "CMD_PTMSEA.R",
    "CMD_KINASELIB_GSEA.R",
    "CMD_MEA.R",
    "CMD_RENDER.R"
  )
  for (script in scripts) {
    expect_true(file.exists(application_file(script)), info = script)
  }

  expect_true(file.exists(application_file("bin/ptm.sh")))
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
    basename(list.files(dirname(application_file("CMD_RENDER.R")), pattern = "^CMD_.*\\.R$"))
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

test_that("application_file tells the caller to reinstall when a file is missing", {
  expect_error(
    application_file("Analysis_does_not_exist.Rmd"),
    "Reinstall the package"
  )
})
