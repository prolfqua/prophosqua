test_that("statistics vignette input is a complete final MuData artifact", {
  path <- system.file(
    "extdata",
    "ptm_results_example.h5mu",
    package = "prophosqua"
  )
  if (!nzchar(path)) {
    path <- test_path("../../inst/extdata/ptm_results_example.h5mu")
  }
  expect_true(file.exists(path))
  result <- read_ptm_h5mu(path, PTM_results)
  expect_s3_class(result, "PTM_results")
  expect_length(result$get_enrichment_documents(), 0L)
  expect_named(
    result$get_tables(),
    c(
      "DPA",
      "DPU",
      "CF",
      "abundances_protein",
      "abundances_site_dpa",
      "abundances_site_cf"
    )
  )
  expect_gt(nrow(result$get_tables()$DPA), 0L)
  expect_gt(nrow(result$get_tables()$DPU), 0L)
  expect_gt(nrow(result$get_tables()$CF), 0L)
})

test_that("statistics vignette has the required shallow tab structure", {
  path <- system.file("doc", "ptm_statistics.qmd", package = "prophosqua")
  if (!nzchar(path)) {
    path <- test_path("../../vignettes/ptm_statistics.qmd")
  }
  source <- readLines(path, warn = FALSE)
  fence <- grepl("^```", source)
  outside_code <- cumsum(fence) %% 2L == 0L & !fence
  top_tabs <- sub(
    "^# ",
    "",
    grep("^# ", source[outside_code], value = TRUE)
  )
  expect_identical(
    top_tabs,
    c("Overview", "DPA", "DPU", "CorrectFirst DPU", "Session Info")
  )
  session_start <- match("# Session Info", source)
  session_source <- source[session_start:length(source)]
  session_fence <- grepl("^```", session_source)
  session_outside_code <- cumsum(session_fence) %% 2L == 0L & !session_fence
  session_tabs <- sub(
    "^## ",
    "",
    grep("^## ", session_source[session_outside_code], value = TRUE)
  )
  expect_identical(session_tabs, c("Report provenance", "R session info"))
  expect_true(any(grepl("read_ptm_h5mu", source, fixed = TRUE)))
  expect_false(any(grepl("read_excel|readRDS|read_xlsx", source)))
})
