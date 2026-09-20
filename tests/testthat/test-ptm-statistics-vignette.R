test_that("statistics vignette input is a complete final MuData artifact", {
  path <- test_path("../../inst/extdata/ptm_results_example.h5mu")
  if (!file.exists(path)) {
    path <- system.file(
      "extdata",
      "ptm_results_example.h5mu",
      package = "prophosqua"
    )
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
  statistics <- result$get_statistics()
  expect_gte(length(statistics$get_inputs()$get_contrasts()), 2L)
  expect_gte(dplyr::n_distinct(result$get_tables()$DPA$site), 60L)
  expect_setequal(unique(result$get_tables()$DPA$modAA), c("S", "T", "Y"))

  dpa_dpu <- statistics$get_dpa_dpu()
  logo_specs <- list(
    DPA = list(
      data = dpa_dpu$combined_site_prot,
      fdr = "FDR.site",
      effect = "diff.site"
    ),
    DPU = list(
      data = dpa_dpu$combined_test_diff,
      fdr = "FDR_I",
      effect = "diff_diff"
    ),
    CF = list(
      data = statistics$get_cf()$results,
      fdr = "FDR.site",
      effect = "diff.site"
    )
  )
  for (analysis in names(logo_specs)) {
    spec <- logo_specs[[analysis]]
    sites <- filter_significant_sites(
      spec$data,
      fdr_col = spec$fdr,
      diff_col = spec$effect,
      fdr_threshold = 0.05,
      fc_threshold = 0.5,
      require_sequence = TRUE
    ) |>
      validate_sequence_window()
    expect_gte(dplyr::n_distinct(sites$contrast), 2L)
    for (contrast in unique(sites$contrast)) {
      expect_setequal(
        sites$regulation[sites$contrast == contrast],
        c("upregulated", "downregulated")
      )
    }
  }
})

test_that("statistics vignette has the required shallow tab structure", {
  path <- test_path("../../vignettes/ptm_statistics.qmd")
  if (!file.exists(path)) {
    path <- system.file("doc", "ptm_statistics.qmd", package = "prophosqua")
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
  sequence_starts <- which(source == "## Sequence logos")
  expect_length(sequence_starts, 3L)
  for (start in sequence_starts) {
    following_section <- which(
      seq_along(source) > start & grepl("^## ", source)
    )
    end <- if (length(following_section)) following_section[[1L]] - 1L else length(source)
    child_tabs <- sub(
      "^### ",
      "",
      grep("^### ", source[start:end], value = TRUE)
    )
    expect_identical(source[[start + 2L]], "::: {.panel-tabset}")
    expect_true(any(source[start:end] == ":::"))
    expect_identical(child_tabs, c("Plot", "Table"))
  }
  expect_equal(sum(source == "## Difference logos"), 3L)
  expect_true(any(grepl("plot_diff_logo", source, fixed = TRUE)))
  expect_true(any(grepl("read_ptm_h5mu", source, fixed = TRUE)))
  expect_false(any(grepl("read_excel|readRDS|read_xlsx", source)))
})
