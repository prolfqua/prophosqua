test_that("single-workbook statistics use the MuData tables", {
  paths <- anndata_pair_fixture()
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(paths$site),
    anndataR::read_h5ad(paths$protein),
    parameters = list(
      run_kinase = FALSE,
      analyses = list(dpa = list(subdir = "DPA"), dpu = list(subdir = "DPU"), cf = list(subdir = "CF"))
    )
  )
  statistics <- inputs$build(DPA_DPU)$build(PTM_statistics, cf = suppressWarnings(inputs$build(CF)))
  before <- statistics$get_dpa_dpu()
  tables <- .ptm_workbook_tables(PTM_results$new(statistics, list()))
  expect_setequal(names(tables), c(names(statistics$get_tables()), "CF_intensities", "CF_sample_annotation"))
  expect_equal(tables[names(statistics$get_tables())], statistics$get_tables())
  expect_equal(statistics$get_dpa_dpu(), before)
})

test_that("final H5MU exports one workbook with all nine enrichment tables", {
  input <- test_path("../../inst/extdata/ptm_results_example.h5mu")
  if (!file.exists(input)) {
    input <- system.file("extdata", "ptm_results_example.h5mu", package = "prophosqua")
  }
  expect_true(file.exists(input))
  result <- read_ptm_h5mu(input, PTM_results)
  output <- tempfile()
  workbook <- export_ptm_h5mu(input, output)
  expect_identical(workbook, file.path(output, "PTM_results.xlsx"))
  expect_identical(list.files(output, recursive = TRUE), "PTM_results.xlsx")

  sheets <- readxl::excel_sheets(workbook)
  expected <- c(names(result$get_tables()), "CF_intensities", "CF_sample_annotation")
  for (branch in result$get_enrichments()) {
    method <- class(branch)[1L]
    name <- paste(branch$get_analysis(), method, sep = "_")
    expected <- c(expected, name)
    if (method %in% c("KinaseGSEA", "MEA")) {
      expected <- c(expected, paste(name, "summary", sep = "_"))
    }
  }
  expect_identical(sheets, expected)
  expect_length(sheets, 23L)
  for (branch in result$get_enrichments()) {
    method <- class(branch)[1L]
    name <- paste(branch$get_analysis(), method, sep = "_")
    field <- c(PTMSEA = "all_clean", KinaseGSEA = "all_results", MEA = "mea_clean")[[method]]
    table <- readxl::read_xlsx(workbook, sheet = name)
    expect_equal(nrow(table), nrow(branch$get_results()[[field]]))
    expect_false(any(c("core_enrichment", "Leading.substrates") %in% names(table)))
  }
})
