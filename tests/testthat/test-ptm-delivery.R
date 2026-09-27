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
  statistics <- suppressWarnings(PTM_statistics$new(inputs))
  before <- statistics$get_dpa_dpu()
  tables <- .ptm_workbook_tables(PTM_results$new(statistics))
  expect_setequal(
    names(tables),
    c(names(statistics$get_tables()), "estimate_counts", "CF_intensities", "CF_sample_annotation")
  )
  expect_equal(tables[names(statistics$get_tables())], statistics$get_tables())
  expect_equal(statistics$get_dpa_dpu(), before)
})

test_that("tables keep observed site estimates unless all are asked for, and count both", {
  paths <- anndata_pair_fixture()
  inputs <- DEA_enriched_total$new(anndataR::read_h5ad(paths$site), anndataR::read_h5ad(paths$protein))
  computed <- suppressWarnings(PTM_statistics$new(inputs))
  dpa_dpu <- computed$get_dpa_dpu()
  cf <- computed$get_cf()
  imputed_site <- dpa_dpu$combined_site_prot$site[1L]
  for (table in c("combined_site_prot", "combined_test_diff")) {
    rows <- dpa_dpu[[table]]$site == imputed_site
    dpa_dpu[[table]]$estimate_type.site[rows] <- "lod_imputed"
  }
  reported <- cf$variants$correct_first_protein_imputed$results
  cf$variants$correct_first_protein_imputed$results$estimate_type[reported$site == imputed_site] <- "lod_imputed"
  statistics <- PTM_statistics$new(inputs, dpa_dpu = dpa_dpu, cf = cf)

  observed <- statistics$get_tables()
  all <- statistics$get_tables("all")
  for (analysis in c("DPA", "DPU", "CF")) {
    expect_true(all(observed[[analysis]]$estimate_type.site == "observed"))
    expect_true(any(all[[analysis]]$estimate_type.site == "lod_imputed"))
    expect_false(imputed_site %in% observed[[analysis]]$site)
  }

  expect_setequal(all$CF$site, cf$variants$correct_first_protein_imputed$results$site)
  expect_setequal(
    statistics$get_tables()$abundances_site_cf$site,
    colnames(cf$variants$correct_first_protein_imputed$abundances)
  )
  expect_setequal(statistics$get_cf_reported()$ptm_data$data_long()$site, all$CF$site)

  counts <- statistics$get_estimate_counts()
  for (analysis in c("DPA", "DPU", "CF")) {
    slice <- counts[counts$analysis == analysis, ]
    expect_identical(sum(slice$total), nrow(all[[analysis]]))
    expect_identical(sum(slice$observed), nrow(observed[[analysis]]))
    expect_identical(slice$total, slice$observed + slice$lod_imputed)
  }
  expect_error(observed_site_estimates(data.frame(site = "a")), "estimate_type.site")
})

test_that("final H5MU exports one workbook with all nine enrichment tables", {
  input <- test_path("../../inst/extdata/ptm_results_example/PTM_results.h5mu")
  if (!file.exists(input)) {
    input <- system.file("extdata", "ptm_results_example", "PTM_results.h5mu", package = "prophosqua")
  }
  expect_true(file.exists(input))
  result <- read_ptm_h5mu(input, PTM_results)
  output <- tempfile()
  workbook <- export_ptm_h5mu(input, output)
  expect_identical(workbook, file.path(output, "PTM_results.xlsx"))
  expect_identical(list.files(output, recursive = TRUE), "PTM_results.xlsx")

  tables <- .ptm_workbook_tables(result)
  sheets <- names(tables)
  expected <- c(names(result$get_tables()), "estimate_counts", "CF_intensities", "CF_sample_annotation")
  for (branch in result$get_enrichments()) {
    method <- class(branch)[1L]
    name <- paste(branch$get_analysis(), method, sep = "_")
    expected <- c(expected, name)
    if (method %in% c("KinaseGSEA", "MEA")) {
      expected <- c(expected, paste(name, "summary", sep = "_"))
    }
  }
  expect_identical(sheets, expected)
  expect_length(sheets, 24L)
  for (branch in result$get_enrichments()) {
    method <- class(branch)[1L]
    name <- paste(branch$get_analysis(), method, sep = "_")
    field <- c(PTMSEA = "all_clean", KinaseGSEA = "all_results", MEA = "mea_clean")[[method]]
    table <- tables[[name]]
    expect_equal(nrow(table), nrow(branch$get_results()[[field]]))
    expect_false(any(c("core_enrichment", "Leading.substrates") %in% names(table)))
  }
})
