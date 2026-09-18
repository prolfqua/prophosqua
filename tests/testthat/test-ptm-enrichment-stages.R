test_that("enabled enrichment collection determines final completeness", {
  paths <- anndata_pair_fixture()
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(paths$site),
    anndataR::read_h5ad(paths$protein),
    parameters = list(
      run_kinase = FALSE,
      analyses = list(dpa = list(subdir = "PTM_DPA"), dpu = list(subdir = "PTM_DPU"), cf = list(subdir = "PTM_CF_DPU"))
    )
  )
  statistics <- inputs$build(DPA_DPU)$build(PTM_statistics, cf = suppressWarnings(inputs$build(CF)))
  final <- PTM_results$new(statistics, list())
  path <- tempfile(fileext = ".h5mu")
  final$write_h5mu(path)
  restored <- read_ptm_h5mu(path, PTM_results)
  expect_length(restored$get_enrichments(), 0L)
  expect_length(restored$get_enrichment_documents(), 0L)
  expect_error(restored$get_enrichment_document("PTMSEA", "DPA"), "not enabled")
  expect_equal(restored$get_tables(), statistics$get_tables())
  expect_error(KinaseAssignments$new(statistics, "DPA", list(term2gene = data.frame())), "requires KinaseInputs")
})
