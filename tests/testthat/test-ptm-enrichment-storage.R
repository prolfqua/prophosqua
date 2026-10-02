test_that("the pipeline layout names one file per stage and analysis", {
  parameters <- list(
    run_kinase = TRUE,
    analyses = list(dpa = list(subdir = "PTM_DPA"), cf = list(subdir = "PTM_CF_DPU"))
  )
  files <- .ptm_enrichment_files(parameters, "out")
  expect_length(files, 10L)
  expect_identical(files[["PTMSEA__DPA"]], file.path("out", "PTM_DPA", "result_ptm_sea.json.gz"))
  expect_identical(
    files[["KinaseAssignments__CF"]],
    file.path("out", "PTM_CF_DPU", "intermediate_kinase_assignments.json.gz")
  )
  expect_length(.ptm_enrichment_files(modifyList(parameters, list(run_kinase = FALSE)), "out"), 0L)
})

test_that("assembly reads the enrichment files laid out beside the final MuData", {
  fixture <- ptm_enrichment_fixture()
  output <- file.path(fixture$root, "PTM_results_again.h5mu")
  expect_s3_class(assemble_ptm_results(fixture$statistics_path, output), "PTM_results")
  restored <- read_ptm_h5mu(output, PTM_results)
  expect_identical(names(restored$get_enrichments()), names(fixture$final$get_enrichments()))
  for (key in c("PTMSEA__DPA", "KinaseGSEA__DPU", "MEA__CF")) {
    payload <- jsonlite::fromJSON(protsea::read_gsea_json_text(fixture$files[[key]]), simplifyVector = FALSE)
    expect_named(payload, c("data", "rank_lists"))
  }
  elsewhere <- file.path(tempfile(), "PTM_results.h5mu")
  dir.create(dirname(elsewhere))
  expect_error(assemble_ptm_results(fixture$statistics_path, elsewhere), "Missing enrichment file")
})

test_that("the MEA is not computed in R", {
  fixture <- ptm_enrichment_fixture()
  expect_error(
    compute_ptm_enrichment(fixture$statistics_path, tempfile(), "MEA", "DPA"),
    "Unsupported computed enrichment stage"
  )
})
