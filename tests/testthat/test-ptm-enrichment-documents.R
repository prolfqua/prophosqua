enrichment_keys <- c(
  "PTMSEA__DPA",
  "PTMSEA__DPU",
  "PTMSEA__CF",
  "KinaseGSEA__DPA",
  "KinaseGSEA__DPU",
  "KinaseGSEA__CF",
  "MEA__DPA",
  "MEA__DPU",
  "MEA__CF"
)

test_that("the nine enrichments are protsea files beside the MuData, which only names them", {
  fixture <- ptm_enrichment_fixture()
  result <- fixture$final
  expect_identical(names(result$get_enrichments()), enrichment_keys)
  expect_identical(
    result$get_enrichment_document("PTMSEA", "dpa"),
    protsea::read_gsea_json_text(fixture$files[["PTMSEA__DPA"]])
  )
  expect_error(result$get_enrichment_document("STRING", "DPA"), "not enabled")
  expect_error(result$get_enrichment_document("PTMSEA", "unknown"), "not enabled")
  expect_error(result$get_enrichment_document("KinaseInputs", "DPA"), "not enabled")

  container <- prolfquapp::read_h5mu(fixture$path)
  namespace <- container$uns$prophosqua
  expect_setequal(names(namespace$enrichment_files), names(fixture$files))
  expect_setequal(names(namespace$enrichment_sha256), names(fixture$files))
  for (modality in container$modalities) {
    expect_false(any(c("enrichment_documents", "completed_stages") %in% names(modality$uns$prophosqua)))
  }
})

test_that("enrichments read back from their protsea files equal the computed ones", {
  skip_if_not_installed("enrichplot")
  fixture <- ptm_enrichment_fixture()
  restored <- fixture$final$get_enrichments()
  gs_info <- utils::getFromNamespace("gsInfo", "enrichplot")
  for (key in enrichment_keys) {
    before <- fixture$stages[[key]]$get_results()
    after <- restored[[key]]$get_results()
    expect_named(after, names(before))
    objects <- .PTM_RESULTS[[class(restored[[key]])[1L]]]$objects
    for (field in setdiff(names(before), objects)) {
      expect_equal(after[[field]], before[[field]], ignore_attr = TRUE, info = paste(key, field))
    }
    original <- before[[objects]]$a_vs_b
    decoded <- after[[objects]]$a_vs_b
    expect_equal(decoded@result, original@result, info = key)
    expect_equal(decoded@geneList, original@geneList, info = key)
    expect_equal(decoded@geneSets, original@geneSets, info = key)
    expect_equal(gs_info(decoded, geneSetID = 1), gs_info(original, geneSetID = 1), info = key)
    expect_silent(enrichplot::gseaplot2(decoded, geneSetID = 1))
  }
})

test_that("the final MuData refuses missing or changed enrichment files", {
  fixture <- ptm_enrichment_fixture()
  copy <- tempfile("ptm_results_copy_")
  dir.create(copy)
  file.copy(fixture$root, copy, recursive = TRUE)
  root <- file.path(copy, basename(fixture$root))
  path <- file.path(root, "PTM_results.h5mu")
  expect_s3_class(read_ptm_h5mu(path, PTM_results), "PTM_results")

  mea <- file.path(root, "DPA", "result_mea.json.gz")
  connection <- gzfile(mea, "w")
  writeLines("{}", connection)
  close(connection)
  expect_error(read_ptm_h5mu(path, PTM_results), "changed since the final MuData was written: MEA__DPA")
  unlink(mea)
  expect_error(read_ptm_h5mu(path, PTM_results), "missing: .*result_mea.json.gz")
})

test_that("kinase preparations must come from the statistics they are restored on", {
  fixture <- ptm_enrichment_fixture()
  inputs <- fixture$files[["KinaseInputs__DPA"]]
  statistics_hash <- .ptm_file_sha256(fixture$statistics_path)
  expect_s3_class(.read_ptm_preparation(inputs, "KinaseInputs", "DPA", statistics_hash)$seqwindows, "data.frame")
  expect_error(.read_ptm_preparation(inputs, "KinaseInputs", "DPA", strrep("0", 64)), "different statistics")
  expect_error(.read_ptm_preparation(inputs, "KinaseAssignments", "DPA", statistics_hash), "Wrong PTM CBOR stage")
  results <- PTM_results$new(fixture$final$get_statistics(), fixture$files, strrep("0", 64))
  expect_error(results$get_enrichments(), "different statistics")
})

test_that("enabled empty enrichments remain complete protsea documents", {
  fixture <- ptm_enrichment_fixture(empty = TRUE)
  enrichments <- fixture$final$get_enrichments()
  expect_length(enrichments, 9L)
  expect_equal(nrow(enrichments$PTMSEA__DPA$get_results()$results$a_vs_b@result), 0L)
  expect_equal(nrow(enrichments$KinaseGSEA__DPA$get_results()$gsea_results$a_vs_b@result), 0L)
  expect_equal(nrow(enrichments$MEA__DPA$get_results()$mea_clean), 0L)
})
