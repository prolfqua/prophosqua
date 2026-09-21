test_that("compact CBOR stages assemble the same nine final JSON documents", {
  fixture <- ptm_enrichment_fixture()
  final <- fixture$final
  statistics_path <- tempfile(fileext = ".h5mu")
  output_path <- tempfile(fileext = ".h5mu")
  final$get_statistics()$write_h5mu(statistics_path)
  hash <- .ptm_statistics_hash(statistics_path)
  paths <- character()
  for (analysis in c("DPA", "DPU", "CF")) {
    branches <- final$get_enrichments()
    ptmsea <- branches[[paste0("PTMSEA__", analysis)]]
    kinase <- branches[[paste0("KinaseGSEA__", analysis)]]
    mea <- branches[[paste0("MEA__", analysis)]]
    assignments <- kinase$get_source()
    inputs <- assignments$get_source()
    motif <- mea$get_source()
    for (stage in list(ptmsea, inputs, assignments, motif, kinase, mea)) {
      path <- tempfile(fileext = ".cbor")
      .write_ptm_cbor(stage, path, hash)
      paths <- c(paths, path)
    }
  }
  assembled <- assemble_ptm_cbor(statistics_path, paths, output_path)
  expect_s3_class(assembled, "PTM_results")
  expect_identical(names(assembled$get_enrichment_documents()), names(final$get_enrichment_documents()))
  for (key in names(final$get_enrichment_documents())) {
    actual <- assembled$get_enrichment_documents()[[key]]
    expected <- final$get_enrichment_documents()[[key]]
    expect_setequal(names(actual), names(expected))
    for (field in names(expected)) {
      expect_identical(actual[[field]], expected[[field]])
    }
  }
  expect_error(assemble_ptm_cbor(statistics_path, paths[-1], output_path), "every enabled CBOR stage")
  bad_path <- tempfile(fileext = ".cbor")
  artifact <- secretbase::cbordec(readBin(paths[[1]], what = "raw", n = file.info(paths[[1]])$size))
  artifact$statistics_sha256 <- paste0("0", substring(artifact$statistics_sha256, 2))
  writeBin(secretbase::cborenc(artifact), bad_path)
  expect_error(assemble_ptm_cbor(statistics_path, c(bad_path, paths[-1]), output_path), "different statistics")
})
