test_that("final PTM results store all nine enrichments as JSON documents", {
  fixture <- ptm_enrichment_fixture()
  result <- fixture$final
  expected <- c(
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
  expect_identical(names(result$get_enrichment_documents()), expected)
  expect_identical(result$get_enrichment_document("PTMSEA", "dpa"), result$get_enrichment_documents()$PTMSEA__DPA)
  expect_error(result$get_enrichment_document("STRING", "DPA"), "Unknown enrichment method")
  expect_error(result$get_enrichment_document("PTMSEA", "unknown"), "Unknown PTM analysis")

  container <- result$as_container()
  expect_setequal(
    names(container$modalities$enriched$uns$prophosqua$enrichment_documents),
    expected[grepl("DPA$", expected)]
  )
  expect_setequal(
    names(container$modalities$cf$uns$prophosqua$enrichment_documents),
    expected[!grepl("DPA$", expected)]
  )
  expect_null(container$modalities$total$uns$prophosqua$enrichment_documents)
  for (modality in c("enriched", "cf")) {
    namespace <- container$modalities[[modality]]$uns$prophosqua
    expect_false(.contains_packed_gsea_result(namespace$completed_stages))
  }
  for (document in result$get_enrichment_documents()) {
    expect_named(document, c("format", "version", "json", "sha256"))
    expect_identical(document$version, "1.1.0")
    expect_false(.contains_packed_gsea_result(jsonlite::fromJSON(document$json, simplifyVector = FALSE)))
  }
})

test_that("final MuData reconstructs GSEA and MEA results from JSON", {
  skip_if_not_installed("enrichplot")
  fixture <- ptm_enrichment_fixture()
  original <- fixture$final
  restored <- read_ptm_h5mu(fixture$path, PTM_results)
  for (key in names(original$get_enrichment_documents())) {
    before_document <- original$get_enrichment_documents()[[key]]
    after_document <- restored$get_enrichment_documents()[[key]]
    for (field in names(before_document)) {
      expect_identical(after_document[[field]], before_document[[field]], info = paste(key, field))
    }
  }
  expect_identical(names(restored$get_enrichments()), names(original$get_enrichments()))

  gs_info <- utils::getFromNamespace("gsInfo", "enrichplot")
  for (method in c("PTMSEA", "KinaseGSEA")) {
    key <- paste0(method, "__DPA")
    field <- c(PTMSEA = "results", KinaseGSEA = "gsea_results")[[method]]
    before <- original$get_enrichments()[[key]]$get_results()[[field]]$a_vs_b
    after <- restored$get_enrichments()[[key]]$get_results()[[field]]$a_vs_b
    expect_equal(after@result, before@result)
    expect_equal(after@geneList, before@geneList)
    expect_equal(after@geneSets, before@geneSets)
    expect_equal(after@params, before@params)
    expect_equal(gs_info(after, geneSetID = 1), gs_info(before, geneSetID = 1))
    expect_silent(enrichplot::gseaplot2(after, geneSetID = 1))
  }

  before_mea <- original$get_enrichments()$MEA__DPA$get_results()
  after_mea <- restored$get_enrichments()$MEA__DPA$get_results()
  expect_equal(after_mea$mea_clean, before_mea$mea_clean, ignore_attr = TRUE)
  expect_equal(after_mea$summary_df, before_mea$summary_df, ignore_attr = TRUE)
  expect_equal(ptm_enrichment_report_data(fixture$path, "MEA", "DPA"), after_mea)
})

test_that("invalid or incomplete enrichment documents fail validation", {
  fixture <- ptm_enrichment_fixture()
  result <- fixture$final
  documents <- result$get_enrichment_documents()
  branches <- result$get_enrichments()
  statistics <- result$get_statistics()

  expect_error(
    PTM_results$new(statistics, branches, documents[-1]),
    "every enabled enrichment document"
  )
  unexpected <- documents
  unexpected$PTMSEA__OTHER <- unexpected$PTMSEA__DPA
  expect_error(
    PTM_results$new(statistics, branches, unexpected),
    "every enabled enrichment document"
  )
  corrupt <- documents
  corrupt$PTMSEA__DPA$sha256 <- paste0("0", substring(corrupt$PTMSEA__DPA$sha256, 2))
  expect_error(PTM_results$new(statistics, branches, corrupt), "checksum mismatch")
  malformed <- documents
  malformed$PTMSEA__DPA$json <- "{"
  malformed$PTMSEA__DPA$sha256 <- digest::digest(malformed$PTMSEA__DPA$json, algo = "sha256", serialize = FALSE)
  expect_error(PTM_results$new(statistics, branches, malformed), "Invalid enrichment JSON")
  wrong_version <- documents
  wrong_version$PTMSEA__DPA$version <- "1.0.0"
  expect_error(PTM_results$new(statistics, branches, wrong_version), "Unsupported enrichment document")
  extra_wrapper_field <- documents
  extra_wrapper_field$PTMSEA__DPA$extra <- "not allowed"
  expect_error(PTM_results$new(statistics, branches, extra_wrapper_field), "unexpected fields")

  missing_container <- result$as_container()
  missing_container$modalities$enriched$uns$prophosqua$enrichment_documents$PTMSEA__DPA <- NULL
  expect_error(.load_ptm_results(missing_container), "every enabled enrichment document")

  misplaced_container <- result$as_container()
  misplaced <- misplaced_container$modalities$enriched$uns$prophosqua$enrichment_documents$PTMSEA__DPA
  misplaced_container$modalities$enriched$uns$prophosqua$enrichment_documents$PTMSEA__DPA <- NULL
  misplaced_container$modalities$cf$uns$prophosqua$enrichment_documents$PTMSEA__DPA <- misplaced
  expect_error(.load_ptm_results(misplaced_container), "wrong modality")
})

test_that("enabled empty enrichments remain complete JSON documents", {
  fixture <- ptm_enrichment_fixture(empty = TRUE)
  restored <- read_ptm_h5mu(fixture$path, PTM_results)
  expect_length(restored$get_enrichment_documents(), 9L)
  expect_equal(nrow(restored$get_enrichments()$PTMSEA__DPA$get_results()$results$a_vs_b@result), 0L)
  expect_equal(nrow(restored$get_enrichments()$KinaseGSEA__DPA$get_results()$gsea_results$a_vs_b@result), 0L)
  expect_equal(nrow(restored$get_enrichments()$MEA__DPA$get_results()$mea_clean), 0L)
})
