ptm_varm_frame <- function(adata, key) {
  namespace <- adata$uns[["prophosqua"]]
  values <- as.matrix(adata$varm[[key]])
  colnames(values) <- as.character(namespace$varm_columns[[key]])
  as.data.frame(values)
}

ptm_result_row <- function(adata, key, site) {
  feature <- match(site, as.character(as.data.frame(adata$var)$site))
  expect_false(is.na(feature))
  ptm_varm_frame(adata, key)[feature, , drop = FALSE]
}

test_that("PTM result H5AD preserves the complete site input", {
  fixture <- ptm_result_fixture()
  site <- anndataR::read_h5ad(fixture$site)
  result <- anndataR::read_h5ad(fixture$output)

  expect_equal(
    unname(tools::md5sum(c(fixture$site, fixture$protein))),
    fixture$input_hashes
  )
  expect_equal(as.data.frame(result$obs), as.data.frame(site$obs))
  expect_equal(as.data.frame(result$var), as.data.frame(site$var))
  expect_equal(as.matrix(result$X), as.matrix(site$X))
  expect_setequal(result$layers_keys(), site$layers_keys())
  for (key in site$layers_keys()) {
    expect_equal(as.matrix(result$layers[[key]]), as.matrix(site$layers[[key]]), info = key)
  }
  for (key in site$varm_keys()) {
    expect_true(key %in% result$varm_keys(), info = key)
    expect_equal(as.matrix(result$varm[[key]]), as.matrix(site$varm[[key]]), info = key)
  }
  expect_null(site$uns[["prophosqua"]])
})

test_that("PTM result H5AD exposes a versioned, complete result contract", {
  fixture <- ptm_result_fixture()
  result <- anndataR::read_h5ad(fixture$output)
  namespace <- result$uns[["prophosqua"]]
  expected_keys <- c(
    "dpa__a_vs_b",
    "dpu__a_vs_b",
    "dpu_unmoderated__a_vs_b",
    "correct_first__a_vs_b"
  )

  expect_equal(namespace$artifact_type, "ptm_results")
  expect_equal(namespace$schema_version, "1.0.0")
  expect_equal(namespace$source_software, "prophosqua")
  expect_setequal(names(namespace$varm_columns), expected_keys)
  expect_setequal(names(namespace$varm_annotations), expected_keys)
  expect_setequal(names(namespace$varm_present), expected_keys)
  expect_setequal(unlist(namespace$result_keys, use.names = FALSE), expected_keys)
  expect_setequal(result$varm_keys(), c("dea__a_vs_b", expected_keys))

  expect_equal(namespace$site_input$path, normalizePath(fixture$site))
  expect_equal(namespace$protein_input$path, normalizePath(fixture$protein))
  expect_equal(namespace$site_input$hash, fixture$input_hashes[[1]])
  expect_equal(namespace$protein_input$hash, fixture$input_hashes[[2]])
  expect_equal(namespace$site_input$schema_version, "1.0.0")
  expect_equal(namespace$protein_input$schema_version, "1.0.0")
})

test_that("PTM result matrices retain statistics, annotations, and alignment", {
  fixture <- ptm_result_fixture()
  pair <- read_ptm_anndata_pair(fixture$site, fixture$protein)
  dpa_dpu <- suppressMessages(.compute_dpa_dpu_from_pair(pair))
  correct_first <- suppressWarnings(.compute_cf_dea_from_pair(
    pair,
    readr::read_tsv(fixture$annot_file, show_col_types = FALSE),
    basename(fixture$annot_file)
  ))
  result <- anndataR::read_h5ad(fixture$output)
  site <- "P1~S10"

  expected_dpa <- dpa_dpu$combined_site_prot[
    dpa_dpu$combined_site_prot$site == site,
    ,
    drop = FALSE
  ]
  actual_dpa <- ptm_result_row(result, "dpa__a_vs_b", site)
  expect_equal(actual_dpa$diff.site, expected_dpa$diff.site)
  expect_equal(actual_dpa$diff.protein, expected_dpa$diff.protein)

  expected_dpu <- dpa_dpu$combined_test_diff[
    dpa_dpu$combined_test_diff$site == site,
    ,
    drop = FALSE
  ]
  actual_dpu <- ptm_result_row(result, "dpu__a_vs_b", site)
  expect_equal(actual_dpu$diff_diff, expected_dpu$diff_diff)
  expect_equal(actual_dpu$FDR_I, expected_dpu$FDR_I)
  expect_equal(
    as.character(result$uns$prophosqua$varm_annotations$dpu__a_vs_b$measured_In)[
      match(site, as.data.frame(result$var)$site)
    ],
    expected_dpu$measured_In
  )

  expected_unmoderated <- dpa_dpu$combined_test_diff_unmoderated[
    dpa_dpu$combined_test_diff_unmoderated$site == site,
    ,
    drop = FALSE
  ]
  actual_unmoderated <- ptm_result_row(result, "dpu_unmoderated__a_vs_b", site)
  expect_equal(actual_unmoderated$diff_diff, expected_unmoderated$diff_diff)
  expect_equal(actual_unmoderated$FDR_I, expected_unmoderated$FDR_I)

  expected_cf <- correct_first$results[
    correct_first$results$site == site,
    ,
    drop = FALSE
  ]
  actual_cf <- ptm_result_row(result, "correct_first__a_vs_b", site)
  expect_equal(actual_cf$diff.site, expected_cf$diff.site)
  expect_equal(actual_cf$FDR.site, expected_cf$FDR.site)
  expect_true(result$uns$prophosqua$varm_present$dpa__a_vs_b[
    match(site, as.data.frame(result$var)$site)
  ])
})

test_that("PTM result alignment distinguishes absent rows from missing statistics", {
  fixture <- anndata_pair_fixture()
  pair <- read_ptm_anndata_pair(fixture$site, fixture$protein)
  dpa_dpu <- suppressMessages(.compute_dpa_dpu_from_pair(pair))
  omitted_site <- dpa_dpu$combined_site_prot$site[[1]]
  incomplete <- dpa_dpu$combined_site_prot[-1, , drop = FALSE]

  payload <- .ptm_analysis_payload(incomplete, pair$site$var, "dpa", "a_vs_b")
  feature <- match(omitted_site, pair$site$var$site)

  expect_false(payload$present$dpa__a_vs_b[[feature]])
  expect_true(all(is.na(payload$values$dpa__a_vs_b[feature, ])))
})

test_that("PTM result alignment rejects unknown and duplicate site identities", {
  fixture <- anndata_pair_fixture()
  pair <- read_ptm_anndata_pair(fixture$site, fixture$protein)
  dpa_dpu <- suppressMessages(.compute_dpa_dpu_from_pair(pair))

  unknown <- dpa_dpu$combined_site_prot
  unknown$site[[1]] <- "UNKNOWN~S1"
  expect_error(
    .ptm_analysis_payload(unknown, pair$site$var, "dpa", "a_vs_b"),
    "absent from the site AnnData axis"
  )

  duplicate <- rbind(
    dpa_dpu$combined_site_prot,
    dpa_dpu$combined_site_prot[1, , drop = FALSE]
  )
  expect_error(
    .ptm_analysis_payload(duplicate, pair$site$var, "dpa", "a_vs_b"),
    "duplicate rows"
  )
})

test_that("PTM result writer validates its destination before computation", {
  fixture <- anndata_pair_fixture()
  expect_error(
    compute_ptm_results_h5ad(
      fixture$site,
      fixture$protein,
      fixture$annot_file,
      fixture$site
    ),
    "must not overwrite either input"
  )
  expect_error(
    compute_ptm_results_h5ad(
      fixture$site,
      fixture$protein,
      fixture$annot_file,
      file.path(tempfile(), "PTM_results.h5ad")
    ),
    "output directory does not exist"
  )
  expect_error(
    compute_ptm_results_h5ad(
      fixture$site,
      fixture$protein,
      tempfile(fileext = ".tsv"),
      tempfile(fileext = ".h5ad")
    ),
    "Annotation file not found"
  )
})
