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

test_that("PTM result MuData preserves both complete inputs", {
  fixture <- ptm_result_fixture()
  site <- anndataR::read_h5ad(fixture$site)
  container <- prolfquapp::read_h5mu(fixture$output)
  result <- container$modalities$enriched

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
  total <- anndataR::read_h5ad(fixture$protein)
  expect_equal(container$modalities$total$X, total$X)
  expect_equal(container$modalities$total$var, total$var)
  expect_equal(container$uns$prophosqua$schema_version, "2.0.0")
  expect_equal(container$uns$prophosqua$stage, "PTM_statistics")
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
  container <- prolfquapp::read_h5mu(fixture$output)
  result <- container$modalities$enriched
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
  actual_dpu <- ptm_result_row(container$modalities$cf, "dpu__a_vs_b", site)
  expect_equal(actual_dpu$diff_diff, expected_dpu$diff_diff)
  expect_equal(actual_dpu$FDR_I, expected_dpu$FDR_I)
  expect_equal(
    as.character(container$modalities$cf$uns$prophosqua$varm_annotations$dpu__a_vs_b$measured_In)[
      match(site, as.data.frame(container$modalities$cf$var)$site)
    ],
    expected_dpu$measured_In
  )

  expected_unmoderated <- dpa_dpu$combined_test_diff_unmoderated[
    dpa_dpu$combined_test_diff_unmoderated$site == site,
    ,
    drop = FALSE
  ]
  actual_unmoderated <- ptm_result_row(container$modalities$cf, "dpu_unmoderated__a_vs_b", site)
  expect_equal(actual_unmoderated$diff_diff, expected_unmoderated$diff_diff)
  expect_equal(actual_unmoderated$FDR_I, expected_unmoderated$FDR_I)

  expected_cf <- correct_first$results[
    correct_first$results$site == site,
    ,
    drop = FALSE
  ]
  actual_cf <- ptm_result_row(container$modalities$cf, "correct_first__a_vs_b", site)
  expect_equal(actual_cf$diff.site, expected_cf$diff.site)
  expect_equal(actual_cf$FDR.site, expected_cf$FDR.site)
  expect_true(as.vector(result$varm$dpa__a_vs_b__present)[
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

test_that("contrast encoding and presence distinguish missing statistics from absent rows", {
  fixture <- anndata_pair_fixture()
  pair <- read_ptm_anndata_pair(fixture$site, fixture$protein)
  data <- .compute_dpa_dpu_from_pair(pair)$combined_site_prot
  first <- data[1, , drop = FALSE]
  first$contrast <- "a/b"
  numeric <- vapply(first, is.numeric, logical(1))
  first[, numeric] <- NA_real_
  second <- data[2, , drop = FALSE]
  second$contrast <- "a%2Fb"
  payload <- .ptm_analysis_payload(rbind(first, second), pair$site$var, "dpa", c("a/b", "a%2Fb"))
  expect_equal(unname(payload$keys), c("dpa__a%2Fb", "dpa__a%252Fb"))
  feature <- match(first$site, pair$site$var$site)
  expect_true(payload$present[[1]][feature])
  expect_false(payload$present[[2]][feature])
  expect_true(all(is.na(payload$values[[1]][feature, ])))
})
