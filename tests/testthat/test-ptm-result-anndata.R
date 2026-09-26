ptm_result_row <- function(adata, key, site) {
  frame <- adata$varm[[key]]
  frame[frame$site %in% site, , drop = FALSE]
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
    expect_equal(result$varm[[key]], site$varm[[key]], info = key)
  }
  expect_null(site$uns[["prophosqua"]])
  total <- anndataR::read_h5ad(fixture$protein)
  expect_equal(container$modalities$total$X, total$X)
  expect_equal(container$modalities$total$var, total$var)
  expect_equal(container$uns$prophosqua$schema_version, "2.0.0")
  expect_equal(container$uns$prophosqua$stage, "PTM_statistics")
})

test_that("PTM result tables retain statistics, annotations, and alignment", {
  fixture <- ptm_result_fixture()
  pair <- .read_dea_pair(fixture$dirs$phospho, fixture$dirs$protein)
  dpa_dpu <- suppressMessages(.compute_dpa_dpu_from_pair(pair))
  correct_first <- suppressWarnings(.compute_cf_dea_from_pair(pair))
  container <- prolfquapp::read_h5mu(fixture$output)
  enriched <- container$modalities$enriched
  site <- "P1~S10"
  expected <- function(table) table[table$site == site, , drop = FALSE]

  actual_dpa <- ptm_result_row(enriched, "dpa__a_vs_b", site)
  expect_equal(nrow(actual_dpa), 1L)
  expect_equal(actual_dpa$diff.site, expected(dpa_dpu$combined_site_prot)$diff.site)
  expect_equal(actual_dpa$diff.protein, expected(dpa_dpu$combined_site_prot)$diff.protein)

  actual_dpu <- ptm_result_row(enriched, "dpu__a_vs_b", site)
  expect_equal(actual_dpu$diff_diff, expected(dpa_dpu$combined_test_diff)$diff_diff)
  expect_equal(actual_dpu$FDR_I, expected(dpa_dpu$combined_test_diff)$FDR_I)
  expect_equal(as.character(actual_dpu$measured_In), expected(dpa_dpu$combined_test_diff)$measured_In)

  actual_unmoderated <- ptm_result_row(enriched, "dpu_unmoderated__a_vs_b", site)
  expect_equal(actual_unmoderated$diff_diff, expected(dpa_dpu$combined_test_diff_unmoderated)$diff_diff)
  expect_equal(actual_unmoderated$FDR_I, expected(dpa_dpu$combined_test_diff_unmoderated)$FDR_I)

  actual_cf <- ptm_result_row(container$modalities$enriched_CF, "correct_first__a_vs_b", site)
  expect_equal(actual_cf$diff.site, expected(correct_first$results)$diff.site)
  expect_equal(actual_cf$FDR.site, expected(correct_first$results)$FDR.site)
})

test_that("a PTM result frame has a row for every site, keyed only where a result exists", {
  fixture <- anndata_pair_fixture()
  pair <- .read_dea_pair(fixture$dirs$phospho, fixture$dirs$protein)
  dpa_dpu <- suppressMessages(.compute_dpa_dpu_from_pair(pair))
  omitted_site <- dpa_dpu$combined_site_prot$site[[1]]
  incomplete <- dpa_dpu$combined_site_prot[-1, , drop = FALSE]

  frame <- .ptm_analysis_payload(incomplete, pair$site$var, "dpa", "a_vs_b")$values$dpa__a_vs_b

  expect_equal(nrow(frame), nrow(pair$site$var))
  expect_false(omitted_site %in% frame$site)
  expect_setequal(frame$site[!is.na(frame$site)], incomplete$site)
})

test_that("PTM result alignment rejects duplicate site identities", {
  fixture <- anndata_pair_fixture()
  pair <- .read_dea_pair(fixture$dirs$phospho, fixture$dirs$protein)
  dpa_dpu <- suppressMessages(.compute_dpa_dpu_from_pair(pair))

  duplicate <- rbind(
    dpa_dpu$combined_site_prot,
    dpa_dpu$combined_site_prot[1, , drop = FALSE]
  )
  expect_error(
    .ptm_analysis_payload(duplicate, pair$site$var, "dpa", "a_vs_b"),
    "duplicate rows"
  )
})

test_that("contrast encoding keeps a present row whose statistics are missing", {
  fixture <- anndata_pair_fixture()
  pair <- .read_dea_pair(fixture$dirs$phospho, fixture$dirs$protein)
  data <- .compute_dpa_dpu_from_pair(pair)$combined_site_prot
  first <- data[1, , drop = FALSE]
  first$contrast <- "a/b"
  numeric <- vapply(first, is.numeric, logical(1))
  first[, numeric] <- NA_real_
  second <- data[2, , drop = FALSE]
  second$contrast <- "a%2Fb"

  frames <- .ptm_analysis_payload(rbind(first, second), pair$site$var, "dpa", c("a/b", "a%2Fb"))$values

  expect_equal(names(frames), c("dpa__a%2Fb", "dpa__a%252Fb"))
  row <- frames[[1]][frames[[1]]$site %in% first$site, , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_true(is.na(row$diff.site))
  expect_false(first$site %in% frames[[2]]$site)
})
