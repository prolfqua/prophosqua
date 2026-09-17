test_that("paired inputs build independent complete analyses without mutation", {
  paths <- anndata_pair_fixture()
  inputs <- DEA_enriched_total$new(anndataR::read_h5ad(paths$site), anndataR::read_h5ad(paths$protein))
  before <- inputs$get_enriched()
  expect_false("get_cf" %in% names(inputs))
  expect_error(CF$new(inputs, result = list()), "incomplete")
  dpa_dpu <- inputs$build(DPA_DPU)
  cf <- suppressWarnings(inputs$build(CF))
  statistics <- dpa_dpu$build(PTM_statistics, cf = cf)
  expect_true(inherits(statistics, "PTM_statistics"))
  expect_equal(inputs$get_enriched()$X, before$X)
  expect_equal(inputs$get_enriched()$uns, before$uns)
  expect_equal(inputs$get_pair()$site$configuration$get_response(), "normalized_abundance")
  expect_equal(cf$get_results()$contrasts, inputs$get_contrasts())
  expect_same_common_columns(
    statistics$get_dpa_dpu()$combined_test_diff,
    compute_dpa_dpu(paths$legacy$phospho, paths$legacy$protein)$combined_test_diff
  )
})

test_that("complete statistics round-trip through MuData without refitting", {
  paths <- anndata_pair_fixture()
  input <- tempfile(fileext = ".h5mu")
  output <- tempfile(fileext = ".h5mu")
  import_ptm_h5mu(paths$site, paths$protein, input)
  original <- suppressWarnings(compute_ptm_results_h5mu(input, output))
  restored <- read_ptm_h5mu(output, PTM_statistics)
  expect_equal(restored$get_tables(), original$get_tables())
  expect_equal(restored$get_cf()$ctr$get_contrasts(), original$get_cf()$ctr$get_contrasts())
  expect_error(read_ptm_h5mu(input, PTM_statistics), "Expected stage")
  container <- prolfquapp::read_h5mu(output)
  expect_equal(names(container$modalities), c("enriched", "total", "cf"))
  expect_true(all(c("dpa__a_vs_b", "dpa__a_vs_b__present") %in% container$modalities$enriched$varm_keys()))
  expect_true(all(c("dpu__a_vs_b", "correct_first__a_vs_b") %in% container$modalities$cf$varm_keys()))
  mask <- container$modalities$cf$varm$dpu__a_vs_b__present
  container$modalities$cf$varm$dpu__a_vs_b__present <- NULL
  prolfquapp::write_h5mu(container$modalities, output, container$obs, container$uns)
  expect_error(read_ptm_h5mu(output), "PTM result matrices is incomplete")
  container$modalities$cf$varm$dpu__a_vs_b__present <- mask
  container$modalities$cf$uns$prophosqua$report_data <- NULL
  prolfquapp::write_h5mu(container$modalities, output, container$obs, container$uns)
  expect_error(read_ptm_h5mu(output), "CF metadata is incomplete")
})

test_that("portable report values preserve empty results, named ranks and missingness", {
  value <- list(
    empty = data.frame(term = character(), p = double()),
    ranks = c(site_a = 1.3, site_b = -0.2),
    annotation = c("a", NA_character_),
    shape = matrix(c(1, NA, 3, 4), 2, dimnames = list(c("a", "b"), c("x", "y")))
  )
  expect_equal(.copy_ptm_value(value), value)
})
