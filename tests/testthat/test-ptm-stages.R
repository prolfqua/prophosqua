test_that("paired inputs compute complete statistics without mutation", {
  paths <- anndata_pair_fixture()
  inputs <- DEA_enriched_total$new(anndataR::read_h5ad(paths$site), anndataR::read_h5ad(paths$protein))
  before <- inputs$get_enriched()
  expect_false("get_cf" %in% names(inputs))
  expect_error(PTM_statistics$new(inputs, cf = list()), "CF is incomplete")
  statistics <- suppressWarnings(PTM_statistics$new(inputs))
  expect_equal(inputs$get_enriched()$X, before$X)
  expect_equal(inputs$get_enriched()$uns, before$uns)
  expect_equal(inputs$get_pair()$site$configuration$get_response(), "normalized_abundance")
  expect_equal(statistics$get_cf()$contrasts, inputs$get_contrasts())
  expect_same_common_columns(
    statistics$get_dpa_dpu()$combined_test_diff,
    compute_dpa_dpu(paths$dirs$phospho, paths$dirs$protein)$combined_test_diff
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
  sheet <- suppressMessages(original$get_tables())$abundances_site_cf
  cf_wide <- original$get_cf()$wide_data
  expect_equal(
    as.matrix(sheet[, -1]),
    as.matrix(cf_wide[match(sheet$site, cf_wide$site), names(sheet)[-1]]),
    ignore_attr = TRUE
  )
  expect_error(read_ptm_h5mu(input, PTM_statistics), "Expected stage")
  container <- prolfquapp::read_h5mu(output)
  expect_equal(names(container$modalities), c("enriched", "total", "enriched_CF"))
  expect_true("dpa__a_vs_b" %in% container$modalities$enriched$varm_keys())
  expect_true("correct_first__a_vs_b" %in% container$modalities$enriched_CF$varm_keys())
  expect_true("dpu__a_vs_b" %in% container$modalities$enriched$varm_keys())
  dpu <- container$modalities$enriched$varm$dpu__a_vs_b
  container$modalities$enriched$varm$dpu__a_vs_b <- NULL
  prolfquapp::write_h5mu(container$modalities, output, container$obs, container$uns)
  expect_error(read_ptm_h5mu(output), "PTM result table is missing")
  container$modalities$enriched$varm$dpu__a_vs_b <- dpu
  container$modalities$enriched_CF$uns$prophosqua$cf <- NULL
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
  expect_equal(.unpack_ptm_value(.pack_ptm_value(value)), value)
})

test_that("the peptide-level total DEA travels through every stage as a passive modality", {
  paths <- anndata_pair_fixture()
  input <- tempfile(fileext = ".h5mu")
  output <- tempfile(fileext = ".h5mu")
  # The protein DEA stands in for a peptide-level one: same samples, same contrasts.
  import_ptm_h5mu(paths$site, paths$protein, input, total_peptide_h5ad = paths$protein)
  inputs <- read_ptm_h5mu(input, DEA_enriched_total)
  expect_equal(names(inputs$get_provenance()$paths), c("enriched", "total", "total_peptide"))
  suppressWarnings(compute_ptm_results_h5mu(input, output))
  container <- prolfquapp::read_h5mu(output)
  expect_equal(
    names(container$modalities),
    c("enriched", "total", "total_peptide", "enriched_CF")
  )
  expect_equal(
    container$modalities$total_peptide$X,
    anndataR::read_h5ad(paths$protein)$X,
    ignore_attr = TRUE
  )
})
