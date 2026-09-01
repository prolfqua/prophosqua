sort_ptm_result <- function(data) {
  keys <- intersect(c("protein_Id", "site", "protein_Id_site", "contrast"), names(data))
  data <- as.data.frame(data)
  data <- data[do.call(order, data[keys]), , drop = FALSE]
  rownames(data) <- NULL
  data
}

expect_same_common_columns <- function(actual, expected, tolerance = 1e-10) {
  common <- intersect(names(expected), names(actual))
  expect_equal(
    sort_ptm_result(actual)[, common, drop = FALSE],
    sort_ptm_result(expected)[, common, drop = FALSE],
    tolerance = tolerance,
    ignore_attr = TRUE
  )
}

test_that("paired AnnData reader returns consumer-owned experiment records", {
  paths <- anndata_pair_fixture()
  pair <- read_ptm_anndata_pair(paths$site, paths$protein)

  expect_s3_class(pair, "prophosqua_anndata_pair")
  expect_s3_class(pair$site, "prophosqua_dea_experiment")
  expect_s3_class(pair$protein, "prophosqua_dea_experiment")
  expect_equal(pair$site$schema_version, "1.0.0")
  expect_equal(pair$site$feature_keys, c("protein_Id", "site"))
  expect_equal(pair$protein$feature_keys, "protein_Id")
  expect_equal(nrow(pair$site$normalized_abundances), 96)
  expect_equal(nrow(pair$protein$normalized_abundances), 48)
  expect_equal(nrow(pair$site$differential_results), 16)
  expect_equal(nrow(pair$protein$differential_results), 8)
  expect_equal(nrow(pair$site$site_info), 16)
})

test_that("paired AnnData reader aligns protein samples to the site order", {
  paths <- anndata_pair_fixture()
  reversed <- rewrite_fixture_h5ad(paths$protein, function(adata) {
    adata[rev(seq_len(adata$n_obs())), ]
  })

  pair <- read_ptm_anndata_pair(paths$site, reversed)

  expect_equal(
    pair$protein$obs[[pair$protein$sample_key]],
    pair$site$obs[[pair$site$sample_key]]
  )
  expect_equal(
    unique(pair$protein$normalized_abundances[[pair$protein$sample_key]]),
    pair$site$obs[[pair$site$sample_key]]
  )
})

test_that("paired AnnData reader rejects swapped experiment roles", {
  paths <- anndata_pair_fixture()
  expect_error(
    read_ptm_anndata_pair(paths$protein, paths$site),
    "inputs may be swapped"
  )
})

test_that("paired AnnData reader rejects unsupported schemas", {
  paths <- anndata_pair_fixture()
  unsupported <- rewrite_fixture_h5ad(paths$site, function(adata) {
    namespace <- adata$uns[["prolfquapp"]]
    namespace$schema_version <- "2.0.0"
    adata$uns[["prolfquapp"]] <- namespace
    adata
  })

  expect_error(
    read_ptm_anndata_pair(unsupported, paths$protein),
    "Unsupported prolfquapp DEA-results schema '2.0.0'"
  )
})

test_that("paired AnnData reader rejects missing hierarchy roles", {
  paths <- anndata_pair_fixture()
  missing_site <- rewrite_fixture_h5ad(paths$site, function(adata) {
    namespace <- adata$uns[["prolfquapp"]]
    namespace$feature_keys <- "protein_Id"
    adata$uns[["prolfquapp"]] <- namespace
    adata
  })

  expect_error(
    read_ptm_anndata_pair(missing_site, paths$protein),
    "missing feature role.*site"
  )
})

test_that("paired AnnData reader rejects mismatched and duplicate samples", {
  paths <- anndata_pair_fixture()
  missing_sample <- rewrite_fixture_h5ad(paths$protein, function(adata) {
    adata[seq.int(2L, adata$n_obs()), ]
  })
  expect_error(
    read_ptm_anndata_pair(paths$site, missing_sample),
    "sample sets differ"
  )

  duplicate_sample <- rewrite_fixture_h5ad(paths$protein, function(adata) {
    obs <- as.data.frame(adata$obs)
    obs[["Name"]][[2]] <- obs[["Name"]][[1]]
    adata$obs <- obs
    adata
  })
  expect_error(
    read_ptm_anndata_pair(paths$site, duplicate_sample),
    "sample-key names must be present, non-empty, and unique"
  )
})

test_that("paired AnnData reader rejects inconsistent sample designs", {
  paths <- anndata_pair_fixture()
  mismatched_design <- rewrite_fixture_h5ad(paths$protein, function(adata) {
    obs <- as.data.frame(adata$obs)
    obs[["G_"]][[1]] <- "different"
    adata$obs <- obs
    adata
  })

  expect_error(
    read_ptm_anndata_pair(paths$site, mismatched_design),
    "disagree on design factor 'G_'"
  )
})

test_that("paired AnnData reader rejects malformed files and missing layers", {
  paths <- anndata_pair_fixture()
  malformed <- tempfile(fileext = ".h5ad")
  writeLines("not an HDF5 file", malformed)
  expect_error(
    read_ptm_anndata_pair(malformed, paths$protein),
    "Cannot read prolfquapp AnnData file"
  )

  missing_raw <- rewrite_fixture_h5ad(paths$site, function(adata) {
    adata$layers[["raw"]] <- NULL
    adata
  })
  expect_error(
    read_ptm_anndata_pair(missing_raw, paths$protein),
    "missing required layer.*raw"
  )
})

test_that("AnnData DPA and DPU equal the legacy DEA-directory path", {
  paths <- anndata_pair_fixture()
  legacy <- suppressMessages(
    compute_dpa_dpu(paths$legacy$phospho, paths$legacy$protein)
  )
  h5ad <- compute_dpa_dpu_h5ad(paths$site, paths$protein)

  expect_same_common_columns(h5ad$combined_site_prot, legacy$combined_site_prot)
  expect_same_common_columns(h5ad$combined_test_diff, legacy$combined_test_diff)
  expect_same_common_columns(
    h5ad$combined_test_diff_unmoderated,
    legacy$combined_test_diff_unmoderated
  )
  expect_equal(h5ad$match_rates, legacy$match_rates)
  expect_equal(
    h5ad$n_unmoderated_untestable,
    legacy$n_unmoderated_untestable
  )
})

test_that("AnnData DPA and DPU reject missing required statistics", {
  paths <- anndata_pair_fixture()
  missing_statistic <- rewrite_fixture_h5ad(paths$site, function(adata) {
    namespace <- adata$uns[["prolfquapp"]]
    for (key in grep("^dea__", adata$varm_keys(), value = TRUE)) {
      columns <- namespace$varm_columns[[key]]
      keep <- columns != "std.error.unmoderated"
      adata$varm[[key]] <- as.matrix(adata$varm[[key]])[, keep, drop = FALSE]
      namespace$varm_columns[[key]] <- columns[keep]
    }
    adata$uns[["prolfquapp"]] <- namespace
    adata
  })

  expect_error(
    compute_dpa_dpu_h5ad(missing_statistic, paths$protein),
    "site DEA results is missing required column.*std.error.unmoderated"
  )
})

test_that("AnnData CorrectFirst equals the legacy DEA-directory path", {
  paths <- anndata_pair_fixture()
  legacy <- suppressWarnings(compute_cf_dea(
    paths$legacy$phospho,
    paths$legacy$protein,
    paths$annot_file
  ))
  h5ad <- suppressWarnings(compute_cf_dea_h5ad(
    paths$site,
    paths$protein,
    paths$annot_file
  ))

  expect_same_common_columns(h5ad$results, legacy$results)
  expect_equal(h5ad$model_counts, legacy$model_counts)
  expect_equal(h5ad$n_before, legacy$n_before)
  expect_equal(h5ad$n_protein_measurements, legacy$n_protein_measurements)
  expect_equal(h5ad$n_site_measurements, legacy$n_site_measurements)
  expect_equal(h5ad$n_merged_measurements, legacy$n_merged_measurements)
  expect_equal(h5ad$n_models, legacy$n_models)
  expect_equal(h5ad$n_site_contrast, legacy$n_site_contrast)
})
