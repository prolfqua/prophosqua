paired_readers <- function() {
  paths <- anndata_pair_fixture()
  list(
    site = prolfquapp::DEAResultReader$new(paths$site),
    protein = prolfquapp::DEAResultReader$new(paths$protein)
  )
}

test_that("the DEA pair reader returns records keyed by the feature keys", {
  paths <- anndata_pair_fixture()
  pair <- .read_dea_pair(paths$dirs$phospho, paths$dirs$protein)

  expect_equal(pair$site$feature_keys, c("protein_Id", "site"))
  expect_equal(pair$protein$feature_keys, "protein_Id")
  expect_equal(nrow(pair$site$normalized_abundances), 96)
  expect_equal(nrow(pair$protein$normalized_abundances), 48)
  expect_equal(nrow(pair$site$differential_results), 16)
  expect_equal(nrow(pair$protein$differential_results), 8)
  expect_equal(nrow(pair$site$site_info), 16)
  expect_identical(pair$site$contrasts, c(a_vs_b = "G_a - G_b"))
  for (table in list(
    pair$site$normalized_abundances,
    pair$site$imputed_abundances,
    pair$site$imputation,
    pair$site$differential_results
  )) {
    expect_true(all(c("protein_Id", "site") %in% names(table)))
  }
  # The site annotation reaches the result rows by the join on the keys.
  expect_false(anyNA(pair$site$differential_results$SequenceWindow))
})

test_that("the DEA pair reader rejects swapped experiment roles", {
  paths <- anndata_pair_fixture()
  expect_error(.read_dea_pair(paths$dirs$protein, paths$dirs$phospho), "inputs may be swapped")
})

test_that("the DEA pair reader rejects missing hierarchy roles", {
  readers <- paired_readers()
  readers$site$subject_id <- "protein_Id"

  # Keyed by protein alone, the site results join their annotation many to
  # many before the role check stops the pairing.
  expect_error(suppressWarnings(.ptm_pair(readers$site, readers$protein)), "missing feature role.*site")
})

test_that("the DEA pair reader rejects mismatched samples", {
  readers <- paired_readers()
  readers$protein$samples <- readers$protein$samples[-1, ]
  expect_error(.ptm_pair(readers$site, readers$protein), "sample sets differ")
})

test_that("the DEA pair reader rejects inconsistent sample designs", {
  readers <- paired_readers()
  readers$protein$samples$G_[[1]] <- "different"
  expect_error(.ptm_pair(readers$site, readers$protein), "disagree on design factor 'G_'")
})

test_that("remove_contaminants drops the features the DEAs flag, from every table", {
  readers <- paired_readers()
  readers$site$annotation$CON <- readers$site$annotation$protein_Id == "P1"
  readers$protein$annotation$CON <- readers$protein$annotation$protein_Id == "P1"
  kept <- .ptm_pair(readers$site, readers$protein)
  removed <- .ptm_pair(readers$site, readers$protein, remove_contaminants = TRUE)

  expect_true("P1" %in% kept$site$differential_results$protein_Id)
  expect_true("P1" %in% kept$protein$differential_results$protein_Id)
  for (side in c("site", "protein")) {
    for (table in c("var", "normalized_abundances", "imputed_abundances", "imputation", "differential_results")) {
      expect_identical(
        nrow(removed[[side]][[table]]),
        sum(kept[[side]][[table]]$protein_Id != "P1"),
        label = paste(side, table)
      )
    }
  }
  expect_false("P1" %in% removed$site$site_info$protein_Id)
})

test_that("DPU rejects a DEA without the statistics it tests on", {
  paths <- anndata_pair_fixture()
  pair <- .read_dea_pair(paths$dirs$phospho, paths$dirs$protein)
  pair$site$differential_results$std.error.unmoderated <- NULL
  expect_error(.compute_dpa_dpu_from_pair(pair), "DPU site result is missing required column.*std.error.unmoderated")
})

test_that("compute_cf_dea fits the contrasts the DEAs recorded", {
  paths <- anndata_pair_fixture()
  result <- suppressWarnings(compute_cf_dea(paths$dirs$phospho, paths$dirs$protein))
  expect_identical(result$contrasts, c(a_vs_b = "G_a - G_b"))
  expect_setequal(unique(result$results$contrast), "a_vs_b")
})
