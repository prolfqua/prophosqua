# The two reports that read a finished DEA run must be renderable without one,
# which is what these example helpers are for. Testing them also gives
# compute_cf_dea() its only end-to-end exercise.

test_that("the example DEA pair has the shape prolfquapp writes", {
  dirs <- example_dea_pair()

  for (dea_dir in c(dirs$phospho, dirs$protein)) {
    expect_true(file.exists(get_dea_xlsx(dea_dir)))
    expect_true(file.exists(get_dea_parquet(dea_dir)))
    expect_true(file.exists(get_dea_yaml(dea_dir)))
  }
  expect_true(file.exists(dirs$annot_file))
})


test_that("the example phospho run carries the site annotation the compute step requires", {
  dirs <- example_dea_pair()
  de <- readxl::read_xlsx(get_dea_xlsx(dirs$phospho), sheet = "diff_exp_analysis")

  expect_true(all(c("posInProtein", "modAA", "SequenceWindow") %in% colnames(de)))
})


test_that("compute_dpa_dpu_example returns what the report reads", {
  objects <- compute_dpa_dpu_example()

  expect_named(
    objects,
    c(
      "match_rates",
      "n_dpa_rows",
      "n_dpu_rows",
      "phospho_dea_dir",
      "protein_dea_dir",
      "dpa_xlsx",
      "dpu_xlsx"
    )
  )
  expect_gt(objects$n_dpa_rows, 0)
  expect_gt(objects$n_dpu_rows, 0)
  expect_true(file.exists(objects$dpa_xlsx))
  expect_true(file.exists(objects$dpu_xlsx))
})


test_that("compute_cf_dea_example models the corrected usage", {
  cf <- suppressWarnings(compute_cf_dea_example())

  expect_gt(cf$n_models, 0)
  expect_gt(nrow(cf$results), 0)
  expect_true(all(c("diff.site", "FDR.site", "statistic.site") %in% colnames(cf$results)))
  # The two wide tables are dropped, as the pipeline drops them.
  expect_false(any(c("wide_data", "wide_annotation") %in% names(cf)))
})


test_that("the example site abundances are missing where the signal is faint", {
  dirs <- example_dea_pair()
  long <- arrow::read_parquet(get_dea_parquet(dirs$phospho))

  # An imputing aggregator fits a dropout model on this; a complete matrix
  # leaves it nothing to fit.
  expect_true(any(is.na(long$normalized_abundance)))
  expect_false(all(is.na(long$normalized_abundance)))
})
