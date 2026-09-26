# The example DEA pair every other test reads, made by prolfquapp's own DEA.

example_reader <- function(dea_dir) {
  prolfquapp::DEAResultReader$new(get_dea_file(dea_dir, "AnnData.h5ad"))
}

test_that("the example DEA pair holds the artifact prolfquapp writes", {
  dirs <- example_dea_pair()

  for (dea_dir in c(dirs$phospho, dirs$protein)) {
    expect_true(file.exists(get_dea_file(dea_dir, "AnnData.h5ad")))
  }
})

test_that("the example phospho run carries the site annotation the compute step requires", {
  dirs <- example_dea_pair()
  annotation <- example_reader(dirs$phospho)$annotation

  expect_true(all(c("posInProtein", "modAA", "SequenceWindow") %in% colnames(annotation)))
})

test_that("the example site abundances are missing where the signal is faint", {
  dirs <- example_dea_pair()
  lfq <- example_reader(dirs$phospho)$lfq_transformed
  values <- lfq$data_long()[[lfq$response()]]

  # An imputing model fits its dropout on this; a complete matrix leaves it
  # nothing to fit.
  expect_true(any(is.na(values)))
  expect_false(all(is.na(values)))
})
