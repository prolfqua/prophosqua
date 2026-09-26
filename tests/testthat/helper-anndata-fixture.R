# The paired DEA artifacts the tests read: what prolfquapp's own DEA writes for
# the synthetic example pair.

make_anndata_pair <- function() {
  dirs <- example_dea_pair()
  list(
    site = get_dea_file(dirs$phospho, "AnnData.h5ad"),
    protein = get_dea_file(dirs$protein, "AnnData.h5ad"),
    dirs = dirs
  )
}

anndata_pair_fixture <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      value <<- make_anndata_pair()
    }
    value
  }
})

ptm_result_fixture <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      paths <- anndata_pair_fixture()
      input_hashes <- tools::md5sum(c(paths$site, paths$protein))
      output <- tempfile(fileext = ".h5mu")
      input <- tempfile(fileext = ".h5mu")
      import_ptm_h5mu(paths$site, paths$protein, input)
      suppressWarnings(compute_ptm_results_h5mu(input, output))
      value <<- c(
        paths,
        list(output = output, input_hashes = unname(input_hashes))
      )
    }
    value
  }
})

sort_ptm_result <- function(data) {
  keys <- intersect(c("protein_Id", "site", "contrast"), names(data))
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
