# Rscript tests/integration/compare_mudata_delivery.R baseline/PTM_results.rds new/PTM_results.rds
# Checks every cell, including identifiers, row order, columns and missingness.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
expected <- readRDS(args[[1L]])
actual <- readRDS(args[[2L]])
stopifnot(identical(names(expected), names(actual)))
for (sheet in names(expected)) {
  left <- expected[[sheet]]
  right <- actual[[sheet]]
  stopifnot(identical(names(left), names(right)), nrow(left) == nrow(right))
  for (column in names(left)) {
    x <- left[[column]]
    y <- right[[column]]
    stopifnot(identical(is.na(x), is.na(y)))
    keep <- !is.na(x)
    if (is.numeric(x) && is.numeric(y)) {
      equal <- x[keep] == y[keep] | abs(x[keep] - y[keep]) <= 1e-12 + 1e-8 * abs(x[keep])
      if (!all(equal)) stop(sheet, "/", column, ": numeric mismatch")
    } else if (!identical(as.character(x[keep]), as.character(y[keep]))) {
      stop(sheet, "/", column, ": value mismatch")
    }
  }
  cat(sheet, ":", nrow(left), "rows verified\n")
}
