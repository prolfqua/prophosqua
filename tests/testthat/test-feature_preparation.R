# Tests for feature preparation functions

test_that("filter_significant_sites filters by FDR and FC", {
  data <- data.frame(
    FDR.site = c(0.01, 0.03, 0.08, 0.02),
    diff.site = c(1.2, -0.8, 0.5, -1.5)
  )

  result <- filter_significant_sites(data, fdr_threshold = 0.05, fc_threshold = 0.6)

  # Should keep rows 1, 2, 4 (FDR < 0.05 and |FC| > 0.6)
  expect_equal(nrow(result), 3)
  expect_true(all(result$FDR.site < 0.05))
  expect_true(all(abs(result$diff.site) > 0.6))
})


test_that("filter_significant_sites adds regulation column", {
  data <- data.frame(
    FDR.site = c(0.01, 0.02, 0.03),
    diff.site = c(1.2, -0.8, -1.5)
  )

  result <- filter_significant_sites(data)

  expect_true("regulation" %in% colnames(result))
  expect_equal(result$regulation[result$diff.site > 0], "upregulated")
  expect_true(all(result$regulation[result$diff.site < 0] == "downregulated"))
})


test_that("filter_significant_sites works with custom column names", {
  data <- data.frame(
    my_fdr = c(0.01, 0.02, 0.08),
    my_fc = c(1.2, -0.8, 0.5)
  )

  result <- filter_significant_sites(
    data,
    fdr_col = "my_fdr",
    diff_col = "my_fc",
    fdr_threshold = 0.05,
    fc_threshold = 0.6
  )

  expect_equal(nrow(result), 2)
  expect_true("regulation" %in% colnames(result))
})


test_that("filter_significant_sites validates required columns", {
  data <- data.frame(
    FDR.site = c(0.01, 0.02),
    other_col = c(1.2, -0.8)
  )

  expect_error(
    filter_significant_sites(data),
    "Missing required columns"
  )
})


test_that("filter_significant_sites handles require_sequence", {
  data <- data.frame(
    FDR.site = c(0.01, 0.02, 0.03, 0.04),
    diff.site = c(1.2, -0.8, 1.0, -1.5),
    SequenceWindow = c("AAASAAAA", NA, "_BBBSBBB", "CCCSCCCC_")
  )

  # Without require_sequence
  result1 <- filter_significant_sites(data, require_sequence = FALSE)
  expect_equal(nrow(result1), 4)

  # With require_sequence - should filter out NA and underscore-bounded
 result2 <- filter_significant_sites(data, require_sequence = TRUE)
  expect_equal(nrow(result2), 1)
  expect_equal(result2$SequenceWindow, "AAASAAAA")
})


test_that("filter_significant_sites returns empty df when no sites pass", {
  data <- data.frame(
    FDR.site = c(0.1, 0.2, 0.3),
    diff.site = c(0.1, 0.2, 0.3)
  )

  result <- filter_significant_sites(data, fdr_threshold = 0.05, fc_threshold = 0.6)
  expect_equal(nrow(result), 0)
  expect_true("regulation" %in% colnames(result))
})


test_that("summarize_significant_sites counts by group", {
  data <- data.frame(
    contrast = c("A", "A", "A", "B", "B"),
    regulation = c("upregulated", "upregulated", "downregulated",
                   "downregulated", "downregulated")
  )

  result <- summarize_significant_sites(data)

  expect_equal(nrow(result), 2)
  expect_true("upregulated" %in% colnames(result))
  expect_true("downregulated" %in% colnames(result))

  a_row <- result[result$contrast == "A", ]
  expect_equal(a_row$upregulated, 2)
  expect_equal(a_row$downregulated, 1)
})


test_that("summarize_significant_sites validates regulation column", {
  data <- data.frame(
    contrast = c("A", "B"),
    other = c(1, 2)
  )

  expect_error(
    summarize_significant_sites(data),
    "regulation"
  )
})


test_that("summarize_significant_sites handles multiple group columns", {
  data <- data.frame(
    contrast = c("A", "A", "A", "A"),
    modAA = c("S", "S", "T", "T"),
    regulation = c("upregulated", "downregulated", "upregulated", "upregulated")
  )

  result <- summarize_significant_sites(data, group_cols = c("contrast", "modAA"))

  expect_equal(nrow(result), 2)  # A-S and A-T
})
