# Tests for prepare_gsea_ranks()

test_that("prepare_gsea_ranks creates sorted named vectors", {
  data <- data.frame(
    contrast = rep(c("A_vs_B", "C_vs_D"), each = 3),
    SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCSCCCC",
                       "DDDSDDDD", "EEESEEEE", "FFFSFFF"),
    statistic.site = c(2.5, 1.0, -0.5, 1.8, -1.2, 0.3),
    stringsAsFactors = FALSE
  )

  ranks <- prepare_gsea_ranks(data, stat_col = "statistic.site")

  # Should return a list with two contrasts

  expect_type(ranks, "list")
  expect_equal(length(ranks), 2)
  expect_equal(names(ranks), c("A_vs_B", "C_vs_D"))

  # Each element should be a named numeric vector
  expect_type(ranks[[1]], "double")
  expect_true(!is.null(names(ranks[[1]])))

  # Should be sorted descending
  expect_true(all(diff(ranks[[1]]) <= 0))
  expect_true(all(diff(ranks[[2]]) <= 0))
})


test_that("prepare_gsea_ranks handles duplicates by keeping first", {
  data <- data.frame(
    contrast = rep("A_vs_B", 4),
    SequenceWindow = c("AAASAAAA", "AAASAAAA", "BBBSBBB", "BBBSBBB"),
    statistic.site = c(2.5, 1.0, -0.5, 0.8),
    stringsAsFactors = FALSE
  )

  ranks <- prepare_gsea_ranks(data, stat_col = "statistic.site")

  # Should have only 2 unique sequences
  expect_equal(length(ranks[["A_vs_B"]]), 2)

  # Should keep first occurrence (2.5 for AAASAAAA, -0.5 for BBBSBBB)
  expect_equal(unname(ranks[["A_vs_B"]]["AAASAAAA"]), 2.5)
  expect_equal(unname(ranks[["A_vs_B"]]["BBBSBBB"]), -0.5)
})


test_that("prepare_gsea_ranks adds suffix correctly", {
  data <- data.frame(
    contrast = "A_vs_B",
    SequenceWindow = c("AAASAAAA", "BBBSBBB"),
    statistic.site = c(2.5, 1.0),
    stringsAsFactors = FALSE
  )

  ranks <- prepare_gsea_ranks(data, stat_col = "statistic.site", add_suffix = "-p")

  expect_true(all(grepl("-p$", names(ranks[[1]]))))
  expect_equal(names(ranks[[1]])[1], "AAASAAAA-p")
})


test_that("prepare_gsea_ranks converts to uppercase", {
  data <- data.frame(
    contrast = "A_vs_B",
    SequenceWindow = c("aaasaaaa", "BbbSbbb"),
    statistic.site = c(2.5, 1.0),
    stringsAsFactors = FALSE
  )

  ranks <- prepare_gsea_ranks(data, stat_col = "statistic.site")
  expect_equal(names(ranks[[1]]), c("AAASAAAA", "BBBSBBB"))

  ranks_no_upper <- prepare_gsea_ranks(data, stat_col = "statistic.site", to_uppercase = FALSE)
  expect_equal(names(ranks_no_upper[[1]]), c("aaasaaaa", "BbbSbbb"))
})


test_that("prepare_gsea_ranks removes NA values", {
  data <- data.frame(
    contrast = "A_vs_B",
    SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCSCCCC", NA),
    statistic.site = c(2.5, NA, 1.0, 0.5),
    stringsAsFactors = FALSE
  )

  ranks <- prepare_gsea_ranks(data, stat_col = "statistic.site")

  # Should only have 2 valid entries (AAASAAAA and CCCSCCCC)
  expect_equal(length(ranks[[1]]), 2)
  expect_false("BBBSBBB" %in% names(ranks[[1]]))
})


test_that("prepare_gsea_ranks validates required columns", {
  data <- data.frame(
    contrast = "A_vs_B",
    SequenceWindow = "AAASAAAA",
    stringsAsFactors = FALSE
  )

  expect_error(
    prepare_gsea_ranks(data, stat_col = "nonexistent"),
    "Missing required columns"
  )
})


test_that("prepare_gsea_ranks works with custom column names", {
  data <- data.frame(
    my_contrast = "A_vs_B",
    my_seq = c("AAASAAAA", "BBBSBBB"),
    my_stat = c(2.5, 1.0),
    stringsAsFactors = FALSE
  )

  ranks <- prepare_gsea_ranks(
    data,
    stat_col = "my_stat",
    seq_col = "my_seq",
    contrast_col = "my_contrast"
  )

  expect_equal(names(ranks), "A_vs_B")
  expect_equal(length(ranks[[1]]), 2)
})
