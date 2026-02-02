test_that("prepare_enrichment_data adds computed columns", {
  df <- data.frame(NES = c(1, -1, 0.5), FDR = c(0.01, 0.2, 0.05))
  result <- prepare_enrichment_data(df)

  expect_true("neg_log_fdr" %in% names(result))
  expect_true("direction" %in% names(result))
  expect_true("significant" %in% names(result))
  expect_equal(result$direction, c("Up", "Down", "Up"))
  expect_equal(result$significant, c(TRUE, FALSE, TRUE))
})

test_that("prepare_enrichment_data respects custom fdr_col", {
  df <- data.frame(NES = c(1, -1), p.adjust = c(0.01, 0.2))
  result <- prepare_enrichment_data(df, fdr_col = "p.adjust")

  expect_equal(result$significant, c(TRUE, FALSE))
})

test_that("prepare_enrichment_data respects custom fdr_threshold", {
  df <- data.frame(NES = c(1, -1), FDR = c(0.04, 0.06))
  result <- prepare_enrichment_data(df, fdr_threshold = 0.05)

  expect_equal(result$significant, c(TRUE, FALSE))
})

test_that("plot_enrichment_dotplot returns ggplot", {
  df <- data.frame(
    kinase = c("A", "B", "C"),
    NES = c(1.5, -1.2, 0.8),
    FDR = c(0.01, 0.05, 0.2)
  )
  p <- plot_enrichment_dotplot(df, n_top = 3)
  expect_s3_class(p, "ggplot")
})

test_that("plot_enrichment_volcano returns ggplot", {
  df <- data.frame(
    kinase = c("A", "B", "C"),
    NES = c(1.5, -1.2, 0.8),
    FDR = c(0.01, 0.05, 0.2),
    contrast = c("X", "X", "X")
  )
  p <- plot_enrichment_volcano(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_enrichment_heatmap returns ggplot", {
  df <- data.frame(
    ID = c("P1", "P2", "P1", "P2"),
    NES = c(1.5, -1.2, 0.8, -0.5),
    p.adjust = c(0.01, 0.05, 0.2, 0.08),
    contrast = c("A", "A", "B", "B")
  )
  p <- plot_enrichment_heatmap(df, n_top = 2)
  expect_s3_class(p, "ggplot")
})

test_that("summarize_enrichment_results handles data frame input", {
  df <- data.frame(
    contrast = c("A", "A", "A", "B", "B"),
    FDR = c(0.01, 0.08, 0.15, 0.03, 0.12)
  )
  result <- summarize_enrichment_results(df)

  expect_equal(nrow(result), 2)
  expect_true("total" %in% names(result))
  expect_true("FDR < 0.1" %in% names(result))
  expect_true("FDR < 0.05" %in% names(result))
})

test_that("extract_gsea_results handles empty results gracefully", {
  result <- extract_gsea_results(setNames(list(), character(0)))
  expect_equal(nrow(result), 0)
  expect_true("ID" %in% names(result))
  expect_true("NES" %in% names(result))
  expect_true("contrast" %in% names(result))
})
