test_that("kinase GSEA retains substrate sets larger than the PTM-SEA limit", {
  ids <- sprintf("SITE%04d", seq_len(1500L))
  data <- data.frame(
    SequenceWindow = ids,
    statistic.site = seq(15, -15, length.out = length(ids)),
    contrast = "A_vs_B"
  )
  term2gene <- data.frame(
    term = rep(c("top", "bottom"), each = 700L),
    gene = c(ids[c(1:600, 701:800)], ids[c(701:800, 901:1500)])
  )

  result <- suppressWarnings(.compute_kinase_tables(
    data,
    term2gene,
    "DPA",
    min_size = 10L,
    max_size = 5000L,
    n_perm = 1000L
  ))

  expect_setequal(result$all_results$kinase, c("top", "bottom"))
  expect_true(all(result$all_results$setSize == 700L))
  expect_true(all(result$all_results$FDR < 0.25))
})
