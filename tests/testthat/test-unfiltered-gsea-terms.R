test_that("PTM-SEA and Kinase GSEA retain tested terms above the FDR cutoff", {
  set.seed(481)
  amino_acids <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1L]]
  sequences <- replicate(180L, paste0(c(sample(amino_acids, 7L), "S", sample(amino_acids, 7L)), collapse = ""))
  expect_length(unique(sequences), length(sequences))
  data <- data.frame(
    SequenceWindow = sequences,
    statistic.site = seq(1, -1, length.out = length(sequences)),
    contrast = "A_vs_B"
  )
  sets <- setNames(replicate(6L, sample(sequences, 35L), simplify = FALSE), paste0("set_", seq_len(6L)))

  ptmsea <- suppressWarnings(.compute_ptmsea_tables(
    data,
    lapply(sets, paste0, "-p"),
    "DPA",
    "statistic.site",
    trim_to = 15L,
    min_size = 10L,
    max_size = 500L,
    n_perm = 1000L
  ))
  expect_setequal(ptmsea$all_clean$ID, names(sets))
  expect_true(any(ptmsea$all_clean$p.adjust > 0.25))

  term2gene <- data.frame(
    term = rep(names(sets), lengths(sets)),
    gene = unlist(sets, use.names = FALSE)
  )
  kinase <- suppressWarnings(.compute_kinase_tables(
    data,
    term2gene,
    "DPA",
    min_size = 10L,
    max_size = 5000L,
    n_perm = 1000L
  ))
  expect_setequal(kinase$all_results$kinase, names(sets))
  expect_true(any(kinase$all_results$FDR > 0.25))
})
