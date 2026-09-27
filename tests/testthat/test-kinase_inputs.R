test_that("filter_sequence_windows drops windows a motif scan cannot use", {
  data <- data.frame(
    SequenceWindow = c(
      "AAAAAAASAAAAAAA", # usable
      "XXXAAAASAAAAAAA", # padded at the N terminus
      "AAAAAAASAAAAXXX", # padded at the C terminus
      NA_character_ # no window: the protein is not in the FASTA
    )
  )
  out <- filter_sequence_windows(data)

  expect_equal(out$SequenceWindow, "AAAAAAASAAAAAAA")
})

test_that("rank_sites_for_mea selects one contrast and orders it descending", {
  data <- data.frame(
    SequenceWindow = c("W1", "W2", "W3"),
    contrast = c("a_vs_b", "a_vs_b", "c_vs_b"),
    statistic.site = c(1, 3, 9)
  )
  out <- rank_sites_for_mea(data, "statistic.site", "a_vs_b")

  expect_equal(out$SequenceWindow, c("W2", "W1"))
  expect_equal(out$statistic.site, c(3, 1))
})

test_that("rank_sites_for_mea keeps the most extreme statistic of a repeated window", {
  # Two sites share a flanking sequence and disagree. Averaging them would
  # cancel to nearly nothing; the enrichment should see the stronger signal.
  data <- data.frame(
    SequenceWindow = c("W1", "W1"),
    contrast = "a_vs_b",
    statistic.site = c(1.5, -4)
  )
  out <- rank_sites_for_mea(data, "statistic.site", "a_vs_b")

  expect_equal(nrow(out), 1)
  expect_equal(out$statistic.site, -4)
})

test_that("rank_sites_for_mea drops sites without a statistic", {
  data <- data.frame(
    SequenceWindow = c("W1", "W2"),
    contrast = "a_vs_b",
    statistic.site = c(2, NA)
  )
  expect_equal(rank_sites_for_mea(data, "statistic.site", "a_vs_b")$SequenceWindow, "W1")
})

test_that("rank_sites_for_mea ranks on the requested column", {
  data <- data.frame(
    SequenceWindow = c("W1", "W2"),
    contrast = "a_vs_b",
    statistic.site = c(1, 2),
    other_stat = c(9, 1)
  )
  out <- rank_sites_for_mea(data, "other_stat", "a_vs_b")
  expect_equal(out$SequenceWindow, c("W1", "W2"))
  expect_equal(out$statistic.site, c(9, 1))
})
