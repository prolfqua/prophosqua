test_that(".rank_windows ranks each contrast's windows in descending order", {
  data <- data.frame(
    contrast = rep(c("A_vs_B", "C_vs_D"), each = 3),
    SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCSCCCC", "DDDSDDDD", "EEESEEEE", "FFFSFFF"),
    statistic.site = c(2.5, 1.0, -0.5, 1.8, -1.2, 0.3)
  )
  ranks <- .rank_windows(data, "statistic.site")
  expect_named(ranks, c("A_vs_B", "C_vs_D"))
  expect_true(all(vapply(ranks, function(x) all(diff(x) <= 0), logical(1))))
  expect_identical(names(ranks$A_vs_B), c("AAASAAAA", "BBBSBBB", "CCCSCCCC"))
})

test_that(".rank_windows keeps the first statistic of a repeated window", {
  data <- data.frame(
    contrast = "A_vs_B",
    SequenceWindow = c("AAASAAAA", "AAASAAAA", "BBBSBBB", "BBBSBBB"),
    statistic.site = c(2.5, 1.0, -0.5, 0.8)
  )
  expect_equal(.rank_windows(data, "statistic.site")$A_vs_B, c(AAASAAAA = 2.5, BBBSBBB = -0.5))
})

test_that(".rank_windows drops sites without a statistic or a window", {
  data <- data.frame(
    contrast = "A_vs_B",
    SequenceWindow = c("AAASAAAA", "BBBSBBB", NA),
    statistic.site = c(2.5, NA, 1)
  )
  expect_equal(.rank_windows(data, "statistic.site")$A_vs_B, c(AAASAAAA = 2.5))
})

test_that(".rank_windows trims windows to the PTMsigDB width", {
  data <- data.frame(contrast = "A_vs_B", SequenceWindow = "ABCDEFGSIJKLMNO", statistic.site = 1)
  expect_named(.rank_windows(data, "statistic.site", trim_to = 11L)$A_vs_B, "CDEFGSIJKLM")
})

test_that(".rank_windows names the missing columns", {
  expect_error(.rank_windows(data.frame(contrast = "a"), "statistic.site"), "SequenceWindow, statistic.site")
})

test_that("PTMsigDB sites match ranked windows once their suffixes are dropped", {
  pathways <- list(set = c("AAAAAAASAAAAAAA-p;u", "AAAAAAASAAAAAAA-p;d", "BBBBBBBSBBBBBBB-p"))
  expect_equal(.ptmsigdb_windows(pathways), list(set = c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB")))
  expect_equal(
    trim_ptmsigdb_pathways(pathways, "11"),
    list(set = c("AAAAASAAAAA-p;u", "AAAAASAAAAA-p;d", "BBBBBSBBBBB-p"))
  )
})
