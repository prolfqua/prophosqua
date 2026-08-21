test_that("canonicalize_sequence_window accepts the flanking column under either name", {
  already <- data.frame(SequenceWindow = "AAASAAA")
  expect_identical(canonicalize_sequence_window(already), already)

  other <- data.frame(PTM_FlankingRegion = "AAASAAA")
  expect_equal(
    suppressMessages(canonicalize_sequence_window(other))$SequenceWindow,
    "AAASAAA"
  )
})

test_that("canonicalize_sequence_window errors when neither name is present", {
  expect_error(
    canonicalize_sequence_window(data.frame(site = "P1~S1")),
    "No SequenceWindow or PTM_FlankingRegion"
  )
})

test_that("filter_sequence_windows drops windows a motif scan cannot use", {
  data <- data.frame(
    SequenceWindow = c(
      "AAAAAAASAAAAAAA", # usable
      "_AAAAAASAAAAAAA", # padded at the N terminus
      "AAAAAAASAAAAAA_", # padded at the C terminus
      "AASAAA", # too short
      "", # empty
      NA_character_
    )
  )
  out <- filter_sequence_windows(data)

  expect_equal(out$SequenceWindow, "AAAAAAASAAAAAAA")
})

test_that("filter_sequence_windows upper-cases what it keeps", {
  data <- data.frame(SequenceWindow = "aaaaaaasaaaaaaa")
  expect_equal(filter_sequence_windows(data)$SequenceWindow, "AAAAAAASAAAAAAA")
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

test_that("prep_kinaselib_inputs writes one seqwindow list and one rnk per contrast", {
  data <- data.frame(
    SequenceWindow = c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB", "AAAAAAASAAAAAAA"),
    contrast = c("a_vs_b", "a_vs_b", "c/b"),
    statistic.site = c(2, -3, 1)
  )
  xlsx <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(DPA = data), xlsx)
  out_dir <- tempfile("KinaseLib")

  written <- suppressMessages(prep_kinaselib_inputs(
    xlsx,
    out_dir,
    analysis_type = "DPA",
    sheet = "DPA",
    stat_column = "statistic.site"
  ))

  expect_true(file.exists(file.path(out_dir, "DPA_seqwindows.tsv")))
  # The contrast name carries a slash, which cannot go into a file name.
  expect_true(file.exists(file.path(out_dir, "DPA_MEA_c_b.rnk")))
  expect_true(file.exists(file.path(out_dir, "DPA_MEA_a_vs_b.rnk")))
  expect_length(written, 3)

  seqwindows <- utils::read.delim(file.path(out_dir, "DPA_seqwindows.tsv"))
  expect_equal(seqwindows$SequenceWindow, c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB"))
})

test_that("prep_kinaselib_inputs names the available columns when the statistic is absent", {
  xlsx <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(
    list(DPA = data.frame(SequenceWindow = "AAAAAAASAAAAAAA", contrast = "a_vs_b")),
    xlsx
  )
  expect_error(
    suppressMessages(prep_kinaselib_inputs(
      xlsx,
      tempfile(),
      "DPA",
      "DPA",
      "statistic.site"
    )),
    "Statistic column 'statistic.site' not found"
  )
})
