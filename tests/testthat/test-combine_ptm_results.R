test_that("standardize_ptm_results renames the DPU statistics to the shared names", {
  dpu <- data.frame(
    protein_Id = "P1",
    site = "P1~S1",
    contrast = "a_vs_b",
    gene_name.site = "GENE",
    diff_diff = 1.5,
    FDR_I = 0.01,
    tstatistic_I = 4
  )
  out <- standardize_ptm_results(dpu, "dpu")

  expect_equal(out$diff.site, 1.5)
  expect_equal(out$FDR.site, 0.01)
  expect_equal(out$statistic.site, 4)
  expect_equal(out$gene_name, "GENE")
  expect_false(any(c("diff_diff", "FDR_I", "tstatistic_I") %in% names(out)))
})

test_that("standardize_ptm_results keeps the DPA statistics of both levels", {
  dpa <- data.frame(
    protein_Id = "P1",
    site = "P1~S1",
    contrast = "a_vs_b",
    gene_name.site = "GENE",
    diff.site = 1,
    FDR.site = 0.01,
    statistic.site = 3,
    diff.protein = 0.5,
    FDR.protein = 0.02,
    statistic.protein = 2
  )
  out <- standardize_ptm_results(dpa, "dpa")

  expect_true(all(c("diff.site", "diff.protein") %in% names(out)))
  expect_equal(out$gene_name, "GENE")
})

test_that("standardize_ptm_results suffixes the single CF estimate type", {
  cf <- data.frame(
    protein_Id = "P1",
    site = "P1~S1",
    contrast = "a_vs_b",
    gene_name = "GENE",
    diff.site = 1,
    FDR.site = 0.01,
    statistic.site = 3,
    estimate_type = "observed"
  )
  out <- standardize_ptm_results(cf, "cf")

  expect_equal(out$estimate_type.site, "observed")
  expect_false("estimate_type" %in% names(out))
})

test_that("standardize_ptm_results is case-insensitive in the analysis type", {
  dpu <- data.frame(protein_Id = "P1", diff_diff = 1)
  expect_identical(
    standardize_ptm_results(dpu, "DPU"),
    standardize_ptm_results(dpu, "dpu")
  )
})

test_that("standardize_ptm_results drops a column the analysis did not produce", {
  # A DEA run without flanking sequences still yields a usable table; the motif
  # reports are what go missing, and this silence is why.
  without_window <- data.frame(
    protein_Id = "P1",
    site = "P1~S1",
    contrast = "a_vs_b",
    diff.site = 1,
    FDR.site = 0.01,
    statistic.site = 3
  )
  out <- standardize_ptm_results(without_window, "cf")

  expect_false("SequenceWindow" %in% names(out))
  expect_equal(nrow(out), 1)
})

test_that("standardize_ptm_results names the types it accepts", {
  expect_error(
    standardize_ptm_results(data.frame(protein_Id = "P1"), "ptmsea"),
    "Must be one of: dpa, dpu, cf"
  )
})

test_that("standardize_ptm_results keeps the column order each sheet is read in", {
  # The three sheets are read by name downstream, but a reordering still shows
  # up in every workbook a person opens, so the order is pinned here.
  full <- data.frame(
    protein_Id = "P1",
    site = "P1~S1",
    contrast = "a_vs_b",
    posInProtein = 1L,
    modAA = "S",
    SequenceWindow = "AAASAAA",
    protein_length = 300L,
    gene_name = "GENE",
    gene_name.site = "GENE",
    diff.site = 1,
    FDR.site = 0.01,
    statistic.site = 3,
    diff.protein = 1,
    FDR.protein = 0.01,
    statistic.protein = 3,
    diff_diff = 1,
    FDR_I = 0.01,
    tstatistic_I = 3,
    estimate_type = "observed",
    estimate_type.site = "observed",
    estimate_type.protein = "observed"
  )

  expect_named(
    standardize_ptm_results(full, "dpa"),
    c(
      "protein_Id",
      "site",
      "contrast",
      "posInProtein",
      "modAA",
      "SequenceWindow",
      "protein_length",
      "diff.site",
      "FDR.site",
      "statistic.site",
      "diff.protein",
      "FDR.protein",
      "statistic.protein",
      "estimate_type.site",
      "estimate_type.protein",
      "gene_name"
    )
  )
  expect_named(
    standardize_ptm_results(full, "dpu"),
    c(
      "protein_Id",
      "site",
      "contrast",
      "posInProtein",
      "modAA",
      "SequenceWindow",
      "protein_length",
      "estimate_type.site",
      "estimate_type.protein",
      "gene_name",
      "diff.site",
      "FDR.site",
      "statistic.site"
    )
  )
  expect_named(
    standardize_ptm_results(full, "cf"),
    c(
      "protein_Id",
      "site",
      "contrast",
      "posInProtein",
      "modAA",
      "SequenceWindow",
      "gene_name",
      "protein_length",
      "diff.site",
      "FDR.site",
      "statistic.site",
      "estimate_type.site"
    )
  )
})
