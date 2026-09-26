test_that("DPA pairs every tested site with its protein", {
  res <- suppressMessages(.compute_dpa_dpu_from_pair(dea_result_pair()))

  # Two proteins carry a site, two contrasts each.
  expect_equal(nrow(res$combined_site_prot), 4)
  expect_true(all(c("diff.site", "diff.protein") %in% names(res$combined_site_prot)))

  # The site annotation of the phospho run reached the result.
  expect_true(all(
    c("posInProtein", "modAA", "SequenceWindow") %in%
      names(res$combined_site_prot)
  ))
})

test_that("DPA reports the match rate per contrast", {
  res <- suppressMessages(.compute_dpa_dpu_from_pair(dea_result_pair()))

  expect_equal(res$match_rates$contrast, c("a_vs_b", "c_vs_b"))
  expect_equal(res$match_rates$total_sites, c(2, 2))
  expect_equal(res$match_rates$matched_sites, c(2, 2))
  expect_equal(res$match_rates$match_rate, c(100, 100))
})

test_that("DPA leaves an unmatched site without a protein estimate", {
  # A site on a protein the total-proteome run did not quantify.
  pair <- dea_result_pair(
    site_dea_table(protein_ids = c("P1", "P9")),
    protein_dea_table(protein_ids = "P1")
  )

  res <- suppressMessages(.compute_dpa_dpu_from_pair(pair))

  unmatched <- res$combined_site_prot[
    res$combined_site_prot$protein_Id == "P9",
  ]
  expect_equal(nrow(unmatched), 2)
  expect_true(all(is.na(unmatched$diff.protein)))
  expect_equal(res$match_rates$matched_sites, c(1, 1))
})

test_that("DPU is the usage difference of a matched pair", {
  res <- suppressMessages(.compute_dpa_dpu_from_pair(dea_result_pair()))

  paired <- res$combined_test_diff[res$combined_test_diff$measured_In == "both", ]
  expect_true(nrow(paired) > 0)
  expect_equal(paired$diff_diff, paired$diff.site - paired$diff.protein)
  expect_equal(paired$SE_I, sqrt(paired$std.error.site^2 + paired$std.error.protein^2))
  expect_equal(
    res$combined_test_diff_unmoderated$SE_I[
      res$combined_test_diff_unmoderated$measured_In == "both"
    ],
    sqrt(paired$std.error.unmoderated.site^2 + paired$std.error.unmoderated.protein^2)
  )
  expect_equal(res$n_unmoderated_untestable, 0)
})

test_that("DPU counts paired rows with invalid raw degrees of freedom", {
  site <- site_dea_table(protein_ids = "P1", contrasts = c("a_vs_b", "c_vs_b"))
  protein <- protein_dea_table(protein_ids = "P1", contrasts = c("a_vs_b", "c_vs_b"))
  protein$df.unmoderated[protein$contrast == "a_vs_b"] <- 0

  res <- suppressMessages(.compute_dpa_dpu_from_pair(dea_result_pair(site, protein)))
  raw <- res$combined_test_diff_unmoderated

  expect_equal(res$n_unmoderated_untestable, 1)
  expect_true(is.na(raw$pValue_I[raw$contrast == "a_vs_b"]))
  expect_false(is.na(raw$pValue_I[raw$contrast == "c_vs_b"]))
})
