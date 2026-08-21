test_that("site_column finds either spelling of the site identifier", {
  expect_equal(site_column(data.frame(site = "a")), "site")
  expect_equal(site_column(data.frame(protein_Id_site = "a")), "protein_Id_site")
})

test_that("site_column prefers `site` when a table carries both", {
  both <- data.frame(protein_Id_site = "a", site = "a")
  expect_equal(site_column(both), "site")
})

test_that("site_column names what it looked for when neither is present", {
  expect_error(
    site_column(data.frame(protein_Id = "P1")),
    "no site identifier column found"
  )
})

test_that("compute_dpa_dpu pairs every tested site with its protein", {
  dirs <- make_dea_pair()
  res <- suppressMessages(
    compute_dpa_dpu(dirs$phospho, dirs$protein)
  )

  # Two proteins carry a site, two contrasts each.
  expect_equal(nrow(res$combined_site_prot), 4)
  expect_true(all(c("diff.site", "diff.protein") %in% names(res$combined_site_prot)))

  # The site annotation of the phospho run reached the result.
  expect_true(all(
    c("posInProtein", "modAA", "SequenceWindow") %in%
      names(res$combined_site_prot)
  ))
})

test_that("compute_dpa_dpu reports the match rate per contrast", {
  dirs <- make_dea_pair()
  res <- suppressMessages(compute_dpa_dpu(dirs$phospho, dirs$protein))

  expect_equal(res$match_rates$contrast, c("a_vs_b", "c_vs_b"))
  expect_equal(res$match_rates$total_sites, c(2, 2))
  expect_equal(res$match_rates$matched_sites, c(2, 2))
  expect_equal(res$match_rates$match_rate, c(100, 100))
})

test_that("compute_dpa_dpu leaves an unmatched site without a protein estimate", {
  # A site on a protein the total-proteome run did not quantify.
  phospho <- make_dea_output(
    site_dea_table(protein_ids = c("P1", "P9")),
    site_annotation_table(c("P1", "P9")),
    name = "DEA_phospho_extra"
  )
  protein <- make_dea_output(
    protein_dea_table(protein_ids = "P1"),
    name = "DEA_protein_one"
  )

  res <- suppressMessages(compute_dpa_dpu(phospho, protein))

  unmatched <- res$combined_site_prot[
    res$combined_site_prot$protein_Id == "P9",
  ]
  expect_equal(nrow(unmatched), 2)
  expect_true(all(is.na(unmatched$diff.protein)))
  expect_equal(res$match_rates$matched_sites, c(1, 1))
})

test_that("compute_dpa_dpu computes the usage difference of a matched pair", {
  dirs <- make_dea_pair()
  res <- suppressMessages(compute_dpa_dpu(dirs$phospho, dirs$protein))

  paired <- res$combined_test_diff[res$combined_test_diff$measured_In == "both", ]
  expect_true(nrow(paired) > 0)
  expect_equal(paired$diff_diff, paired$diff.site - paired$diff.protein)
  expect_equal(paired$SE_I, sqrt(paired$std.error.site^2 + paired$std.error.protein^2))
})
