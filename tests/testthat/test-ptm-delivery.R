test_that("terminal statistics exports retain the legacy metadata column order", {
  paths <- anndata_pair_fixture()
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(paths$site),
    anndataR::read_h5ad(paths$protein),
    parameters = list(
      analyses = list(dpa = list(subdir = "DPA"), dpu = list(subdir = "DPU"), cf = list(subdir = "CF"))
    )
  )
  statistics <- inputs$build(DPA_DPU)$build(PTM_statistics, cf = suppressWarnings(inputs$build(CF)))
  before <- statistics$get_dpa_dpu()
  output <- tempfile()
  .export_ptm_statistics(statistics, output)
  tables <- list(
    readxl::read_xlsx(file.path(output, "DPA", "Result_DPA.xlsx")),
    readRDS(file.path(output, "DPU", "combined_test_diff.rds"))
  )
  for (table in tables) {
    expect_identical(
      intersect(names(table), c("estimate_type.site", "contrast", "diff.site")),
      c("estimate_type.site", "contrast", "diff.site")
    )
    expect_lt(match("estimate_type.protein", names(table)), match("diff.protein", names(table)))
  }
  expect_equal(statistics$get_dpa_dpu(), before)
})

test_that("terminal DPA/DPU annotations retain legacy blank-cell missingness", {
  source <- data.frame(
    gene_name.site = c("", "GENE1", NA_character_),
    diff.site = c(1, 2, NA_real_),
    gene_name.protein = c("GENE2", "", NA_character_),
    diff.protein = c(3, NA_real_, 4)
  )
  output <- .ptm_dpa_dpu_delivery(source)
  expect_identical(output$gene_name.site, c(NA_character_, "GENE1", NA_character_))
  expect_identical(output$gene_name.protein, c("GENE2", NA_character_, NA_character_))
  expect_identical(output$diff.site, source$diff.site)
  expect_identical(output$diff.protein, source$diff.protein)
  expect_identical(source$gene_name.site[[1L]], "")
})
