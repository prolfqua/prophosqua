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
