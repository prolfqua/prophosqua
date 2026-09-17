#' Export completed MuData as terminal delivery files
#'
#' Run only after reports have finished. No computation or report consumes the
#' workbooks or RDS files produced here.
#' @param input_h5mu Final PTM_results artifact.
#' @param output_dir Delivery directory.
#' @return Paths written, invisibly.
#' @export
export_ptm_h5mu <- function(input_h5mu, output_dir) {
  result <- read_ptm_h5mu(input_h5mu, PTM_results)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- result$get_tables()
  .write_ptm_delivery(tables, file.path(output_dir, "PTM_results"))
  .export_ptm_statistics(result$get_statistics(), output_dir)
  exporters <- list(PTMSEA = .export_ptmsea, KinaseGSEA = .export_kinase_gsea, MEA = .export_mea)
  for (branch in result$get_enrichments()) {
    exporters[[class(branch)[1L]]](branch, output_dir)
  }
  invisible(list.files(output_dir, pattern = "\\.(xlsx|rds)$", recursive = TRUE, full.names = TRUE))
}

.write_ptm_delivery <- function(tables, stem, object = tables) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(tables, paste0(stem, ".xlsx"))
  saveRDS(object, paste0(stem, ".rds"))
}

.export_ptm_statistics <- function(statistics, output_dir) {
  parameters <- statistics$get_inputs()$get_parameters()
  directories <- lapply(parameters$analyses, function(value) file.path(output_dir, value$subdir))
  dpa_dpu <- statistics$get_dpa_dpu()
  dpa_dpu$combined_site_prot <- .ptm_dpa_dpu_delivery(dpa_dpu$combined_site_prot)
  dpa_dpu$combined_test_diff <- .ptm_dpa_dpu_delivery(dpa_dpu$combined_test_diff)
  cf <- statistics$get_cf()
  cf$wide_annotation <- cf$wide_annotation |>
    dplyr::relocate(tidyselect::any_of("CONTROL"), .after = "G_")
  .write_ptm_delivery(
    list(combinedSiteProteinData = dpa_dpu$combined_site_prot),
    file.path(directories$dpa, "Result_DPA")
  )
  .write_ptm_delivery(list(combinedStats = dpa_dpu$combined_test_diff), file.path(directories$dpu, "Result_DPU"))
  .write_ptm_delivery(list(results = cf$results), file.path(directories$cf, "CorrectFirst_PTM_usage_results"))
  saveRDS(dpa_dpu$combined_test_diff, file.path(directories$dpu, "combined_test_diff.rds"))
  writexl::write_xlsx(cf$wide_data, file.path(directories$cf, "CorrectFirst_intensities_wide.xlsx"))
  writexl::write_xlsx(cf$wide_annotation, file.path(directories$cf, "CorrectFirst_intensities_file_annotation.xlsx"))
}

.ptm_dpa_dpu_delivery <- function(table) {
  table |>
    dplyr::relocate(tidyselect::any_of(c("modelName.site", "estimate_type.site", "contrast")), .before = "diff.site") |>
    dplyr::relocate(tidyselect::any_of(c("modelName.protein", "estimate_type.protein")), .before = "diff.protein")
}

.enrichment_delivery_stem <- function(branch, output_dir, subdir, name) {
  analysis <- branch$get_analysis()
  settings <- branch$get_statistics()$get_inputs()$get_parameters()$analyses[[tolower(analysis)]]
  file.path(output_dir, settings$subdir, subdir, name)
}

.export_ptmsea <- function(branch, output_dir) {
  result <- branch$get_results()
  data <- result$all_clean |> dplyr::arrange(.data$contrast, .data$pvalue)
  sheets <- list(all_clean = data)
  for (contrast in unique(data$contrast)) {
    key <- gsub("[^a-zA-Z0-9_]", "_", substr(contrast, 1, 31))
    sheets[[key]] <- data[data$contrast == contrast, , drop = FALSE]
  }
  sheets$significant_FDR10 <- data |> dplyr::filter(.data$p.adjust < 0.1)
  stem <- .enrichment_delivery_stem(branch, output_dir, "PTMSEA", paste0("PTMSEA_", branch$get_analysis(), "_results"))
  .write_ptm_delivery(sheets, stem, result)
}

.export_kinase_gsea <- function(branch, output_dir) {
  result <- branch$get_results()
  sheets <- list(
    all_results = result$all_results |> dplyr::arrange(.data$contrast, .data$FDR),
    significant = result$all_results |> dplyr::filter(.data$FDR < 0.1) |> dplyr::arrange(.data$contrast, .data$FDR),
    summary = result$gsea_info
  )
  stem <- .enrichment_delivery_stem(branch, output_dir, "KinaseLib", paste0("KinaseLib_GSEA_", branch$get_analysis()))
  .write_ptm_delivery(sheets, stem, result)
}

.export_mea <- function(branch, output_dir) {
  result <- branch$get_results()
  columns <- c("contrast", "kinase", "NES", "pvalue", "FDR", "n_leading", "set_size")
  data <- result$mea_clean |> dplyr::select(tidyselect::all_of(columns)) |> dplyr::arrange(.data$contrast, .data$FDR)
  sheets <- list(all_results = data, significant = data |> dplyr::filter(.data$FDR < 0.1), summary = result$summary_df)
  stem <- .enrichment_delivery_stem(branch, output_dir, "KinaseLib", paste0("MEA_", branch$get_analysis(), "_results"))
  .write_ptm_delivery(sheets, stem, result)
}
