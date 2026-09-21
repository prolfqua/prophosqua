#' Export completed MuData as one terminal workbook
#'
#' Run only after reports have finished. No computation or report consumes the
#' workbook produced here. The final H5MU retains the complete enrichment JSON.
#' @param input_h5mu Final PTM_results artifact.
#' @param output_dir Delivery directory.
#' @return Workbook path, invisibly.
#' @export
export_ptm_h5mu <- function(input_h5mu, output_dir) {
  result <- read_ptm_h5mu(input_h5mu, PTM_results)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, "PTM_results.xlsx")
  writexl::write_xlsx(.ptm_workbook_tables(result), path)
  invisible(path)
}

.ptm_workbook_tables <- function(result) {
  tables <- result$get_tables()
  statistics <- result$get_statistics()
  cf <- statistics$get_cf()
  tables$CF_intensities <- cf$wide_data
  tables$CF_sample_annotation <- cf$wide_annotation |>
    dplyr::relocate(tidyselect::any_of("CONTROL"), .after = "G_")
  for (branch in result$get_enrichments()) {
    analysis <- branch$get_analysis()
    method <- class(branch)[1L]
    payload <- branch$get_results()
    table <- switch(method, PTMSEA = payload$all_clean, KinaseGSEA = payload$all_results, MEA = payload$mea_clean)
    tables[[paste(analysis, method, sep = "_")]] <- dplyr::select(
      table,
      -tidyselect::any_of(c("core_enrichment", "Leading.substrates"))
    )
    if (identical(method, "KinaseGSEA")) {
      tables[[paste(analysis, method, "summary", sep = "_")]] <- payload$gsea_info
    }
    if (identical(method, "MEA")) {
      tables[[paste(analysis, method, "summary", sep = "_")]] <- payload$summary_df
    }
  }
  tables
}
