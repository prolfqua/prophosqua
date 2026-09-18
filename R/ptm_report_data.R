#' Read DPA/DPU report inputs from final MuData
#' @param path Final PTM_results artifact.
#' @return Complete list used by the integration report.
#' @export
ptm_dpa_dpu_report_data <- function(path) {
  statistics <- read_ptm_h5mu(path, PTM_results)$get_statistics()
  result <- statistics$get_dpa_dpu()
  provenance <- statistics$get_inputs()$get_provenance()
  list(
    match_rates = result$match_rates,
    n_dpa_rows = nrow(result$combined_site_prot),
    n_dpu_rows = nrow(result$combined_test_diff),
    phospho_dea_dir = dirname(provenance$paths[[1L]]),
    protein_dea_dir = dirname(provenance$paths[[2L]]),
    dpa_xlsx = "Result_DPA.xlsx",
    dpu_xlsx = "Result_DPU.xlsx"
  )
}

#' Read a completed enrichment report from final MuData
#' @param path Final PTM_results artifact.
#' @param method PTMSEA, KinaseGSEA or MEA.
#' @param analysis DPA, DPU or CF.
#' @return Complete report result; absent analyses raise an error.
#' @export
ptm_enrichment_report_data <- function(path, method, analysis) {
  result <- read_ptm_h5mu(path, PTM_results)
  key <- .ptm_varm_key(method, analysis)
  document <- result$get_enrichment_document(method, analysis)
  .restore_ptm_enrichment_result(document, key) # nolint: object_usage_linter.
}
