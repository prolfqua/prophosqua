#' Complete PTM statistics: DPA, DPU and CorrectFirst
#'
#' Computed from the paired DEA experiments, or restored from stored results.
#' @export
PTM_statistics <- R6::R6Class(
  "PTM_statistics",
  private = list(inputs = NULL, dpa_dpu = NULL, cf = NULL),
  public = list(
    #' @description Compute DPA, moderated and unmoderated DPU, and CorrectFirst
    #'   with its imputed variants, or restore them.
    #' @param inputs Complete DEA_enriched_total stage.
    #' @param dpa_dpu Stored DPA/DPU result, when restoring.
    #' @param cf Stored CorrectFirst result, when restoring.
    initialize = function(
      inputs,
      dpa_dpu = .compute_dpa_dpu_from_pair(inputs$get_pair()),
      cf = .compute_cf_dea_from_pair(inputs$get_pair())
    ) {
      stopifnot(inherits(inputs, "DEA_enriched_total"))
      .require_ptm_fields(dpa_dpu, c(.DPA_DPU_METHODS, .DPA_DPU_SUMMARIES), "DPA/DPU")
      .require_ptm_fields(cf, c("results", "ptm_data", "wide_data", "wide_annotation", "variants", .CF_METADATA), "CF")
      private$inputs <- inputs
      private$dpa_dpu <- dpa_dpu
      private$cf <- cf
    },
    #' @description Return the paired inputs.
    get_inputs = function() private$inputs,
    #' @description Return the DPA/DPU results.
    get_dpa_dpu = function() private$dpa_dpu,
    #' @description Return the CorrectFirst results; the imputed variants are
    #'   under `variants`.
    get_cf = function() private$cf,
    #' @description Return the six standard delivery tables in memory.
    get_tables = function() .ptm_delivery_tables(self),
    #' @description Return this complete statistics component.
    get_statistics = function() self,
    #' @description Return a detached storage representation.
    as_container = function() .statistics_container(self),
    #' @description Write this complete stage atomically.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(self$as_container(), path)
  )
)

# The DPA/DPU tables, by the method name each is stored under.
.DPA_DPU_METHODS <- c(
  dpa = "combined_site_prot",
  dpu = "combined_test_diff",
  dpu_unmoderated = "combined_test_diff_unmoderated"
)

.DPA_DPU_SUMMARIES <- c("n_unmoderated_untestable", "match_rates")

# What CorrectFirst stores as metadata in uns; the rest of its result is
# rebuilt from the enriched_CF modality.
.CF_METADATA <- c(
  "contrasts",
  "model_counts",
  "n_before",
  "n_protein_measurements",
  "n_site_measurements",
  "n_merged_measurements",
  "n_models",
  "n_site_contrast"
)
