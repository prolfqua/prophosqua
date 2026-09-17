#' Complete CorrectFirst analysis
#' @export
CF <- R6::R6Class(
  "CF",
  private = list(inputs = NULL, result = NULL),
  public = list(
    #' @description Compute or restore a complete CF result from paired inputs.
    #' @param inputs Complete paired DEA stage.
    #' @param result Completed computation, supplied when restoring a stage.
    initialize = function(inputs, result = .compute_cf_stage(inputs)) {
      stopifnot(inherits(inputs, "DEA_enriched_total"))
      .validate_cf_result(result)
      private$inputs <- inputs
      private$result <- .pack_cf_result(result)
    },
    #' @description Return the paired inputs.
    get_inputs = function() private$inputs,
    #' @description Return complete report and computation results.
    get_results = function() .unpack_cf_result(private$result),
    #' @description Build a complete next stage.
    #' @param Type Target R6 class.
    #' @param ... Target arguments.
    build = function(Type, ...) Type$new(self, ...),
    #' @description Return a detached storage representation.
    as_container = function() .cf_container(private$inputs, self$get_results()),
    #' @description Write this complete stage atomically.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(self$as_container(), path)
  )
)

#' Complete DPA and DPU analysis
#' @export
DPA_DPU <- R6::R6Class(
  "DPA_DPU",
  private = list(inputs = NULL, result = NULL),
  public = list(
    #' @description Compute or restore DPA and moderated/unmoderated DPU.
    #' @param inputs Complete paired DEA stage.
    #' @param result Completed computation, supplied when restoring a stage.
    initialize = function(inputs, result = .compute_dpa_dpu_from_pair(inputs$get_pair())) {
      stopifnot(inherits(inputs, "DEA_enriched_total"))
      .require_ptm_fields(
        result,
        c(
          "combined_site_prot",
          "combined_test_diff",
          "combined_test_diff_unmoderated",
          "n_unmoderated_untestable",
          "match_rates"
        ),
        "DPA/DPU"
      )
      private$inputs <- inputs
      private$result <- .pack_ptm_value(result)
    },
    #' @description Return the paired inputs.
    get_inputs = function() private$inputs,
    #' @description Return completed DPA/DPU results.
    get_results = function() .unpack_ptm_value(private$result),
    #' @description Build a complete next stage.
    #' @param Type Target R6 class.
    #' @param ... Target arguments.
    build = function(Type, ...) Type$new(self, ...),
    #' @description Return a detached storage representation.
    as_container = function() {
      container <- private$inputs$as_container()
      container$uns$prophosqua$stage <- "DPA_DPU"
      container$modalities$enriched$uns$prophosqua <- list(dpa_dpu = private$result)
      container
    },
    #' @description Write this complete stage atomically.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(self$as_container(), path)
  )
)

#' Complete PTM statistics composed from independent analyses
#' @export
PTM_statistics <- R6::R6Class(
  "PTM_statistics",
  private = list(dpa_dpu = NULL, cf = NULL),
  public = list(
    #' @description Combine completed DPA/DPU and CF analyses from the same inputs.
    #' @param dpa_dpu Complete DPA_DPU stage.
    #' @param cf Complete CF stage.
    initialize = function(dpa_dpu, cf) {
      stopifnot(inherits(dpa_dpu, "DPA_DPU"), inherits(cf, "CF"))
      .require_same_ptm_inputs(dpa_dpu$get_inputs(), cf$get_inputs())
      private$dpa_dpu <- dpa_dpu
      private$cf <- cf
    },
    #' @description Return the paired inputs.
    get_inputs = function() private$cf$get_inputs(),
    #' @description Return completed DPA/DPU report inputs.
    get_dpa_dpu = function() private$dpa_dpu$get_results(),
    #' @description Return completed CF report inputs.
    get_cf = function() private$cf$get_results(),
    #' @description Return the six standard delivery tables in memory.
    get_tables = function() .ptm_delivery_tables(self),
    #' @description Return this complete statistics component.
    get_statistics = function() self,
    #' @description Build a complete next stage.
    #' @param Type Target R6 class.
    #' @param ... Target arguments.
    build = function(Type, ...) Type$new(self, ...),
    #' @description Return a detached storage representation.
    as_container = function() .statistics_container(self),
    #' @description Write this complete stage atomically.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(self$as_container(), path)
  )
)

.compute_cf_stage <- function(inputs) {
  .compute_cf_dea_from_pair(inputs$get_pair(), inputs$get_design(), "MuData obs", inputs$get_contrasts())
}

.validate_cf_result <- function(result) {
  .require_ptm_fields(
    result,
    c(
      "results",
      "ptm_data",
      "ctr",
      "annot",
      "contrasts",
      "wide_data",
      "wide_annotation",
      "model_counts",
      "n_before",
      "n_protein_measurements",
      "n_site_measurements",
      "n_merged_measurements",
      "n_models",
      "n_site_contrast"
    ),
    "CF"
  )
}

.pack_cf_result <- function(result) {
  result$ptm_config <- prolfqua::R6_extract_values(result$ptm_data$get_config())
  result$ptm_long <- result$ptm_data$data_long()
  result$contrast_subject <- result$ctr$subject_id
  result$contrast_table <- result$ctr$get_contrasts()
  result$ptm_data <- NULL
  result$ctr <- NULL
  .pack_ptm_value(result)
}

.unpack_cf_result <- function(packed) {
  result <- .unpack_ptm_value(packed)
  result$ptm_data <- prolfqua::LFQData$new(result$ptm_long, prolfqua::list_to_AnalysisConfiguration(result$ptm_config))
  result$ctr <- prolfqua::ContrastsTable$new(result$contrast_table, subject_id = result$contrast_subject)
  result[c("ptm_config", "ptm_long", "contrast_subject", "contrast_table")] <- NULL
  .validate_cf_result(result)
  result
}

.require_same_ptm_inputs <- function(left, right) {
  getters <- c("get_provenance", "get_parameters", "get_resources", "get_contrasts", "get_design")
  same <- vapply(getters, function(getter) identical(left[[getter]](), right[[getter]]()), logical(1))
  if (!all(same)) {
    stop("PTM stages have different source inputs.")
  }
  for (getter in c("get_enriched", "get_total")) {
    a <- left[[getter]]()
    b <- right[[getter]]()
    if (!identical(a$var_names, b$var_names) || !isTRUE(all.equal(a$X, b$X))) {
      stop("PTM stages have different source measurements.")
    }
  }
}
