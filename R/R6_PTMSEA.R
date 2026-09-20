#' Complete PTMSEA stage
#'
#' Construct through `source$build(Type, analysis)` and persist with `write_h5mu(path)`.
#' @export
PTMSEA <- R6::R6Class(
  "PTMSEA",
  private = list(source = NULL, analysis = NULL, result = NULL),
  public = list(
    #' @description Construct a complete result or restore an existing one.
    #' @param source Required complete preceding stage.
    #' @param analysis DPA, DPU or CF.
    #' @param result Completed computation, supplied during restoration.
    initialize = function(source, analysis, result = .compute_ptmsea_stage(source, analysis)) {
      if (!inherits(source, "PTM_statistics")) {
        stop("PTMSEA requires PTM_statistics")
      }
      .validate_enrichment_analysis(analysis)
      .require_ptm_fields(
        result,
        c(
          "results",
          "ranks",
          "all_clean",
          "pathways",
          "data_info",
          "ptmsigdb_summary",
          "overlap_stats",
          "n_overlap",
          "n_our_sites",
          "prep_info",
          "results_info",
          "has_results",
          "analysis_inputs"
        ),
        "PTMSEA"
      )
      private$source <- source
      private$analysis <- analysis
      private$result <- .pack_ptm_value(result)
    },
    #' @description Return the complete preceding stage.
    get_source = function() private$source,
    #' @description Return the complete statistics component.
    get_statistics = function() private$source$get_statistics(),
    #' @description Return the analysis name.
    get_analysis = function() private$analysis,
    #' @description Return the complete result.
    get_results = function() .unpack_ptm_value(private$result),
    #' @description Build a complete next stage.
    #' @param Type Target R6 class.
    #' @param ... Target arguments.
    build = function(Type, ...) Type$new(self, ...),
    #' @description Return a detached storage representation.
    as_container = function() .enrichment_container(self),
    #' @description Write this complete stage atomically.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(self$as_container(), path)
  )
)
