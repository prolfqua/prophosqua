#' Complete MotifEnrichment stage
#'
#' Construct through `source$build(Type, analysis)` and persist with `write_h5mu(path)`.
#' @export
MotifEnrichment <- R6::R6Class(
  "MotifEnrichment",
  private = list(source = NULL, analysis = NULL, result = NULL),
  public = list(
    #' @description Construct a complete result or restore an existing one.
    #' @param source Required complete preceding stage.
    #' @param analysis DPA, DPU or CF.
    #' @param result Completed computation, supplied during restoration.
    initialize = function(source, analysis, result = .compute_motifenrichment_stage(source, analysis)) {
      if (!inherits(source, "KinaseAssignments")) {
        stop("MotifEnrichment requires KinaseAssignments")
      }
      .validate_enrichment_analysis(analysis)
      if (!identical(analysis, source$get_analysis())) {
        stop("Source analysis differs from target analysis.")
      }
      .require_ptm_fields(result, "mea_results", "MotifEnrichment")
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
