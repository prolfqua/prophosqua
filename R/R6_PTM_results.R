#' Final PTM results: the statistics and the enrichment of every enabled analysis
#'
#' The enrichment is kept in files beside the final MuData, laid out as the
#' pipeline lays out each analysis directory: protsea documents for the
#' completed PTM-SEA, Kinase GSEA and MEA, gzipped JSON for the kinase-library
#' preparations. The MuData holds only their names and checksums; the results
#' are read from the files when first asked for.
#' @export
PTM_results <- R6::R6Class(
  "PTM_results",
  private = list(
    statistics = NULL,
    files = NULL,
    statistics_sha256 = NULL,
    enrichments = NULL,
    container = function(path) {
      container <- private$statistics$as_container()
      namespace <- container$uns$prophosqua
      namespace$stage <- "PTM_results"
      if (length(private$files)) {
        namespace$enrichment_files <- as.list(.ptm_relative_paths(private$files, path))
        namespace$enrichment_sha256 <- as.list(vapply(private$files, .ptm_file_sha256, character(1)))
        namespace$statistics_sha256 <- private$statistics_sha256
      }
      container$uns$prophosqua <- namespace
      container
    }
  ),
  public = list(
    #' @description Pair the statistics with the enrichment files of every enabled analysis.
    #' @param statistics Complete PTM_statistics stage.
    #' @param files Every stage file of every enabled analysis, keyed by stage
    #'   and analysis, as the pipeline lays them out.
    #' @param statistics_sha256 Checksum of the statistics file the kinase
    #'   preparations were computed from.
    initialize = function(statistics, files = character(), statistics_sha256 = NULL) {
      stopifnot(inherits(statistics, "PTM_statistics"))
      expected <- names(.ptm_enrichment_files(statistics$get_inputs()$get_parameters(), "."))
      if (anyDuplicated(names(files)) || !setequal(names(files), expected)) {
        stop("Final PTM results require every enabled enrichment file.")
      }
      missing <- files[!file.exists(files)]
      if (length(missing)) {
        stop("Missing enrichment file(s): ", paste(missing, collapse = ", "), call. = FALSE)
      }
      private$statistics <- statistics
      private$files <- files[expected]
      private$statistics_sha256 <- statistics_sha256
    },
    #' @description Return the complete statistics component.
    get_statistics = function() private$statistics,
    #' @description Return the completed PTMSEA, KinaseGSEA and MEA of every enabled analysis.
    get_enrichments = function() {
      if (is.null(private$enrichments)) {
        private$enrichments <- .ptm_restore_enrichments(private$statistics, private$files, private$statistics_sha256)
      }
      private$enrichments
    },
    #' @description Return the protsea document of one completed enrichment, as JSON text.
    #' @param method PTMSEA, KinaseGSEA or MEA.
    #' @param analysis DPA, DPU or CF.
    get_enrichment_document = function(method, analysis) {
      key <- .ptm_key(method, toupper(analysis))
      if (!method %in% .ptm_json_enrichment_methods || !key %in% names(private$files)) {
        stop("Enrichment document is not enabled: ", key)
      }
      protsea::read_gsea_json_text(private$files[[key]])
    },
    #' @description Return standard statistics and abundance tables.
    #' @param estimates `"observed"` keeps the DPA, DPU and CF rows whose site
    #'   estimate is observed; `"all"` keeps every row.
    get_tables = function(estimates = c("observed", "all")) private$statistics$get_tables(estimates),
    #' @description Count the DPA, DPU and CF rows of each contrast by site
    #'   estimate type, before the imputed ones are dropped.
    get_estimate_counts = function() private$statistics$get_estimate_counts(),
    #' @description Write the final artifact atomically. The enrichment files
    #'   must lie beside it.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(private$container(path), path)
  )
)

.enabled_ptm_enrichments <- function(parameters) {
  if (!isTRUE(parameters$run_kinase)) {
    return(character())
  }
  .ptm_stage_keys(.ptm_json_enrichment_methods, toupper(names(parameters$analyses)))
}

# Every stage of every analysis, stage by stage.
.ptm_stage_keys <- function(stages, analyses) {
  as.vector(outer(analyses, stages, function(analysis, stage) .ptm_key(stage, analysis)))
}

.load_ptm_results <- function(container, h5mu_path) {
  statistics <- .load_ptm_statistics(container)
  namespace <- container$uns$prophosqua
  relative <- unlist(namespace$enrichment_files)
  files <- stats::setNames(file.path(dirname(h5mu_path), relative), names(relative))
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop(
      "Final PTM results need their enrichment files beside the MuData file; missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  recorded <- unlist(namespace$enrichment_sha256)
  changed <- names(files)[vapply(
    names(files),
    function(key) !identical(.ptm_file_sha256(files[[key]]), as.character(recorded[[key]])),
    logical(1)
  )]
  if (length(changed)) {
    stop(
      "Enrichment files changed since the final MuData was written: ",
      paste(changed, collapse = ", "),
      call. = FALSE
    )
  }
  PTM_results$new(statistics, files, namespace$statistics_sha256)
}
