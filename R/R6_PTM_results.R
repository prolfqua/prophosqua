#' Final PTM results with every enabled enrichment completed
#' @export
PTM_results <- R6::R6Class(
  "PTM_results",
  private = list(statistics = NULL, enrichments = NULL),
  public = list(
    #' @description Assemble all enabled analyses from their completed stages.
    #' @param statistics Complete PTM_statistics stage.
    #' @param enrichments List of complete PTMSEA, KinaseGSEA and MEA objects.
    initialize = function(statistics, enrichments) {
      stopifnot(inherits(statistics, "PTM_statistics"))
      keys <- vapply(enrichments, function(value) .ptm_varm_key(class(value)[1L], value$get_analysis()), character(1))
      expected <- .enabled_ptm_enrichments(statistics$get_inputs()$get_parameters())
      if (anyDuplicated(keys) || !setequal(keys, expected)) {
        stop("Final PTM results require every enabled enrichment.")
      }
      for (value in enrichments) {
        .require_same_ptm_inputs(statistics$get_inputs(), value$get_statistics()$get_inputs())
      }
      private$statistics <- statistics
      private$enrichments <- stats::setNames(enrichments, keys)[expected]
    },
    #' @description Return the complete statistics component.
    get_statistics = function() private$statistics,
    #' @description Return the collection of complete enabled enrichment objects.
    get_enrichments = function() private$enrichments,
    #' @description Return standard statistics and abundance tables.
    get_tables = function() private$statistics$get_tables(),
    #' @description Return a detached storage representation.
    as_container = function() .final_ptm_container(self),
    #' @description Write the final artifact atomically.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(self$as_container(), path)
  )
)

.enabled_ptm_enrichments <- function(parameters) {
  analyses <- toupper(names(parameters$analyses))
  methods <- if (isTRUE(parameters$run_kinase)) c("PTMSEA", "KinaseGSEA", "MEA") else character()
  unlist(
    lapply(methods, function(method) {
      vapply(analyses, function(analysis) .ptm_varm_key(method, analysis), character(1))
    }),
    use.names = FALSE
  )
}

.final_ptm_container <- function(results) {
  container <- results$get_statistics()$as_container()
  container$uns$prophosqua$stage <- "PTM_results"
  enrichments <- results$get_enrichments()
  container$uns$prophosqua$completed_enrichments <- names(enrichments)
  for (key in names(enrichments)) {
    branch <- enrichments[[key]]
    modality <- .ptm_analysis_modality(branch$get_analysis())
    namespace <- container$modalities[[modality]]$uns$prophosqua
    completed <- .completed_enrichment_stages(branch)
    for (stage in names(completed)) {
      namespace$completed_stages[[stage]] <- completed[[stage]]
    }
    namespace$enrichment_documents[[key]] <- .ptm_enrichment_document(branch)
    container$modalities[[modality]]$uns$prophosqua <- namespace
  }
  container
}

.completed_enrichment_stages <- function(stage) {
  source <- stage$get_source()
  preceding <- if (inherits(source, "PTM_statistics")) list() else .completed_enrichment_stages(source)
  key <- .ptm_varm_key(class(stage)[1L], stage$get_analysis())
  preceding[[key]] <- .pack_ptm_value(stage$get_results())
  preceding
}

.ptm_enrichment_document <- function(branch) {
  builders <- list(
    PTMSEA = function(x) gsea_result_data(x$get_results()$results, category = "PTM-SEA"),
    KinaseGSEA = function(x) gsea_result_data(x$get_results()$gsea_results, category = "KinaseLib"),
    MEA = function(x) {
      assignments <- x$get_source()$get_source()
      ranks <- lapply(assignments$get_source()$get_results()$ranks, function(table) {
        stats::setNames(table$statistic.site, table$SequenceWindow)
      })
      mea_gsea_result_data(x$get_results()$mea_clean, ranks, assignments$get_results()$term2gene)
    }
  )
  document <- as.character(jsonlite::toJSON(builders[[class(branch)[1L]]](branch), auto_unbox = TRUE, digits = NA))
  list(
    format = "string_gsea",
    version = "1.1.0",
    json = document,
    sha256 = digest::digest(document, algo = "sha256", serialize = FALSE)
  )
}

.load_ptm_results <- function(container) {
  statistics <- .load_ptm_statistics(container)
  keys <- as.character(container$uns$prophosqua$completed_enrichments)
  readers <- list(PTMSEA = .load_ptmsea, KinaseGSEA = .load_kinase_gsea, MEA = .load_mea)
  branches <- lapply(keys, function(key) {
    method <- sub("__.*$", "", key)
    container$uns$prophosqua$analysis <- utils::URLdecode(sub("^.*?__", "", key))
    reader <- readers[[method]]
    if (is.null(reader)) {
      stop("Unknown completed enrichment: ", key)
    }
    reader(container)
  })
  PTM_results$new(statistics, branches)
}

#' Assemble completed statistics and enrichment MuData artifacts
#' @param statistics_h5mu Complete statistics artifact.
#' @param enrichment_h5mu Paths to every enabled completed enrichment artifact.
#' @param output_h5mu Final artifact.
#' @return The complete final object, invisibly.
#' @export
assemble_ptm_h5mu <- function(statistics_h5mu, enrichment_h5mu, output_h5mu) {
  statistics <- read_ptm_h5mu(statistics_h5mu, PTM_statistics)
  branches <- lapply(enrichment_h5mu, read_ptm_h5mu)
  result <- PTM_results$new(statistics, branches)
  result$write_h5mu(output_h5mu)
  invisible(result)
}
