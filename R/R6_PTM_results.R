#' Final PTM results with every enabled enrichment completed
#' @export
PTM_results <- R6::R6Class(
  "PTM_results",
  private = list(
    statistics = NULL,
    enrichments = NULL,
    enrichment_documents = NULL,
    enrichment_store = NULL
  ),
  public = list(
    #' @description Assemble all enabled analyses from their completed stages.
    #' @param statistics Complete PTM_statistics stage.
    #' @param enrichments List of complete PTMSEA, KinaseGSEA and MEA objects,
    #'   computed from `statistics`.
    #' @param enrichment_documents Stored JSON documents when restoring MuData.
    #' @param enrichment_cbor Where the enrichment payload is kept, as
    #'   `list(paths, statistics_sha256)` with paths relative to the MuData file.
    #'   `NULL` keeps the payload inside the MuData file itself.
    initialize = function(statistics, enrichments, enrichment_documents = NULL, enrichment_cbor = NULL) {
      stopifnot(inherits(statistics, "PTM_statistics"))
      keys <- vapply(enrichments, function(value) .ptm_key(class(value)[1L], value$get_analysis()), character(1))
      expected <- .enabled_ptm_enrichments(statistics$get_inputs()$get_parameters())
      if (anyDuplicated(keys) || !setequal(keys, expected)) {
        stop("Final PTM results require every enabled enrichment.")
      }
      if (!all(vapply(enrichments, function(value) identical(value$get_statistics(), statistics), logical(1)))) {
        stop("PTM enrichments were computed from different statistics.")
      }
      private$statistics <- statistics
      private$enrichments <- stats::setNames(enrichments, keys)[expected]
      if (is.null(enrichment_documents)) {
        enrichment_documents <- lapply(private$enrichments, .ptm_enrichment_document)
      }
      if (anyDuplicated(names(enrichment_documents)) || !setequal(names(enrichment_documents), expected)) {
        stop("Final PTM results require every enabled enrichment document.")
      }
      for (key in expected) {
        .validate_ptm_enrichment_document(enrichment_documents[[key]], key)
      }
      private$enrichment_documents <- enrichment_documents[expected]
      private$enrichment_store <- .validate_ptm_cbor_store(enrichment_cbor)
    },
    #' @description Return where the enrichment payload is kept, or NULL.
    get_enrichment_store = function() private$enrichment_store,
    #' @description Return the complete statistics component.
    get_statistics = function() private$statistics,
    #' @description Return the collection of complete enabled enrichment objects.
    get_enrichments = function() private$enrichments,
    #' @description Return every complete JSON enrichment document.
    get_enrichment_documents = function() private$enrichment_documents,
    #' @description Return one complete JSON enrichment document.
    #' @param method PTMSEA, KinaseGSEA or MEA.
    #' @param analysis DPA, DPU or CF.
    get_enrichment_document = function(method, analysis) {
      key <- .ptm_key(method, toupper(analysis))
      if (!key %in% names(private$enrichment_documents)) {
        stop("Enrichment document is not enabled: ", key)
      }
      private$enrichment_documents[[key]]
    },
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
  if (!isTRUE(parameters$run_kinase)) {
    return(character())
  }
  .ptm_stage_keys(.ptm_json_enrichment_methods, toupper(names(parameters$analyses)))
}

# Every stage of every analysis, stage by stage.
.ptm_stage_keys <- function(stages, analyses) {
  as.vector(outer(analyses, stages, function(analysis, stage) .ptm_key(stage, analysis)))
}

.validate_ptm_cbor_store <- function(store) {
  if (is.null(store)) {
    return(NULL)
  }
  .require_ptm_fields(store, c("paths", "statistics_sha256"), "PTM CBOR store")
  if (!length(store$paths) || is.null(names(store$paths)) || anyDuplicated(names(store$paths))) {
    stop("A PTM CBOR store must name every artifact exactly once.", call. = FALSE)
  }
  list(paths = lapply(store$paths, as.character), statistics_sha256 = as.character(store$statistics_sha256))
}

# Two storage forms, because a final result either owns CBOR artifacts beside it
# or holds its payload itself. An assembled result names its artifacts and keeps
# nothing: its three MEA payloads alone are half a gigabyte, and every byte of
# them is already in the CBOR files. A result built in memory, as the package
# example is, has nowhere else to put the payload and embeds it in the uns of
# each analysis' modality.
.final_ptm_container <- function(results) {
  container <- results$get_statistics()$as_container()
  container$uns$prophosqua$stage <- "PTM_results"
  enrichments <- results$get_enrichments()
  container$uns$prophosqua$completed_enrichments <- names(enrichments)
  store <- results$get_enrichment_store()
  if (!is.null(store)) {
    container$uns$prophosqua$enrichment_cbor <- store$paths
    container$uns$prophosqua$statistics_sha256 <- store$statistics_sha256
    return(container)
  }
  documents <- results$get_enrichment_documents()
  for (key in names(enrichments)) {
    branch <- enrichments[[key]]
    modality <- .ptm_analysis_modality(branch$get_analysis())
    namespace <- container$modalities[[modality]]$uns$prophosqua
    completed <- .completed_enrichment_stages(branch$get_source())
    namespace$completed_stages[names(completed)] <- completed
    namespace$enrichment_documents[[key]] <- documents[[key]]
    container$modalities[[modality]]$uns$prophosqua <- namespace
  }
  container
}

# The kinase preparations a completed enrichment was computed from, packed and
# keyed by stage and analysis.
.completed_enrichment_stages <- function(stage) {
  if (inherits(stage, "PTM_statistics")) {
    return(list())
  }
  preceding <- .completed_enrichment_stages(stage$get_source())
  preceding[[.ptm_key(class(stage)[1L], stage$get_analysis())]] <- .pack_ptm_value(stage$get_results())
  preceding
}

# A completed enrichment as a string_gsea document, the prophosqua extension
# carrying the rest of its result. statistics_sha256 binds a stored document to
# the statistics it was computed from; a document built in memory has no file
# to be paired with and carries no hash.
.ptm_enrichment_document <- function(branch, statistics_hash = NULL) {
  method <- class(branch)[1L]
  result <- branch$get_results()
  object_field <- c(PTMSEA = "results", KinaseGSEA = "gsea_results", MEA = NA_character_)[[method]]
  portable_result <- result
  if (!is.na(object_field)) {
    portable_result[[object_field]] <- NULL
  }
  extension <- list(
    method = method,
    analysis = branch$get_analysis(),
    result_names = names(result),
    result = .pack_ptm_value(portable_result)
  )
  if (!is.null(statistics_hash)) {
    extension$statistics_sha256 <- statistics_hash
  }
  if (identical(method, "MEA")) {
    document <- .append_ptm_json_extension(branch$get_source()$get_results()$gsea_json, extension)
  } else {
    category <- c(PTMSEA = "PTM-SEA", KinaseGSEA = "KinaseLib")[[method]]
    contents <- protsea::gsea_result_data(result[[object_field]], category = category)
    contents$prophosqua <- extension
    document <- protsea::gsea_result_json_text(contents)
  }
  list(
    format = "string_gsea",
    version = "1.2.0",
    json = document,
    sha256 = digest::digest(document, algo = "sha256", serialize = FALSE)
  )
}

# The MEA document is the kinase-library tool's own string_gsea JSON, kept byte
# for byte with the extension spliced in as its last member.
.append_ptm_json_extension <- function(json, extension) {
  json <- trimws(json)
  if (!startsWith(json, "{") || !endsWith(json, "}")) {
    stop("Native MEA JSON must be one JSON object.")
  }
  paste0(substr(json, 1L, nchar(json) - 1L), ',"prophosqua":', protsea::gsea_result_json_text(extension), "}")
}

.validate_ptm_enrichment_document <- function(document, key) {
  fields <- c("format", "version", "json", "sha256")
  .require_ptm_fields(document, fields, paste0("Enrichment document ", key))
  if (length(document) != length(fields) || !setequal(names(document), fields)) {
    stop("Enrichment document wrapper has unexpected fields: ", key)
  }
  if (!identical(document$format, "string_gsea") || !identical(document$version, "1.2.0")) {
    stop("Unsupported enrichment document: ", key)
  }
  if (!identical(document$sha256, digest::digest(document$json, algo = "sha256", serialize = FALSE))) {
    stop("Enrichment document checksum mismatch: ", key)
  }
  invisible(document)
}

.restore_ptm_enrichment_result <- function(document, key) {
  payload <- tryCatch(
    jsonlite::fromJSON(document$json, simplifyVector = FALSE),
    error = function(error) stop("Invalid enrichment JSON for ", key, ": ", conditionMessage(error), call. = FALSE)
  )
  extension <- payload$prophosqua
  .require_ptm_fields(extension, c("method", "analysis", "result_names", "result"), paste0("PTM JSON ", key))
  if (!identical(.ptm_key(extension$method, extension$analysis), key)) {
    stop("Enrichment JSON identity differs from its MuData key: ", key)
  }
  result <- .unpack_ptm_value(extension$result)
  object_field <- c(PTMSEA = "results", KinaseGSEA = "gsea_results", MEA = NA_character_)[[extension$method]]
  if (!is.na(object_field)) {
    category <- c(PTMSEA = "PTM-SEA", KinaseGSEA = "KinaseLib")[[extension$method]]
    result[[object_field]] <- lapply(protsea::decode_gsea_json(document$json), function(contrast) contrast[[category]])
  }
  result[as.character(unlist(extension$result_names, use.names = FALSE))]
}

.load_ptm_results <- function(container, h5mu_path = NULL) {
  statistics <- .load_ptm_statistics(container)
  keys <- as.character(container$uns$prophosqua$completed_enrichments)
  if (!identical(keys, .enabled_ptm_enrichments(statistics$get_inputs()$get_parameters()))) {
    stop("Final PTM results list does not match enabled enrichments.")
  }
  if (!length(container$uns$prophosqua$enrichment_cbor)) {
    artifacts <- .embedded_ptm_artifacts(container)
    return(PTM_results$new(
      statistics,
      .ptm_branches_from_artifacts(statistics, artifacts),
      .ptm_documents_from_artifacts(artifacts)
    ))
  }
  if (is.null(h5mu_path)) {
    stop("Final PTM results keep their enrichment beside the file; read them with read_ptm_h5mu().", call. = FALSE)
  }
  hash <- as.character(container$uns$prophosqua$statistics_sha256)
  paths <- .ptm_resolve_cbor(container, h5mu_path)
  artifacts <- lapply(paths, .read_ptm_cbor, statistics_hash = hash)
  stored <- vapply(artifacts, function(x) .ptm_key(x$stage, x$analysis), character(1), USE.NAMES = FALSE)
  if (!identical(stored, names(paths))) {
    stop("CBOR manifest does not match the artifacts it names.", call. = FALSE)
  }
  PTM_results$new(
    statistics,
    .ptm_branches_from_artifacts(statistics, artifacts),
    .ptm_documents_from_artifacts(artifacts),
    enrichment_cbor = list(paths = container$uns$prophosqua$enrichment_cbor, statistics_sha256 = hash)
  )
}

# The embedded payload in the shape of the CBOR artifacts: a packed result for
# each kinase preparation, a document for each completed enrichment.
.embedded_ptm_artifacts <- function(container) {
  artifacts <- list()
  for (modality in container$modalities) {
    namespace <- modality$uns$prophosqua
    for (key in names(namespace$completed_stages)) {
      artifacts[[key]] <- list(result = namespace$completed_stages[[key]])
    }
    for (key in names(namespace$enrichment_documents)) {
      artifacts[[key]] <- list(document = namespace$enrichment_documents[[key]])
    }
  }
  artifacts
}
