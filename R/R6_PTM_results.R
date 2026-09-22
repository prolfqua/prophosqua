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
    #' @param enrichments List of complete PTMSEA, KinaseGSEA and MEA objects.
    #' @param enrichment_documents Stored JSON documents when restoring MuData.
    #' @param enrichment_cbor Where the enrichment payload is kept, as
    #'   `list(paths, statistics_sha256)` with paths relative to the MuData file.
    #'   `NULL` keeps the payload inside the MuData file itself.
    initialize = function(statistics, enrichments, enrichment_documents = NULL, enrichment_cbor = NULL) {
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
      if (is.null(enrichment_documents)) {
        enrichment_documents <- lapply(private$enrichments, .ptm_enrichment_document)
      }
      private$enrichment_documents <- .validate_ptm_enrichment_documents(enrichment_documents, expected)
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
      if (length(method) != 1L || !method %in% .ptm_json_enrichment_methods) {
        stop("Unknown enrichment method: ", paste(method, collapse = ", "))
      }
      if (length(analysis) != 1L) {
        stop("Analysis must be one of DPA, DPU or CF.")
      }
      analysis <- toupper(analysis)
      .validate_enrichment_analysis(analysis)
      key <- .ptm_varm_key(method, analysis)
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

.ptm_json_enrichment_methods <- c("PTMSEA", "KinaseGSEA", "MEA")

.enabled_ptm_enrichments <- function(parameters) {
  analyses <- toupper(names(parameters$analyses))
  methods <- if (isTRUE(parameters$run_kinase)) c("PTMSEA", "KinaseGSEA", "MEA") else character()
  as.character(unlist(
    lapply(methods, function(method) {
      vapply(analyses, function(analysis) .ptm_varm_key(method, analysis), character(1))
    }),
    use.names = FALSE
  ))
}

.validate_ptm_cbor_store <- function(store) {
  if (is.null(store)) {
    return(NULL)
  }
  .require_ptm_fields(store, c("paths", "statistics_sha256"), "PTM CBOR store")
  if (!length(store$paths) || is.null(names(store$paths)) || anyDuplicated(names(store$paths))) {
    stop("A PTM CBOR store must name every artifact exactly once.", call. = FALSE)
  }
  list(
    paths = lapply(store$paths, as.character),
    statistics_sha256 = as.character(store$statistics_sha256)
  )
}

# Two storage forms, because a final result either owns CBOR artifacts beside it
# or holds its payload itself. An assembled result names its artifacts and keeps
# nothing: its three MEA payloads alone are half a gigabyte, and every byte of
# them is already in the CBOR files. A result built in memory, as the package
# fixtures are, has nowhere else to put the payload and embeds it.
.final_ptm_container <- function(results) {
  container <- results$get_statistics()$as_container()
  container$uns$prophosqua$stage <- "PTM_results"
  enrichments <- results$get_enrichments()
  documents <- results$get_enrichment_documents()
  container$uns$prophosqua$completed_enrichments <- names(enrichments)
  store <- results$get_enrichment_store()
  if (!is.null(store)) {
    container$uns$prophosqua$enrichment_cbor <- store$paths
    container$uns$prophosqua$statistics_sha256 <- store$statistics_sha256
    return(container)
  }
  for (key in names(enrichments)) {
    branch <- enrichments[[key]]
    modality <- .ptm_analysis_modality(branch$get_analysis())
    namespace <- container$modalities[[modality]]$uns$prophosqua
    completed <- .completed_enrichment_stages(branch)
    for (stage in names(completed)) {
      namespace$completed_stages[[stage]] <- completed[[stage]]
    }
    namespace$enrichment_documents[[key]] <- documents[[key]]
    container$modalities[[modality]]$uns$prophosqua <- namespace
  }
  container
}

.completed_enrichment_stages <- function(stage) {
  source <- stage$get_source()
  preceding <- if (inherits(source, "PTM_statistics")) list() else .completed_enrichment_stages(source)
  if (class(stage)[1L] %in% .ptm_json_enrichment_methods) {
    return(preceding)
  }
  key <- .ptm_varm_key(class(stage)[1L], stage$get_analysis())
  preceding[[key]] <- .pack_ptm_value(stage$get_results())
  preceding
}

# statistics_sha256 binds a stored document to the statistics it was computed
# from, the check the CBOR envelope performs for the preparations. A document
# built in memory has no file to be paired with and carries no hash.
.ptm_enrichment_document <- function(branch, statistics_hash = NULL) {
  method <- class(branch)[1L]
  analysis <- branch$get_analysis()
  result <- branch$get_results()
  builders <- list(
    PTMSEA = function(x) gsea_result_data(x$get_results()$results, category = "PTM-SEA"),
    KinaseGSEA = function(x) gsea_result_data(x$get_results()$gsea_results, category = "KinaseLib")
  )
  object_field <- c(PTMSEA = "results", KinaseGSEA = "gsea_results", MEA = NA_character_)[[method]]
  portable_result <- result
  if (!is.na(object_field)) {
    portable_result[[object_field]] <- NULL
  }
  extension <- list(
    method = method,
    analysis = analysis,
    result_names = names(result),
    result = .pack_ptm_value(portable_result)
  )
  if (!is.null(statistics_hash)) {
    extension$statistics_sha256 <- statistics_hash
  }
  if (identical(method, "MEA")) {
    document <- .append_ptm_json_extension(
      branch$get_source()$get_results()$gsea_json,
      extension
    )
  } else {
    contents <- builders[[method]](branch)
    contents$prophosqua <- extension
    document <- as.character(jsonlite::toJSON(contents, auto_unbox = TRUE, digits = NA, na = "null"))
  }
  wrapper <- list(
    format = "string_gsea",
    version = "1.2.0",
    json = document,
    sha256 = digest::digest(document, algo = "sha256", serialize = FALSE)
  )
  .validate_ptm_enrichment_document(wrapper, .ptm_varm_key(method, analysis))
  wrapper
}

.append_ptm_json_extension <- function(json, extension) {
  source <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(error) {
      stop("Invalid native MEA JSON: ", conditionMessage(error), call. = FALSE)
    }
  )
  .require_ptm_fields(source, c("data", "rank_lists"), "Native MEA JSON")
  if ("prophosqua" %in% names(source)) {
    stop("Native MEA JSON already contains a prophosqua extension.")
  }
  json <- trimws(json)
  if (!startsWith(json, "{") || !endsWith(json, "}")) {
    stop("Native MEA JSON must be one JSON object.")
  }
  extension_json <- as.character(jsonlite::toJSON(
    extension,
    auto_unbox = TRUE,
    digits = NA,
    na = "null"
  ))
  paste0(substr(json, 1L, nchar(json) - 1L), ',"prophosqua":', extension_json, "}")
}

.validate_ptm_enrichment_documents <- function(documents, expected) {
  if (!is.list(documents) || is.null(names(documents))) {
    if (!length(documents) && !length(expected)) {
      return(list())
    }
    stop("Enrichment documents must be a named list.")
  }
  if (anyDuplicated(names(documents)) || !setequal(names(documents), expected)) {
    stop("Final PTM results require every enabled enrichment document.")
  }
  documents <- documents[expected]
  for (key in names(documents)) {
    .validate_ptm_enrichment_document(documents[[key]], key)
  }
  documents
}

.validate_ptm_enrichment_document <- function(document, key) {
  payload <- .ptm_enrichment_payload(document, key)
  method <- .validate_ptm_enrichment_identity(payload, key)
  .validate_ptm_enrichment_contrasts(payload, key, method)
  protsea::decode_gsea_json(document$json)
  invisible(payload)
}

.ptm_enrichment_payload <- function(document, key) {
  fields <- c("format", "version", "json", "sha256")
  .require_ptm_fields(document, fields, paste0("Enrichment document ", key))
  if (length(document) != length(fields) || anyDuplicated(names(document)) || !setequal(names(document), fields)) {
    stop("Enrichment document wrapper has unexpected fields: ", key)
  }
  scalar_text <- vapply(document[fields], function(value) is.character(value) && length(value) == 1L, logical(1))
  if (!all(scalar_text) || !identical(document$format, "string_gsea") || !identical(document$version, "1.2.0")) {
    stop("Unsupported enrichment document: ", key)
  }
  checksum <- digest::digest(document$json, algo = "sha256", serialize = FALSE)
  if (!identical(document$sha256, checksum)) {
    stop("Enrichment document checksum mismatch: ", key)
  }
  payload <- tryCatch(
    jsonlite::fromJSON(document$json, simplifyVector = FALSE),
    error = function(error) stop("Invalid enrichment JSON for ", key, ": ", conditionMessage(error), call. = FALSE)
  )
  .require_ptm_fields(payload, c("data", "rank_lists", "prophosqua"), paste0("Enrichment JSON ", key))
  payload
}

.validate_ptm_enrichment_identity <- function(payload, key) {
  extension <- payload$prophosqua
  .require_ptm_fields(extension, c("method", "analysis", "result_names", "result"), paste0("PTM JSON ", key))
  method <- sub("__.*$", "", key)
  analysis <- utils::URLdecode(sub("^.*?__", "", key))
  if (!identical(extension$method, method) || !identical(extension$analysis, analysis)) {
    stop("Enrichment JSON identity differs from its MuData key: ", key)
  }
  if (.contains_packed_gsea_result(extension$result)) {
    stop("Enrichment JSON contains a serialized gseaResult: ", key)
  }
  method
}

.validate_ptm_enrichment_contrasts <- function(payload, key, method) {
  if (!setequal(names(payload$data), names(payload$rank_lists))) {
    stop("Enrichment JSON contrast names differ between data and rank lists: ", key)
  }
  category <- c(PTMSEA = "PTM-SEA", KinaseGSEA = "KinaseLib", MEA = "MEA")[[method]]
  if (is.null(category)) {
    stop("Unknown enrichment document: ", key)
  }
  for (contrast_name in names(payload$data)) {
    contrast <- payload$data[[contrast_name]]
    .require_ptm_fields(contrast, c("contrast", "gene_pool", "categories"), paste0("Contrast ", contrast_name))
    if (!identical(contrast$contrast, contrast_name) || !identical(names(contrast$categories), category)) {
      stop("Enrichment JSON contrast or category identity is invalid: ", key)
    }
    rank_list <- payload$rank_lists[[contrast_name]]
    .require_ptm_fields(rank_list, c("contrast", "entries"), paste0("Rank list ", contrast_name))
    if (!identical(rank_list$contrast, contrast_name)) {
      stop("Enrichment JSON rank-list identity is invalid: ", key)
    }
    category_data <- contrast$categories[[category]]
    .require_ptm_fields(
      category_data,
      c("category", "contrast", "terms", "gsea_result"),
      paste0("Category ", category)
    )
    if (!identical(category_data$category, category) || !identical(category_data$contrast, contrast_name)) {
      stop("Enrichment JSON category identity is invalid: ", key)
    }
    required_term_fields <- c(
      "term_id",
      "category",
      "description",
      "enrichment_score",
      "direction",
      "fdr",
      "method",
      "genes_mapped",
      "genes_in_set",
      "gene_ids",
      "leading_edge_ids"
    )
    for (term in category_data$terms) {
      .require_ptm_fields(term, required_term_fields, paste0("Enrichment term in ", key))
    }
    .validate_ptm_enrichment_traces(
      category_data,
      length(rank_list$entries),
      key
    )
  }
  invisible(payload)
}

.validate_ptm_enrichment_traces <- function(category_data, rank_count, key) {
  native <- category_data$gsea_result
  .require_ptm_fields(
    native,
    c("result", "gene_sets", "params", "running_scores", "hit_indices"),
    paste0("Native GSEA result in ", key)
  )
  .require_ptm_fields(
    native$result,
    c("columns", "types", "row_names", "row_name_type"),
    paste0("Native GSEA table in ", key)
  )
  .require_ptm_fields(native$result$columns, "ID", paste0("Native GSEA table in ", key))
  term_ids <- vapply(
    category_data$terms,
    function(term) as.character(term$term_id),
    character(1)
  )
  result_ids <- as.character(unlist(native$result$columns$ID, use.names = FALSE))
  if (!identical(result_ids, term_ids)) {
    stop("Native GSEA table and term identifiers differ: ", key)
  }
  .validate_ptm_enrichment_trace_names(native, term_ids, key)
  for (term_id in term_ids) {
    .validate_ptm_enrichment_term_trace(native, term_id, rank_count, key)
  }
  invisible(category_data)
}

.validate_ptm_enrichment_trace_names <- function(native, term_ids, key) {
  for (field in c("running_scores", "hit_indices")) {
    value_names <- names(native[[field]])
    if (is.null(value_names)) {
      value_names <- character()
    }
    if (anyDuplicated(value_names) || !setequal(value_names, term_ids)) {
      stop("Native GSEA ", field, " names differ from term identifiers: ", key)
    }
  }
  invisible(native)
}

.validate_ptm_enrichment_term_trace <- function(native, term_id, rank_count, key) {
  running <- as.numeric(unlist(native$running_scores[[term_id]], use.names = FALSE))
  if (length(running) != rank_count || any(!is.finite(running))) {
    stop("Native GSEA running score is invalid for ", term_id, ": ", key)
  }
  hits <- as.numeric(unlist(native$hit_indices[[term_id]], use.names = FALSE))
  invalid_hits <- any(!is.finite(hits)) ||
    any(hits != floor(hits)) ||
    any(hits < 1 | hits > rank_count) ||
    anyDuplicated(hits) ||
    is.unsorted(hits, strictly = TRUE)
  if (invalid_hits) {
    stop("Native GSEA hit positions are invalid for ", term_id, ": ", key)
  }
  invisible(native)
}

.contains_packed_gsea_result <- function(value) {
  if (!is.list(value)) {
    return(FALSE)
  }
  if (identical(value$type, "gseaResult")) {
    return(TRUE)
  }
  any(vapply(value, .contains_packed_gsea_result, logical(1)))
}

.restore_ptm_enrichment_result <- function(document, key) {
  payload <- .validate_ptm_enrichment_document(document, key)
  result <- .unpack_ptm_value(payload$prophosqua$result)
  method <- payload$prophosqua$method
  object_field <- c(PTMSEA = "results", KinaseGSEA = "gsea_results", MEA = NA_character_)[[method]]
  if (!is.na(object_field)) {
    category <- c(PTMSEA = "PTM-SEA", KinaseGSEA = "KinaseLib")[[method]]
    decoded <- protsea::decode_gsea_json(document$json)
    result[[object_field]] <- lapply(decoded, function(contrast) contrast[[category]])
  }
  result_names <- as.character(unlist(payload$prophosqua$result_names, use.names = FALSE))
  result[result_names]
}

.collect_ptm_enrichment_documents <- function(container, expected) {
  documents <- list()
  locations <- list()
  for (modality in names(container$modalities)) {
    current <- container$modalities[[modality]]$uns$prophosqua$enrichment_documents
    if (!is.null(current)) {
      documents <- c(documents, current)
      locations <- c(locations, stats::setNames(as.list(rep(modality, length(current))), names(current)))
    }
  }
  if (anyDuplicated(names(documents))) {
    stop("MuData contains duplicate enrichment document keys.")
  }
  documents <- .validate_ptm_enrichment_documents(documents, expected)
  for (key in names(documents)) {
    analysis <- utils::URLdecode(sub("^.*?__", "", key))
    if (!identical(locations[[key]], .ptm_analysis_modality(analysis))) {
      stop("Enrichment document is stored in the wrong modality: ", key)
    }
  }
  documents
}

.load_ptm_results <- function(container, h5mu_path = NULL) {
  statistics <- .load_ptm_statistics(container)
  keys <- as.character(container$uns$prophosqua$completed_enrichments)
  expected <- .enabled_ptm_enrichments(statistics$get_inputs()$get_parameters())
  if (anyDuplicated(keys) || !identical(keys, expected)) {
    stop("Final PTM results list does not match enabled enrichments.")
  }
  if (length(container$uns$prophosqua$enrichment_cbor)) {
    return(.load_ptm_results_from_cbor(container, statistics, h5mu_path))
  }
  documents <- .collect_ptm_enrichment_documents(container, expected)
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
  PTM_results$new(statistics, branches, documents)
}

.load_ptm_results_from_cbor <- function(container, statistics, h5mu_path) {
  if (is.null(h5mu_path)) {
    stop("Final PTM results keep their enrichment beside the file; read them with read_ptm_h5mu().", call. = FALSE)
  }
  hash <- as.character(container$uns$prophosqua$statistics_sha256)
  paths <- .ptm_resolve_cbor(container, h5mu_path)
  artifacts <- lapply(paths, .read_ptm_cbor, statistics_hash = hash)
  keys <- vapply(
    artifacts,
    function(artifact) .ptm_varm_key(artifact$stage, artifact$analysis),
    character(1)
  )
  if (!identical(unname(keys), names(paths))) {
    stop("CBOR manifest does not match the artifacts it names.", call. = FALSE)
  }
  PTM_results$new(
    statistics,
    .ptm_branches_from_artifacts(statistics, artifacts),
    .ptm_documents_from_artifacts(artifacts),
    enrichment_cbor = list(
      paths = container$uns$prophosqua$enrichment_cbor,
      statistics_sha256 = hash
    )
  )
}
