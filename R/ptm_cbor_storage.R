# Compact Snakemake handoffs, and the only place the enrichment payload is kept.
# The final MuData records where these artifacts are and reads them back; it does
# not copy them, because one MEA payload alone is a quarter of a gigabyte.
#
# Two artifact kinds, because two kinds of thing are being stored. A completed
# PTM-SEA, Kinase GSEA or MEA stage is a string_gsea document, whose format,
# writer and reader belong to protsea; it is written as gzipped JSON through
# protsea and stays readable by anything that knows the format. The kinase
# preparations are packed R structures with no such format, and stay gzipped
# CBOR. Both are gzipped: the payload is text and compresses about 2.5x.
.ptm_cbor_version <- "1.0.0"
.ptm_cbor_stages <- c(
  "PTMSEA",
  "KinaseInputs",
  "KinaseAssignments",
  "MotifEnrichment",
  "KinaseGSEA",
  "MEA"
)

.ptm_statistics_hash <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

.is_ptm_json_artifact <- function(path) grepl("\\.json(\\.gz)?$", path)

.write_ptm_cbor <- function(stage, path, statistics_hash) {
  name <- class(stage)[1L]
  if (!name %in% .ptm_cbor_stages) {
    stop("Unsupported stage artifact: ", name)
  }
  if (.is_ptm_json_artifact(path) != (name %in% .ptm_json_enrichment_methods)) {
    stop("Stage ", name, " does not belong in ", basename(path))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = ".ptm-stage-",
    tmpdir = dirname(path),
    fileext = if (grepl("\\.gz$", path)) ".gz" else ""
  )
  on.exit(unlink(temporary), add = TRUE)
  if (name %in% .ptm_json_enrichment_methods) {
    protsea::write_gsea_json_text(.ptm_enrichment_document(stage, statistics_hash)$json, temporary)
  } else {
    artifact <- list(
      format = "prophosqua_stage",
      version = .ptm_cbor_version,
      stage = name,
      analysis = stage$get_analysis(),
      statistics_sha256 = statistics_hash,
      result = .pack_ptm_value(stage$get_results())
    )
    if (identical(name, "KinaseInputs")) {
      artifact$settings <- stage$get_statistics()$get_inputs()$get_parameters()$kinaselib
    }
    # gzfile writes a real gzip container; memCompress("gzip") would emit a bare
    # zlib stream that no gzip reader accepts, and the .gz name would be a lie.
    connection <- gzfile(temporary, "wb")
    writeBin(secretbase::cborenc(artifact), connection)
    close(connection)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not write stage artifact: ", path)
  }
  invisible(path)
}

# A stored document describes itself: prophosqua's extension names the method,
# the analysis and the statistics it was computed from, so the envelope the
# CBOR artifacts carry is not needed beside it.
.read_ptm_json_artifact <- function(path, statistics_hash) {
  json <- protsea::read_gsea_json_text(path)
  document <- list(
    format = "string_gsea",
    version = "1.2.0",
    json = json,
    sha256 = digest::digest(json, algo = "sha256", serialize = FALSE)
  )
  payload <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(error) stop("Invalid enrichment JSON: ", path, call. = FALSE)
  )
  extension <- payload$prophosqua
  .require_ptm_fields(extension, c("method", "analysis", "statistics_sha256"), paste0("PTM JSON ", path))
  if (!identical(as.character(extension$statistics_sha256), statistics_hash)) {
    stop("PTM stage artifact was produced from different statistics: ", path)
  }
  .validate_enrichment_analysis(extension$analysis)
  .validate_ptm_enrichment_document(document, .ptm_varm_key(extension$method, extension$analysis))
  list(
    format = "prophosqua_stage",
    version = .ptm_cbor_version,
    stage = as.character(extension$method),
    analysis = as.character(extension$analysis),
    statistics_sha256 = statistics_hash,
    document = document
  )
}

.read_ptm_cbor <- function(path, statistics_hash) {
  if (.is_ptm_json_artifact(path)) {
    return(.read_ptm_json_artifact(path, statistics_hash))
  }
  compressed <- readBin(path, what = "raw", n = file.info(path)$size)
  artifact <- secretbase::cbordec(memDecompress(compressed, type = "gzip"))
  .require_ptm_fields(
    artifact,
    c("format", "version", "stage", "analysis", "statistics_sha256"),
    "PTM CBOR artifact"
  )
  if (!identical(artifact$format, "prophosqua_stage") || !identical(artifact$version, .ptm_cbor_version)) {
    stop("Unsupported PTM CBOR artifact: ", path)
  }
  if (!artifact$stage %in% .ptm_cbor_stages) {
    stop("Unknown PTM CBOR stage: ", artifact$stage)
  }
  .validate_enrichment_analysis(artifact$analysis)
  if (!identical(artifact$statistics_sha256, statistics_hash)) {
    stop("PTM CBOR artifact was produced from different statistics: ", path)
  }
  key <- .ptm_varm_key(artifact$stage, artifact$analysis)
  field <- if (artifact$stage %in% .ptm_json_enrichment_methods) "document" else "result"
  .require_ptm_fields(artifact, field, key)
  allowed <- c("format", "version", "stage", "analysis", "statistics_sha256", field)
  if (identical(artifact$stage, "KinaseInputs")) {
    .require_ptm_fields(artifact, "settings", key)
    allowed <- c(allowed, "settings")
  }
  if (length(artifact) != length(allowed) || !setequal(names(artifact), allowed)) {
    stop("Unexpected PTM CBOR fields: ", key)
  }
  if (identical(field, "document")) {
    .validate_ptm_enrichment_document(artifact$document, key)
    artifact$document <- artifact$document[c("format", "version", "json", "sha256")]
  }
  artifact
}

.read_ptm_preparation <- function(paths, stage, analysis, statistics_hash) {
  path <- paths[[stage]]
  if (is.null(path)) {
    stop("Missing ", stage, " CBOR input for ", analysis)
  }
  artifact <- .read_ptm_cbor(path, statistics_hash)
  if (!identical(artifact$stage, stage) || !identical(artifact$analysis, analysis)) {
    stop("Wrong PTM CBOR stage or analysis: ", path)
  }
  .unpack_ptm_value(artifact$result)
}

#' Compute one complete enrichment stage into a compact CBOR artifact
#' @param statistics_h5mu The shared, read-only statistics MuData file.
#' @param output_cbor Output CBOR artifact.
#' @param stage One of PTMSEA, KinaseInputs, KinaseGSEA or MEA.
#' @param analysis DPA, DPU or CF.
#' @param preparation_cbor Named paths to required preparation CBOR artifacts.
#' @return Output path, invisibly.
#' @export
compute_ptm_enrichment_cbor <- function(
  statistics_h5mu,
  output_cbor,
  stage,
  analysis,
  preparation_cbor = list()
) {
  statistics <- read_ptm_h5mu(statistics_h5mu, PTM_statistics)
  .validate_enrichment_analysis(analysis)
  statistics_hash <- .ptm_statistics_hash(statistics_h5mu)
  Type <- list(PTMSEA = PTMSEA, KinaseInputs = KinaseInputs, KinaseGSEA = KinaseGSEA, MEA = MEA)[[stage]]
  if (is.null(Type)) {
    stop("Unsupported computed CBOR stage: ", stage)
  }
  source <- statistics
  if (stage %in% c("KinaseGSEA", "MEA")) {
    kinase_inputs <- KinaseInputs$new(
      statistics,
      analysis,
      .read_ptm_preparation(preparation_cbor, "KinaseInputs", analysis, statistics_hash)
    )
    source <- KinaseAssignments$new(
      kinase_inputs,
      analysis,
      .read_ptm_preparation(preparation_cbor, "KinaseAssignments", analysis, statistics_hash)
    )
  }
  if (identical(stage, "MEA")) {
    source <- MotifEnrichment$new(
      source,
      analysis,
      .read_ptm_preparation(preparation_cbor, "MotifEnrichment", analysis, statistics_hash)
    )
  }
  .write_ptm_cbor(source$build(Type, analysis = analysis), output_cbor, statistics_hash)
}

#' Assemble a final MuData file from compact CBOR stage artifacts
#' @param statistics_h5mu Complete statistics artifact.
#' @param enrichment_cbor Every enabled CBOR artifact, including kinase preparations.
#' @param output_h5mu Final artifact.
#' @return The complete final object, invisibly.
#' @export
assemble_ptm_cbor <- function(statistics_h5mu, enrichment_cbor, output_h5mu) {
  statistics <- read_ptm_h5mu(statistics_h5mu, PTM_statistics)
  statistics_hash <- .ptm_statistics_hash(statistics_h5mu)
  artifacts <- lapply(enrichment_cbor, .read_ptm_cbor, statistics_hash = statistics_hash)
  names(artifacts) <- vapply(
    artifacts,
    function(x) .ptm_varm_key(x$stage, x$analysis),
    character(1)
  )
  enabled <- .enabled_ptm_enrichments(statistics$get_inputs()$get_parameters())
  analyses <- toupper(names(statistics$get_inputs()$get_parameters()$analyses))
  expected <- if (length(enabled)) {
    as.character(unlist(
      lapply(
        .ptm_cbor_stages,
        function(stage) vapply(analyses, function(analysis) .ptm_varm_key(stage, analysis), character(1))
      ),
      use.names = FALSE
    ))
  } else {
    character()
  }
  if (anyDuplicated(names(artifacts)) || !setequal(names(artifacts), expected)) {
    stop("Final PTM results require every enabled CBOR stage exactly once.")
  }
  results <- PTM_results$new(
    statistics,
    .ptm_branches_from_artifacts(statistics, artifacts),
    .ptm_documents_from_artifacts(artifacts),
    enrichment_cbor = list(
      paths = stats::setNames(
        as.list(.ptm_relative_cbor(enrichment_cbor, output_h5mu)),
        names(artifacts)
      ),
      statistics_sha256 = statistics_hash
    )
  )
  results$write_h5mu(output_h5mu)
  invisible(results)
}

.ptm_relative_cbor <- function(paths, output_h5mu) {
  dir.create(dirname(output_h5mu), recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(dirname(output_h5mu), mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  vapply(
    paths,
    function(path) {
      full <- normalizePath(path, mustWork = TRUE)
      if (!startsWith(full, prefix)) {
        stop("CBOR artifact lies outside the final MuData directory: ", path, call. = FALSE)
      }
      substring(full, nchar(prefix) + 1L)
    },
    character(1),
    USE.NAMES = FALSE
  )
}

.ptm_resolve_cbor <- function(container, h5mu_path) {
  manifest <- container$uns$prophosqua$enrichment_cbor
  if (!length(manifest)) {
    return(character())
  }
  resolved <- stats::setNames(
    file.path(dirname(h5mu_path), as.character(unlist(manifest, use.names = FALSE))),
    names(manifest)
  )
  missing <- resolved[!file.exists(resolved)]
  if (length(missing)) {
    stop(
      "Final PTM results need their CBOR artifacts beside the MuData file; missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  resolved
}

# The stored document is the payload, not a rendering of it: restoring a result
# and writing it again does not reproduce the same JSON byte for byte.
.ptm_documents_from_artifacts <- function(artifacts) {
  json <- Filter(function(artifact) artifact$stage %in% .ptm_json_enrichment_methods, artifacts)
  lapply(json, function(artifact) artifact$document)
}

.ptm_branches_from_artifacts <- function(statistics, artifacts) {
  analyses <- unique(vapply(artifacts, function(artifact) artifact$analysis, character(1)))
  unlist(
    lapply(analyses, function(analysis) {
      artifact <- function(stage) {
        key <- .ptm_varm_key(stage, analysis)
        if (is.null(artifacts[[key]])) {
          stop("Missing CBOR stage: ", key, call. = FALSE)
        }
        artifacts[[key]]
      }
      restore <- function(stage) {
        .restore_ptm_enrichment_result(artifact(stage)$document, .ptm_varm_key(stage, analysis))
      }
      inputs <- KinaseInputs$new(statistics, analysis, .unpack_ptm_value(artifact("KinaseInputs")$result))
      assignments <- KinaseAssignments$new(
        inputs,
        analysis,
        .unpack_ptm_value(artifact("KinaseAssignments")$result)
      )
      motif <- MotifEnrichment$new(
        assignments,
        analysis,
        .unpack_ptm_value(artifact("MotifEnrichment")$result)
      )
      list(
        PTMSEA$new(statistics, analysis, restore("PTMSEA")),
        KinaseGSEA$new(assignments, analysis, restore("KinaseGSEA")),
        MEA$new(motif, analysis, restore("MEA"))
      )
    }),
    recursive = FALSE
  )
}
