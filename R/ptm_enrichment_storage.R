# The enrichment of each analysis is kept in files beside the final MuData,
# never inside it. The completed PTM-SEA, Kinase GSEA and MEA results are
# protsea documents, gzipped JSON; the two kinase-library preparations are
# gzipped CBOR envelopes naming the stage, the analysis and the statistics
# they were computed from.

# The file of each stage in an analysis directory; the pipeline's Snakefile
# names the same files.
.ptm_stage_files <- c(
  PTMSEA = "result_ptm_sea.json.gz",
  KinaseInputs = "intermediate_kinase_inputs.cbor.gz",
  KinaseAssignments = "intermediate_kinase_assignments.cbor.gz",
  KinaseGSEA = "result_kinase_gsea.json.gz",
  MEA = "result_mea.json.gz"
)

# Each completed result: its protsea category, the field holding its gseaResult
# objects, the method protsea records, and the tables derived from them.
.PTM_RESULTS <- list(
  PTMSEA = list(category = "PTM-SEA", objects = "results", method = "fgsea", build = function(x) .ptmsea_result(x)),
  KinaseGSEA = list(
    category = "KinaseLib",
    objects = "gsea_results",
    method = "fgsea",
    build = function(x) .kinasegsea_result(x)
  ),
  MEA = list(category = "MEA", objects = "results", method = "gseapy", build = function(x) .mea_result(x))
)

.ptm_cbor_version <- "1.0.0"

.ptm_file_sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

# Every stage file of every enabled analysis, keyed by stage and analysis.
.ptm_enrichment_files <- function(parameters, root) {
  if (!isTRUE(parameters$run_kinase)) {
    return(character())
  }
  files <- lapply(names(parameters$analyses), function(analysis) {
    directory <- file.path(root, .ptm_parameter(parameters, "analyses", analysis, "subdir"))
    stats::setNames(file.path(directory, .ptm_stage_files), .ptm_key(names(.ptm_stage_files), toupper(analysis)))
  })
  unlist(files)
}

# Written beside the destination and renamed, so a failed write leaves no file.
.ptm_write_atomic <- function(path, write) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".ptm-stage-", tmpdir = dirname(path), fileext = paste0(".", basename(path)))
  on.exit(unlink(temporary), add = TRUE)
  write(temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not write enrichment file: ", path)
  }
  invisible(path)
}

.write_ptm_result <- function(stage, path) {
  spec <- .PTM_RESULTS[[class(stage)[1L]]]
  gsea <- protsea::gsea_result_data(stage$get_results()[[spec$objects]], spec$category, spec$method)
  .ptm_write_atomic(path, function(temporary) protsea::write_gsea_result_json(gsea, temporary))
}

.read_ptm_result <- function(path, method) {
  spec <- .PTM_RESULTS[[method]]
  decoded <- protsea::decode_gsea_json(protsea::read_gsea_json_text(path))
  spec$build(lapply(decoded, function(contrast) contrast[[spec$category]]))
}

.write_ptm_preparation <- function(stage, path, statistics_hash) {
  name <- class(stage)[1L]
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
  .ptm_write_atomic(path, function(temporary) {
    connection <- gzfile(temporary, "wb")
    on.exit(close(connection))
    writeBin(secretbase::cborenc(artifact), connection)
  })
}

.read_ptm_preparation <- function(path, stage, analysis, statistics_hash) {
  if (is.null(path)) {
    stop("Missing ", stage, " input for ", analysis)
  }
  compressed <- readBin(path, what = "raw", n = file.info(path)$size)
  artifact <- secretbase::cbordec(memDecompress(compressed, type = "gzip"))
  .require_ptm_fields(artifact, c("format", "version", "stage", "analysis", "statistics_sha256", "result"), path)
  if (!identical(artifact$format, "prophosqua_stage") || !identical(artifact$version, .ptm_cbor_version)) {
    stop("Unsupported PTM CBOR artifact: ", path)
  }
  if (!identical(artifact$stage, stage) || !identical(artifact$analysis, analysis)) {
    stop("Wrong PTM CBOR stage or analysis: ", path)
  }
  if (!identical(artifact$statistics_sha256, statistics_hash)) {
    stop("PTM CBOR artifact was produced from different statistics: ", path)
  }
  .unpack_ptm_value(artifact$result)
}

# A stage in the file the pipeline writes it to; the examples and tests write
# every stage, the pipeline leaves the kinase-library stages to that tool.
.write_ptm_stage_file <- function(stage, path, statistics_hash) {
  if (class(stage)[1L] %in% names(.PTM_RESULTS)) {
    return(.write_ptm_result(stage, path))
  }
  .write_ptm_preparation(stage, path, statistics_hash)
}

# Every enrichment of every analysis, restored from its files on the one
# statistics object it was computed from.
.ptm_restore_enrichments <- function(statistics, files, statistics_hash) {
  if (!length(files)) {
    return(list())
  }
  analyses <- toupper(names(statistics$get_inputs()$get_parameters()$analyses))
  enrichments <- lapply(analyses, function(analysis) {
    file <- function(stage) files[[.ptm_key(stage, analysis)]]
    preparation <- function(stage) .read_ptm_preparation(file(stage), stage, analysis, statistics_hash)
    assignments <- KinaseAssignments$new(
      KinaseInputs$new(statistics, analysis, preparation("KinaseInputs")),
      analysis,
      preparation("KinaseAssignments")
    )
    list(
      PTMSEA$new(statistics, analysis, .read_ptm_result(file("PTMSEA"), "PTMSEA")),
      KinaseGSEA$new(assignments, analysis, .read_ptm_result(file("KinaseGSEA"), "KinaseGSEA")),
      MEA$new(assignments, analysis, .read_ptm_result(file("MEA"), "MEA"))
    )
  })
  enrichments <- Reduce(c, enrichments, list())
  keys <- vapply(enrichments, function(x) .ptm_key(class(x)[1L], x$get_analysis()), character(1))
  stats::setNames(enrichments, keys)[.enabled_ptm_enrichments(statistics$get_inputs()$get_parameters())]
}

#' Compute one enrichment result, or the kinase-library inputs, of one analysis
#'
#' PTM-SEA and Kinase GSEA are written as protsea documents, gzipped JSON; the
#' kinase-library inputs as a gzipped CBOR envelope. The MEA is computed by the
#' kinase-library tool, which writes its own protsea document.
#' @param statistics_h5mu The shared, read-only statistics MuData file.
#' @param output Output file, as the pipeline lays out an analysis directory.
#' @param stage One of PTMSEA, KinaseInputs or KinaseGSEA.
#' @param analysis DPA, DPU or CF.
#' @param preparation Named paths to the KinaseInputs and KinaseAssignments
#'   files, for KinaseGSEA.
#' @return Output path, invisibly.
#' @export
compute_ptm_enrichment <- function(statistics_h5mu, output, stage, analysis, preparation = list()) {
  .validate_enrichment_analysis(analysis)
  statistics <- read_ptm_h5mu(statistics_h5mu, PTM_statistics)
  statistics_hash <- .ptm_file_sha256(statistics_h5mu)
  if (identical(stage, "KinaseInputs")) {
    return(.write_ptm_preparation(KinaseInputs$new(statistics, analysis), output, statistics_hash))
  }
  if (identical(stage, "PTMSEA")) {
    return(.write_ptm_result(PTMSEA$new(statistics, analysis), output))
  }
  if (identical(stage, "KinaseGSEA")) {
    prepared <- function(stage) .read_ptm_preparation(preparation[[stage]], stage, analysis, statistics_hash)
    assignments <- KinaseAssignments$new(
      KinaseInputs$new(statistics, analysis, prepared("KinaseInputs")),
      analysis,
      prepared("KinaseAssignments")
    )
    return(.write_ptm_result(KinaseGSEA$new(assignments, analysis), output))
  }
  stop("Unsupported computed enrichment stage: ", stage)
}

#' Assemble the final MuData from the statistics and the enrichment files
#'
#' The final MuData holds the statistics and the names and checksums of the
#' enrichment files, which must lie beside it as the pipeline lays them out.
#' @param statistics_h5mu Complete statistics artifact.
#' @param output_h5mu Final artifact.
#' @return The complete final object, invisibly.
#' @export
assemble_ptm_results <- function(statistics_h5mu, output_h5mu) {
  statistics <- read_ptm_h5mu(statistics_h5mu, PTM_statistics)
  files <- .ptm_enrichment_files(statistics$get_inputs()$get_parameters(), dirname(output_h5mu))
  results <- PTM_results$new(statistics, files, .ptm_file_sha256(statistics_h5mu))
  results$write_h5mu(output_h5mu)
  invisible(results)
}

.ptm_relative_paths <- function(paths, output_h5mu) {
  prefix <- paste0(normalizePath(dirname(output_h5mu), mustWork = TRUE), .Platform$file.sep)
  vapply(
    paths,
    function(path) {
      full <- normalizePath(path, mustWork = TRUE)
      if (!startsWith(full, prefix)) {
        stop("Enrichment file lies outside the final MuData directory: ", path, call. = FALSE)
      }
      substring(full, nchar(prefix) + 1L)
    },
    character(1)
  )
}
