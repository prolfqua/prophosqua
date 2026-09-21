# Compact Snakemake handoffs. The final enrichment payload is still the
# validated string_gsea JSON document; CBOR is only its on-disk transport.
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

.write_ptm_cbor <- function(stage, path, statistics_hash) {
  name <- class(stage)[1L]
  if (!name %in% .ptm_cbor_stages) {
    stop("Unsupported CBOR stage: ", name)
  }
  artifact <- list(
    format = "prophosqua_stage",
    version = .ptm_cbor_version,
    stage = name,
    analysis = stage$get_analysis(),
    statistics_sha256 = statistics_hash
  )
  if (name %in% .ptm_json_enrichment_methods) {
    artifact$document <- .ptm_enrichment_document(stage)
  } else {
    artifact$result <- .pack_ptm_value(stage$get_results())
  }
  if (identical(name, "KinaseInputs")) {
    artifact$settings <- stage$get_statistics()$get_inputs()$get_parameters()$kinaselib
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".ptm-cbor-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeBin(secretbase::cborenc(artifact), temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not write CBOR artifact: ", path)
  }
  invisible(path)
}

.read_ptm_cbor <- function(path, statistics_hash) {
  artifact <- secretbase::cbordec(readBin(path, what = "raw", n = file.info(path)$size))
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
  container <- statistics$as_container()
  container$uns$prophosqua$stage <- "PTM_results"
  container$uns$prophosqua$completed_enrichments <- enabled
  for (key in names(artifacts)) {
    artifact <- artifacts[[key]]
    modality <- .ptm_analysis_modality(artifact$analysis)
    namespace <- container$modalities[[modality]]$uns$prophosqua
    if (artifact$stage %in% .ptm_json_enrichment_methods) {
      namespace$enrichment_documents[[key]] <- artifact$document
    } else {
      namespace$completed_stages[[key]] <- .pack_ptm_value(.unpack_ptm_value(artifact$result))
    }
    container$modalities[[modality]]$uns$prophosqua <- namespace
  }
  .write_ptm_container(container, output_h5mu)
  invisible(read_ptm_h5mu(output_h5mu, PTM_results))
}
