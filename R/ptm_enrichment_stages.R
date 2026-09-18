.ptm_analysis_modality <- function(analysis) {
  c(DPA = "enriched", DPU = "cf", CF = "cf")[[analysis]]
}

.compute_ptmsea_stage <- function(source, analysis) {
  inputs <- source$get_inputs()
  parameters <- inputs$get_parameters()
  .compute_ptmsea_tables(
    source$get_tables()[[analysis]],
    inputs$get_resources()$ptmsigdb,
    analysis,
    parameters$analyses[[tolower(analysis)]]$stat_column,
    parameters$ptmsigdb$trim_to,
    parameters$gsea$min_size,
    parameters$gsea$max_size,
    parameters$gsea$n_perm
  )
}

.compute_kinaseinputs_stage <- function(source, analysis) {
  data <- filter_sequence_windows(canonicalize_sequence_window(source$get_tables()[[analysis]]))
  stat_column <- source$get_inputs()$get_parameters()$analyses[[tolower(analysis)]]$stat_column
  seqwindows <- dplyr::arrange(dplyr::distinct(dplyr::select(data, "SequenceWindow")), .data$SequenceWindow)
  ranks <- lapply(unique(data$contrast), function(contrast) rank_sites_for_mea(data, stat_column, contrast))
  names(ranks) <- unique(data$contrast)
  list(seqwindows = seqwindows, ranks = ranks)
}

.compute_kinaseassignments_stage <- function(source, analysis) {
  stop("KinaseAssignments requires completed kinase-library scoring results.")
}

.compute_kinasegsea_stage <- function(source, analysis) {
  statistics <- source$get_statistics()
  parameters <- statistics$get_inputs()$get_parameters()$gsea
  .compute_kinase_tables(
    statistics$get_tables()[[analysis]],
    source$get_results()$term2gene,
    analysis,
    parameters$min_size,
    parameters$max_size,
    parameters$n_perm
  )
}

.compute_motifenrichment_stage <- function(source, analysis) {
  stop("MotifEnrichment requires completed kinase-library MEA results.")
}

.compute_mea_stage <- function(source, analysis) {
  result <- source$get_results()$mea_results
  # Preserve the column names formerly produced by the R CSV reader.
  names(result) <- make.names(names(result), unique = TRUE)
  .compute_mea_tables(result, length(unique(result$contrast)))
}

.load_completed_enrichment <- function(container, Type, load_source) {
  analysis <- container$uns$prophosqua$analysis
  modality <- .ptm_analysis_modality(analysis)
  key <- .ptm_varm_key(Type$classname, analysis)
  if (Type$classname %in% c("PTMSEA", "KinaseGSEA", "MEA")) {
    documents <- container$modalities[[modality]]$uns$prophosqua$enrichment_documents
    .require_ptm_fields(documents, key, Type$classname)
    result <- .restore_ptm_enrichment_result(documents[[key]], key) # nolint: object_usage_linter.
    return(Type$new(load_source(container), analysis, result))
  }
  stages <- container$modalities[[modality]]$uns$prophosqua$completed_stages
  .require_ptm_fields(stages, key, Type$classname)
  Type$new(load_source(container), analysis, .unpack_ptm_value(stages[[key]]))
}

.load_ptmsea <- function(container) .load_completed_enrichment(container, PTMSEA, .load_ptm_statistics)

.load_kinase_inputs <- function(container) .load_completed_enrichment(container, KinaseInputs, .load_ptm_statistics)

.load_kinase_assignments <- function(container) {
  .load_completed_enrichment(container, KinaseAssignments, .load_kinase_inputs)
}

.load_kinase_gsea <- function(container) .load_completed_enrichment(container, KinaseGSEA, .load_kinase_assignments)

.load_motif_enrichment <- function(container) {
  .load_completed_enrichment(container, MotifEnrichment, .load_kinase_assignments)
}

.load_mea <- function(container) .load_completed_enrichment(container, MEA, .load_motif_enrichment)

.validate_enrichment_analysis <- function(analysis) {
  if (!analysis %in% c("DPA", "DPU", "CF")) stop("Unknown PTM analysis: ", analysis)
}

.enrichment_container <- function(stage) {
  container <- stage$get_source()$as_container()
  classname <- class(stage)[1L]
  analysis <- stage$get_analysis()
  container$uns$prophosqua$stage <- classname
  container$uns$prophosqua$analysis <- analysis
  modality <- .ptm_analysis_modality(analysis)
  namespace <- container$modalities[[modality]]$uns$prophosqua
  key <- .ptm_varm_key(classname, analysis)
  if (classname %in% c("PTMSEA", "KinaseGSEA", "MEA")) {
    namespace$enrichment_documents[[key]] <- .ptm_enrichment_document(stage)
  } else {
    namespace$completed_stages[[key]] <- .pack_ptm_value(stage$get_results())
  }
  container$modalities[[modality]]$uns$prophosqua <- namespace
  container
}
