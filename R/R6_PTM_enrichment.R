#' A Completed PTM Enrichment Stage
#'
#' One step of the enrichment of one analysis, DPA, DPU or CF, computed from the
#' stage before it or restored from its stored result. The steps chain as
#' PTM_statistics -> [PTMSEA], and PTM_statistics -> [KinaseInputs] ->
#' [KinaseAssignments] -> [KinaseGSEA], with KinaseAssignments ->
#' [MotifEnrichment] -> [MEA]. KinaseAssignments and MotifEnrichment come from
#' the kinase-library tool and are only restored.
#' @export
PTM_enrichment <- R6::R6Class(
  "PTM_enrichment",
  private = list(source = NULL, analysis = NULL, result = NULL),
  public = list(
    #' @description Validate and keep a completed stage.
    #' @param source The stage this one is computed from.
    #' @param analysis DPA, DPU or CF.
    #' @param result The completed result.
    initialize = function(source, analysis, result) {
      stage <- class(self)[[1L]]
      spec <- .ENRICHMENT_STAGES[[stage]]
      if (!inherits(source, spec$source)) {
        stop(stage, " requires ", spec$source)
      }
      .validate_enrichment_analysis(analysis)
      if (inherits(source, "PTM_enrichment") && !identical(analysis, source$get_analysis())) {
        stop("Source analysis differs from target analysis.")
      }
      .require_ptm_fields(result, spec$fields, stage)
      private$source <- source
      private$analysis <- analysis
      # Kept packed, so a computed stage returns exactly what a restored one does.
      private$result <- .pack_ptm_value(result)
    },
    #' @description Return the stage this one was computed from.
    get_source = function() private$source,
    #' @description Return the complete statistics component.
    get_statistics = function() private$source$get_statistics(),
    #' @description Return the analysis name.
    get_analysis = function() private$analysis,
    #' @description Return the complete result.
    get_results = function() .unpack_ptm_value(private$result)
  )
)

#' PTM-SEA of One Analysis against PTMsigDB
#' @export
PTMSEA <- R6::R6Class(
  "PTMSEA",
  inherit = PTM_enrichment,
  public = list(
    #' @description Compute PTM-SEA of one analysis, or restore it.
    #' @param source Complete PTM_statistics stage.
    #' @param analysis DPA, DPU or CF.
    #' @param result Stored result, when restoring.
    initialize = function(source, analysis, result = .compute_ptmsea_stage(source, analysis)) {
      super$initialize(source, analysis, result)
    }
  )
)

#' Kinase-Library Inputs of One Analysis
#'
#' The sequence windows the kinase-library tool scans and the ranked windows of
#' each contrast its motif enrichment walks.
#' @export
KinaseInputs <- R6::R6Class(
  "KinaseInputs",
  inherit = PTM_enrichment,
  public = list(
    #' @description Compute the kinase-library inputs of one analysis, or restore them.
    #' @param source Complete PTM_statistics stage.
    #' @param analysis DPA, DPU or CF.
    #' @param result Stored result, when restoring.
    initialize = function(source, analysis, result = .compute_kinaseinputs_stage(source, analysis)) {
      super$initialize(source, analysis, result)
    }
  )
)

#' Kinase-Library Motif-Scan Assignments of One Analysis
#' @export
KinaseAssignments <- R6::R6Class("KinaseAssignments", inherit = PTM_enrichment)

#' Kinase GSEA of One Analysis against the Motif-Scan Assignments
#' @export
KinaseGSEA <- R6::R6Class(
  "KinaseGSEA",
  inherit = PTM_enrichment,
  public = list(
    #' @description Compute the kinase GSEA of one analysis, or restore it.
    #' @param source Complete KinaseAssignments stage.
    #' @param analysis DPA, DPU or CF.
    #' @param result Stored result, when restoring.
    initialize = function(source, analysis, result = .compute_kinasegsea_stage(source, analysis)) {
      super$initialize(source, analysis, result)
    }
  )
)

#' Kinase-Library Motif Enrichment of One Analysis
#' @export
MotifEnrichment <- R6::R6Class("MotifEnrichment", inherit = PTM_enrichment)

#' Motif Enrichment Analysis Tables of One Analysis
#' @export
MEA <- R6::R6Class(
  "MEA",
  inherit = PTM_enrichment,
  public = list(
    #' @description Tabulate the motif enrichment of one analysis, or restore it.
    #' @param source Complete MotifEnrichment stage.
    #' @param analysis DPA, DPU or CF.
    #' @param result Stored result, when restoring.
    initialize = function(source, analysis, result = .compute_mea_tables(source$get_results()$mea_results)) {
      super$initialize(source, analysis, result)
    }
  )
)

# Each stage: the stage it is computed from and the result fields it keeps.
.ENRICHMENT_STAGES <- list(
  PTMSEA = list(source = "PTM_statistics", fields = c("results", "all_clean")),
  KinaseInputs = list(source = "PTM_statistics", fields = c("seqwindows", "ranks")),
  KinaseAssignments = list(source = "KinaseInputs", fields = "term2gene"),
  MotifEnrichment = list(source = "KinaseAssignments", fields = c("mea_results", "gsea_json")),
  KinaseGSEA = list(source = "KinaseAssignments", fields = c("gsea_results", "all_results", "gsea_info")),
  MEA = list(source = "MotifEnrichment", fields = c("mea_clean", "summary_df"))
)

# The completed enrichments, stored as string_gsea documents.
.ptm_json_enrichment_methods <- c("PTMSEA", "KinaseGSEA", "MEA")

.ptm_analysis_modality <- function(analysis) {
  c(DPA = "enriched", DPU = "enriched", CF = "enriched_CF")[[analysis]]
}

.validate_enrichment_analysis <- function(analysis) {
  if (!analysis %in% c("DPA", "DPU", "CF")) stop("Unknown PTM analysis: ", analysis)
}

.compute_ptmsea_stage <- function(source, analysis) {
  inputs <- source$get_inputs()
  parameters <- inputs$get_parameters()
  .compute_ptmsea_tables(
    source$get_tables()[[analysis]],
    inputs$get_resources()$ptmsigdb,
    parameters$analyses[[tolower(analysis)]]$stat_column,
    parameters$ptmsigdb$trim_to,
    parameters$gsea$min_size,
    parameters$gsea$max_size,
    parameters$gsea$n_perm
  )
}

.compute_kinaseinputs_stage <- function(source, analysis) {
  data <- filter_sequence_windows(source$get_tables()[[analysis]])
  stat_column <- source$get_inputs()$get_parameters()$analyses[[tolower(analysis)]]$stat_column
  contrasts <- unique(data$contrast)
  list(
    seqwindows = dplyr::arrange(dplyr::distinct(dplyr::select(data, "SequenceWindow")), .data$SequenceWindow),
    ranks = stats::setNames(
      lapply(contrasts, function(contrast) rank_sites_for_mea(data, stat_column, contrast)),
      contrasts
    )
  )
}

.compute_kinasegsea_stage <- function(source, analysis) {
  statistics <- source$get_statistics()
  parameters <- statistics$get_inputs()$get_parameters()
  .compute_kinase_tables(
    statistics$get_tables()[[analysis]],
    source$get_results()$term2gene,
    parameters$analyses[[tolower(analysis)]]$stat_column,
    parameters$gsea$min_size,
    parameters$kinaselib$gsea_max_size,
    parameters$gsea$n_perm
  )
}
