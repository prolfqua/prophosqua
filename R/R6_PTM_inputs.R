#' Complete paired enriched and total DEA experiments
#'
#' The first PTM stage: both DEA experiments as prolfquapp wrote them, with the
#' reference resources and analysis parameters imported beside them. MuData is
#' the persistence boundary of every stage.
#' @importFrom R6 R6Class
#' @export
DEA_enriched_total <- R6::R6Class(
  "DEA_enriched_total",
  private = list(modalities = NULL, pair = NULL, metadata = NULL),
  public = list(
    #' @description Validate and retain both complete DEA experiments.
    #' @param enriched,total AnnData experiments written by prolfquapp.
    #' @param resources Imported reference data.
    #' @param parameters Analysis parameters. `remove_contaminants = TRUE`
    #'   drops the features either DEA flags as contaminants (`CON`); they are
    #'   kept otherwise.
    #' @param provenance Source identities recorded during import.
    #' @param total_peptide Optional peptide-level DEA of the total proteome,
    #'   with the same samples and contrasts. No analysis reads it; every stage
    #'   carries it as the modality `total_peptide`.
    initialize = function(
      enriched,
      total,
      resources = list(),
      parameters = list(),
      provenance = list(),
      total_peptide = NULL
    ) {
      pair <- .ptm_pair(
        prolfquapp::DEAResultReader$new(enriched),
        prolfquapp::DEAResultReader$new(total),
        isTRUE(parameters$remove_contaminants)
      )
      by_name <- function(contrasts) contrasts[order(names(contrasts))]
      if (!identical(by_name(pair$site$contrasts), by_name(pair$protein$contrasts))) {
        stop("Paired DEA contrast definitions differ.")
      }
      if (!is.null(total_peptide)) {
        .validate_total_peptide(total_peptide, pair)
      }
      private$pair <- pair
      private$modalities <- c(
        list(enriched = enriched$clone(deep = TRUE), total = total$clone(deep = TRUE)),
        if (!is.null(total_peptide)) list(total_peptide = total_peptide$clone(deep = TRUE))
      )
      private$metadata <- list(
        resources = .pack_ptm_value(resources),
        parameters = .pack_ptm_value(parameters),
        provenance = .pack_ptm_value(provenance)
      )
    },
    #' @description Return an independent enriched experiment.
    get_enriched = function() private$modalities$enriched$clone(deep = TRUE),
    #' @description Return an independent total experiment.
    get_total = function() private$modalities$total$clone(deep = TRUE),
    #' @description Return shared design in enriched sample order.
    get_design = function() as.data.frame(private$modalities$enriched$obs),
    #' @description Return stored named contrast expressions.
    get_contrasts = function() private$pair$site$contrasts,
    #' @description Return the experiments as prolfquapp's reader decoded them.
    get_pair = function() private$pair,
    #' @description Return imported reference resources.
    get_resources = function() .unpack_ptm_value(private$metadata$resources),
    #' @description Return recorded analysis parameters.
    get_parameters = function() .unpack_ptm_value(private$metadata$parameters),
    #' @description Return original source identities.
    get_provenance = function() .unpack_ptm_value(private$metadata$provenance),
    #' @description Write this complete stage atomically.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(self$as_container(), path),
    #' @description Return a detached storage representation.
    as_container = function() {
      list(
        modalities = lapply(private$modalities, function(experiment) experiment$clone(deep = TRUE)),
        obs = self$get_design(),
        uns = list(prophosqua = c(list(schema_version = "2.0.0", stage = "DEA_enriched_total"), private$metadata))
      )
    }
  )
)

# The peptide-level total DEA belongs to the pair when prolfquapp decodes it
# with the pair's samples and contrasts.
.validate_total_peptide <- function(total_peptide, pair) {
  peptide <- .dea_record(prolfquapp::DEAResultReader$new(total_peptide))
  .validate_paired_samples(pair$site, peptide)
  by_name <- function(contrasts) contrasts[order(names(contrasts))]
  if (!identical(by_name(peptide$contrasts), by_name(pair$site$contrasts))) {
    stop("Contrast definitions of total_peptide differ from the paired DEA.", call. = FALSE)
  }
}

.write_ptm_container <- function(container, path) {
  prolfquapp::write_h5mu(container$modalities, path, container$obs, container$uns)
}
