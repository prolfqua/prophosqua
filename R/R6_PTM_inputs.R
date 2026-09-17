#' Complete paired enriched and total DEA experiments
#'
#' Each stage exposes only results that exist at construction. `build(Type, ...)`
#' constructs a new complete stage without modifying its source. MuData is the
#' persistence boundary for every stage.
#' @importFrom R6 R6Class
#' @export
DEA_enriched_total <- R6::R6Class(
  "DEA_enriched_total",
  private = list(modalities = NULL, metadata = NULL),
  public = list(
    #' @description Validate and retain both complete DEA experiments.
    #' @param enriched,total AnnData experiments written by prolfquapp.
    #' @param resources Imported reference data.
    #' @param parameters Analysis parameters.
    #' @param provenance Source identities recorded during import.
    initialize = function(enriched, total, resources = list(), parameters = list(), provenance = list()) {
      .ptm_pair_from_modalities(enriched, total)
      contrasts <- .dea_contrasts(enriched)
      if (!identical(contrasts, .dea_contrasts(total))) {
        stop("Paired DEA contrast definitions differ.")
      }
      private$modalities <- list(enriched = enriched$clone(deep = TRUE), total = total$clone(deep = TRUE))
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
    get_contrasts = function() .dea_contrasts(private$modalities$enriched),
    #' @description Return decoded experiments for the existing computations.
    get_pair = function() .ptm_pair_from_modalities(self$get_enriched(), self$get_total()),
    #' @description Return imported reference resources.
    get_resources = function() .unpack_ptm_value(private$metadata$resources),
    #' @description Return recorded analysis parameters.
    get_parameters = function() .unpack_ptm_value(private$metadata$parameters),
    #' @description Return original source identities.
    get_provenance = function() .unpack_ptm_value(private$metadata$provenance),
    #' @description Build a complete next stage.
    #' @param Type Target R6 class.
    #' @param ... Arguments required by the target stage.
    build = function(Type, ...) Type$new(self, ...),
    #' @description Write this complete stage atomically.
    #' @param path Destination H5MU file.
    write_h5mu = function(path) .write_ptm_container(self$as_container(), path),
    #' @description Return a detached storage representation.
    as_container = function() {
      list(
        modalities = list(enriched = self$get_enriched(), total = self$get_total()),
        obs = self$get_design(),
        uns = list(prophosqua = c(list(schema_version = "2.0.0", stage = "DEA_enriched_total"), private$metadata))
      )
    }
  )
)

.dea_contrasts <- function(adata) {
  value <- adata$uns$prolfquapp$contrasts
  .require_columns(value, c("contrast_name", "contrast"), "Stored DEA contrasts")
  result <- stats::setNames(as.character(value$contrast), as.character(value$contrast_name))
  if (!length(result) || anyNA(result) || anyDuplicated(names(result))) {
    stop("Invalid stored DEA contrasts.")
  }
  result
}

.ptm_pair_from_modalities <- function(enriched, total) {
  site <- .dea_anndata_record(enriched, "enriched")
  protein <- .dea_anndata_record(total, "total")
  .validate_site_experiment(site)
  .validate_protein_experiment(protein)
  site$site_info <- .site_info_from_var(site$var, "enriched")
  protein <- .align_protein_samples(site, protein)
  .validate_shared_design(site, protein)
  list(site = site, protein = protein)
}

.write_ptm_container <- function(container, path) {
  prolfquapp::write_h5mu(container$modalities, path, container$obs, container$uns)
}
