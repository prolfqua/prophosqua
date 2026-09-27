.CF_VARIANTS <- c("correct_first_protein_imputed", "correct_first_site_protein_imputed")

.statistics_container <- function(statistics) {
  inputs <- statistics$get_inputs()
  container <- inputs$as_container()
  container$uns$prophosqua$stage <- "PTM_statistics"
  container$modalities$enriched_CF <- .enriched_cf_modality(container, statistics$get_cf())
  contrasts <- names(inputs$get_contrasts())
  dpa_dpu <- statistics$get_dpa_dpu()
  enriched <- container$modalities$enriched
  enriched$uns$prophosqua <- list(dpa_dpu = .pack_ptm_value(dpa_dpu[.DPA_DPU_SUMMARIES]))
  for (method in names(.DPA_DPU_METHODS)) {
    enriched <- .add_ptm_method(enriched, dpa_dpu[[.DPA_DPU_METHODS[[method]]]], method, contrasts)
  }
  container$modalities$enriched <- enriched
  container
}

# All of CorrectFirst over the sites it corrects, in the enriched order: X holds
# CF, each imputed variant is a layer named by its result key, and varm holds
# the results of every fitted one. The variants' matrices, samples x site ids,
# already span those sites (.cf_sites()), so the modality takes their axis; a
# site a variant does not correct is NA in its matrix and in its results. uns
# keeps the CF metadata only: everything else is rebuilt from X, obs and var.
.enriched_cf_modality <- function(container, result) {
  obs <- container$obs
  enriched_var <- as.data.frame(container$modalities$enriched$var)
  sites <- unique(unlist(lapply(result$variants, function(variant) colnames(variant$abundances))))
  var <- enriched_var[enriched_var$site %in% sites, , drop = FALSE]
  contrasts <- names(result$contrasts)
  on_axes <- function(abundances) {
    abundances <- abundances[rownames(obs), var$site, drop = FALSE]
    dimnames(abundances) <- list(rownames(obs), rownames(var))
    abundances
  }
  cf <- anndataR::AnnData(
    X = on_axes(.cf_abundances(result$wide_data, rownames(obs), var$site)),
    obs = obs,
    var = var
  )
  metadata <- c(
    result[.CF_METADATA],
    list(ptm_config = prolfqua::R6_extract_values(result$ptm_data$get_config()))
  )
  cf$uns$prophosqua <- list(cf = .pack_ptm_value(metadata))
  cf <- .add_ptm_method(cf, result$results, "correct_first", contrasts)
  for (method in names(result$variants)) {
    variant <- result$variants[[method]]
    if (!is.null(variant$results)) {
      cf <- .add_ptm_method(cf, variant$results, method, contrasts)
    }
    cf$layers[[method]] <- on_axes(variant$abundances)
  }
  cf
}

# The corrected values as the LFQData CF fitted: X on the modality's axes.
.cf_ptm_data <- function(cf, config) {
  values <- as.matrix(cf$X)
  dimnames(values) <- list(cf$obs_names, as.data.frame(cf$var)$site)
  .cf_lfqdata(values, cf$var, cf$obs, prolfqua::list_to_AnalysisConfiguration(config))
}

# Corrected values, samples x site ids, as an LFQData: joined to the site keys
# and the design.
.cf_lfqdata <- function(values, var, obs, config) {
  keys <- dplyr::distinct(dplyr::as_tibble(as.data.frame(var)[config$hierarchy_keys()]))
  long <- tidyr::pivot_longer(
    dplyr::as_tibble(t(values), rownames = "site"),
    cols = -"site",
    names_to = config$sample_name,
    values_to = config$get_response()
  ) |>
    dplyr::inner_join(keys, by = "site") |>
    dplyr::filter(any(!is.na(.data[[config$get_response()]])), .by = "site") |>
    dplyr::left_join(dplyr::as_tibble(as.data.frame(obs)), by = config$sample_name)
  prolfqua::LFQData$new(long, config)
}

.load_cf_variants <- function(cf, contrasts) {
  methods <- intersect(.CF_VARIANTS, cf$layers_keys())
  stats::setNames(
    lapply(methods, function(method) {
      abundances <- as.matrix(cf$layers[[method]])
      dimnames(abundances) <- list(cf$obs_names, as.data.frame(cf$var)$site)
      if (.ptm_key(method, contrasts[[1]]) %in% cf$varm_keys()) {
        return(list(results = .ptm_method_table(cf, method, contrasts), abundances = abundances))
      }
      list(abundances = abundances)
    }),
    methods
  )
}

.add_ptm_method <- function(adata, table, method, contrasts) {
  # Protein-only rows of the outer join carry no site and have no feature to
  # align to; every row that names a site is kept.
  rows <- !is.na(table$site) & table$site %in% adata$var$site
  payload <- .ptm_analysis_payload(table[rows, , drop = FALSE], adata$var, method, contrasts)
  for (key in names(payload$values)) {
    adata$varm[[key]] <- payload$values[[key]]
  }
  namespace <- adata$uns$prophosqua
  namespace$varm_order <- c(namespace$varm_order, payload$order)
  adata$uns$prophosqua <- namespace
  adata
}

#' Read a complete typed PTM stage from MuData
#' @param path H5MU artifact.
#' @param expected Required R6 class, when a caller requires a particular stage.
#' @return A complete stage object; malformed or incomplete artifacts fail.
#' @export
read_ptm_h5mu <- function(path, expected = NULL) {
  container <- prolfquapp::read_h5mu(path)
  metadata <- container$uns$prophosqua
  .require_ptm_fields(metadata, c("stage", "schema_version", "resources", "parameters", "provenance"), "PTM metadata")
  if (!identical(metadata$schema_version, "2.0.0")) {
    stop("Unsupported PTM MuData schema.")
  }
  result <- switch(
    metadata$stage,
    DEA_enriched_total = .load_ptm_inputs(container),
    PTM_statistics = .load_ptm_statistics(container),
    PTM_results = .load_ptm_results(container, path),
    stop("Unknown PTM stage: ", metadata$stage)
  )
  if (!is.null(expected) && !inherits(result, expected$classname)) {
    stop("Expected stage ", expected$classname)
  }
  result
}

.load_ptm_inputs <- function(container) {
  .require_ptm_fields(container$modalities, c("enriched", "total"), "PTM modalities")
  metadata <- container$uns$prophosqua
  enriched <- container$modalities$enriched$clone(deep = TRUE)
  # Derived results do not belong to the paired-input component.
  enriched$uns$prophosqua <- NULL
  derived <- paste0("^(", paste(c(names(.DPA_DPU_METHODS), "correct_first"), collapse = "|"), ")__")
  for (key in grep(derived, enriched$varm_keys(), value = TRUE)) {
    enriched$varm[[key]] <- NULL
  }
  DEA_enriched_total$new(
    enriched,
    container$modalities$total,
    resources = .unpack_ptm_value(metadata$resources),
    parameters = .unpack_ptm_value(metadata$parameters),
    provenance = .unpack_ptm_value(metadata$provenance)
  )
}

.load_ptm_cf <- function(container) {
  .require_ptm_fields(container$modalities, "enriched_CF", "CF modalities")
  cf <- container$modalities$enriched_CF
  .require_ptm_fields(cf$uns$prophosqua, "cf", "CF metadata")
  metadata <- .unpack_ptm_value(cf$uns$prophosqua$cf)
  contrasts <- names(metadata$contrasts)
  ptm_data <- .cf_ptm_data(cf, metadata$ptm_config)
  c(
    metadata[.CF_METADATA],
    list(results = .ptm_method_table(cf, "correct_first", contrasts), ptm_data = ptm_data),
    .cf_wide(ptm_data),
    list(variants = .load_cf_variants(cf, contrasts))
  )
}

.load_ptm_dpa_dpu <- function(container, contrasts) {
  enriched <- container$modalities$enriched
  .require_ptm_fields(enriched$uns$prophosqua, "dpa_dpu", "DPA/DPU metadata")
  tables <- lapply(names(.DPA_DPU_METHODS), function(method) .ptm_method_table(enriched, method, contrasts))
  c(stats::setNames(tables, unname(.DPA_DPU_METHODS)), .unpack_ptm_value(enriched$uns$prophosqua$dpa_dpu))
}

.load_ptm_statistics <- function(container) {
  inputs <- .load_ptm_inputs(container)
  PTM_statistics$new(
    inputs,
    dpa_dpu = .load_ptm_dpa_dpu(container, names(inputs$get_contrasts())),
    cf = .load_ptm_cf(container)
  )
}

#' Import paired DEA files into the first complete MuData stage
#'
#' Reads the two DEA files and the reference resources the enabled analyses
#' need, once; every later stage reads MuData only.
#' @param enriched_h5ad,total_h5ad Producer-owned DEA files.
#' @param output_h5mu Destination.
#' @param parameters Analysis parameters, the merged pipeline configuration.
#' @param ptmsigdb Filtered PTMsigDB as `.rds` or `.gmt`. When `NULL` and the
#'   kinase analyses are enabled, PTMsigDB is downloaded and filtered as
#'   `parameters$ptmsigdb` asks.
#' @return Complete paired-input stage, invisibly.
#' @export
import_ptm_h5mu <- function(enriched_h5ad, total_h5ad, output_h5mu, parameters = list(), ptmsigdb = NULL) {
  paths <- c(enriched = enriched_h5ad, total = total_h5ad)
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(enriched_h5ad),
    anndataR::read_h5ad(total_h5ad),
    resources = .import_ptm_resources(parameters, ptmsigdb),
    parameters = parameters,
    provenance = list(paths = normalizePath(paths), md5 = unname(tools::md5sum(paths)))
  )
  inputs$write_h5mu(output_h5mu)
  invisible(inputs)
}

#' Compute all PTM statistics using only MuData
#' @param input_h5mu Paired-input stage.
#' @param output_h5mu Destination statistics stage.
#' @return Complete statistics stage, invisibly.
#' @export
compute_ptm_results_h5mu <- function(input_h5mu, output_h5mu) {
  result <- PTM_statistics$new(read_ptm_h5mu(input_h5mu, DEA_enriched_total))
  result$write_h5mu(output_h5mu)
  invisible(result)
}
