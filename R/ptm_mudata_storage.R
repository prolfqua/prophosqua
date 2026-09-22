.cf_container <- function(inputs, result) {
  container <- inputs$as_container()
  container$uns$prophosqua$stage <- "CF"
  wide <- result$wide_data
  ids <- as.character(wide[[site_column(wide)]])
  var <- inputs$get_enriched()$var
  var <- as.data.frame(var)[match(ids, var$site), , drop = FALSE]
  samples <- rownames(container$obs)
  values <- t(as.matrix(wide[, samples, drop = FALSE]))
  dimnames(values) <- list(samples, rownames(var))
  cf <- anndataR::AnnData(X = values, obs = container$obs, var = var)
  cf$uns$prophosqua <- list(schema_version = "2.0.0", report_data = .pack_cf_result(result))
  cf <- .add_ptm_method(cf, result$results, "correct_first", names(result$contrasts))
  container$modalities$cf <- cf
  container
}

.add_ptm_method <- function(adata, table, method, contrasts) {
  site_col <- site_column(table)
  # Outer-join context is retained in the complete report tables. Only rows
  # belonging to this modality are feature-aligned in varm.
  rows <- !is.na(table[[site_col]]) & table[[site_col]] %in% adata$var$site
  payload <- .ptm_analysis_payload(table[rows, , drop = FALSE], adata$var, method, contrasts)
  namespace <- adata$uns$prophosqua
  namespace$schema_version <- "2.0.0"
  namespace$result_keys[[method]] <- payload$keys
  for (key in names(payload$values)) {
    adata$varm[[key]] <- payload$values[[key]]
    adata$varm[[paste0(key, "__present")]] <- matrix(payload$present[[key]], ncol = 1L)
    namespace$varm_columns[[key]] <- payload$columns[[key]]
    namespace$varm_annotations[[key]] <- payload$annotations[[key]]
  }
  adata$uns$prophosqua <- namespace
  adata
}

.statistics_container <- function(statistics) {
  inputs <- statistics$get_inputs()
  result <- statistics$get_dpa_dpu()
  container <- .cf_container(inputs, statistics$get_cf())
  container$uns$prophosqua$stage <- "PTM_statistics"
  enriched <- container$modalities$enriched
  enriched$uns$prophosqua <- list(dpa_dpu = .pack_ptm_value(result))
  container$modalities$enriched <- .add_ptm_method(
    enriched,
    result$combined_site_prot,
    "dpa",
    names(inputs$get_contrasts())
  )
  cf <- container$modalities$cf
  cf <- .add_ptm_method(cf, result$combined_test_diff, "dpu", names(inputs$get_contrasts()))
  cf <- .add_ptm_method(cf, result$combined_test_diff_unmoderated, "dpu_unmoderated", names(inputs$get_contrasts()))
  container$modalities$cf <- cf
  container
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
  readers <- list(
    DEA_enriched_total = .load_ptm_inputs,
    CF = .load_ptm_cf,
    DPA_DPU = .load_ptm_dpa_dpu,
    PTM_statistics = .load_ptm_statistics,
    PTMSEA = .load_ptmsea,
    KinaseInputs = .load_kinase_inputs,
    KinaseAssignments = .load_kinase_assignments,
    KinaseGSEA = .load_kinase_gsea,
    MotifEnrichment = .load_motif_enrichment,
    MEA = .load_mea,
    PTM_results = function(container) .load_ptm_results(container, path)
  )
  reader <- readers[[metadata$stage]]
  if (is.null(reader)) {
    stop("Unknown PTM stage: ", metadata$stage)
  }
  result <- reader(container)
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
  for (key in grep("^(dpa|dpu|correct_first)__", enriched$varm_keys(), value = TRUE)) {
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
  .require_ptm_fields(container$modalities, "cf", "CF modalities")
  namespace <- container$modalities$cf$uns$prophosqua
  .require_ptm_fields(namespace, c("report_data", "result_keys"), "CF metadata")
  .validate_ptm_method(container$modalities$cf, "correct_first")
  CF$new(.load_ptm_inputs(container), .unpack_cf_result(namespace$report_data))
}

.load_ptm_dpa_dpu <- function(container) {
  namespace <- container$modalities$enriched$uns$prophosqua
  .require_ptm_fields(namespace, "dpa_dpu", "DPA/DPU metadata")
  DPA_DPU$new(.load_ptm_inputs(container), .unpack_ptm_value(namespace$dpa_dpu))
}

.load_ptm_statistics <- function(container) {
  .validate_ptm_method(container$modalities$enriched, "dpa")
  .validate_ptm_method(container$modalities$cf, "dpu")
  .validate_ptm_method(container$modalities$cf, "dpu_unmoderated")
  PTM_statistics$new(.load_ptm_dpa_dpu(container), .load_ptm_cf(container))
}

.validate_ptm_method <- function(adata, method) {
  namespace <- adata$uns$prophosqua
  .require_ptm_fields(namespace$result_keys, method, "PTM result keys")
  keys <- namespace$result_keys[[method]]
  if (!length(keys)) {
    stop("PTM result keys are empty: ", method)
  }
  for (key in keys) {
    mask_key <- paste0(key, "__present")
    .require_ptm_fields(adata$varm, c(key, mask_key), "PTM result matrices")
    .require_ptm_fields(namespace$varm_columns, key, "PTM result columns")
    values <- adata$varm[[key]]
    present <- adata$varm[[mask_key]]
    valid <- c(
      identical(dim(values), c(nrow(adata$var), length(namespace$varm_columns[[key]]))),
      identical(dim(present), c(nrow(adata$var), 1L)),
      is.logical(present),
      !anyNA(present)
    )
    if (!all(valid)) {
      stop("Invalid PTM result matrix or presence mask: ", key)
    }
  }
}

#' Import paired DEA files into the first complete MuData stage
#' @param enriched_h5ad,total_h5ad Producer-owned DEA files.
#' @param output_h5mu Destination.
#' @param resources,parameters Imported reference data and analysis parameters.
#' @return Complete paired-input stage, invisibly.
#' @export
import_ptm_h5mu <- function(enriched_h5ad, total_h5ad, output_h5mu, resources = list(), parameters = list()) {
  paths <- c(enriched = enriched_h5ad, total = total_h5ad)
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(enriched_h5ad),
    anndataR::read_h5ad(total_h5ad),
    resources,
    parameters,
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
  inputs <- read_ptm_h5mu(input_h5mu, DEA_enriched_total)
  dpa_dpu <- inputs$build(DPA_DPU)
  cf <- inputs$build(CF)
  result <- dpa_dpu$build(PTM_statistics, cf = cf)
  result$write_h5mu(output_h5mu)
  invisible(result)
}
