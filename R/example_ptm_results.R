# A deterministic final PTM MuData artifact for package vignettes.

.example_ptm_results_matrix <- function(long, sample_key, feature_key, value_column) {
  wide <- long |>
    dplyr::select(tidyselect::all_of(c(sample_key, feature_key, value_column))) |>
    tidyr::pivot_wider(
      names_from = tidyselect::all_of(feature_key),
      values_from = tidyselect::all_of(value_column)
    ) |>
    as.data.frame()
  rownames(wide) <- wide[[sample_key]]
  as.matrix(wide[, setdiff(names(wide), sample_key), drop = FALSE])
}

.example_ptm_results_varm <- function(results, var, feature_key) {
  numeric_columns <- c(
    "diff",
    "std.error",
    "df",
    "std.error.unmoderated",
    "df.unmoderated",
    "statistic",
    "FDR"
  )
  annotation_columns <- intersect(
    c("contrast", "estimate_type", "modelName"),
    names(results)
  )
  values <- list()
  columns <- list()
  annotations <- list()

  for (contrast in unique(results$contrast)) {
    key <- paste0(
      "constrast_",
      utils::URLencode(contrast, reserved = TRUE, repeated = TRUE)
    )
    table <- results[results$contrast == contrast, , drop = FALSE]
    table <- table[match(var[[feature_key]], table[[feature_key]]), , drop = FALSE]
    values[[key]] <- as.matrix(table[, numeric_columns, drop = FALSE])
    columns[[key]] <- numeric_columns
    annotations[[key]] <- lapply(table[, annotation_columns, drop = FALSE], unname)
  }
  list(values = values, columns = columns, annotations = annotations)
}

.write_example_ptm_dea_h5ad <- function(dea_dir, path, site) {
  long <- arrow::read_parquet(get_dea_parquet(dea_dir))
  config <- yaml::read_yaml(get_dea_yaml(dea_dir))
  results <- suppressMessages(load_and_preprocess_data(
    get_dea_xlsx(dea_dir),
    c("protein_Id", "contrast")
  ))
  sample_key <- config$sample_name
  feature_key <- if (site) "site" else "protein_Id"
  feature_keys <- if (site) c("protein_Id", "site") else "protein_Id"

  obs_columns <- unique(c(
    sample_key,
    config$file_name,
    names(config$factors),
    config$isotope_label
  ))
  obs <- as.data.frame(unique(long[, obs_columns, drop = FALSE]))
  obs <- obs[match(unique(long[[sample_key]]), obs[[sample_key]]), , drop = FALSE]
  rownames(obs) <- obs[[sample_key]]

  annotation_columns <- c(
    "protein_Id",
    "site",
    "description",
    "gene_name",
    "protein_length",
    "posInProtein",
    "modAA",
    "SequenceWindow"
  )
  annotation_columns <- intersect(annotation_columns, names(results))
  var <- as.data.frame(unique(results[, annotation_columns, drop = FALSE]))
  var <- var[!duplicated(var[[feature_key]]), , drop = FALSE]
  feature_ids <- unique(long[[feature_key]])
  var <- var[match(feature_ids, var[[feature_key]]), , drop = FALSE]
  rownames(var) <- feature_ids

  transformed <- .example_ptm_results_matrix(
    long,
    sample_key,
    feature_key,
    "normalized_abundance"
  )
  transformed <- transformed[rownames(obs), rownames(var), drop = FALSE]
  nr_children <- .example_ptm_results_matrix(
    long,
    sample_key,
    feature_key,
    config$nr_children
  )
  nr_children <- nr_children[rownames(obs), rownames(var), drop = FALSE]
  varm <- .example_ptm_results_varm(results, var, feature_key)

  namespace <- list(
    artifact_type = "dea_results",
    schema_version = "2.0.0",
    source_software = "prophosqua-example",
    analysis_configuration = config,
    layer_names = c("rawData", "transformedData", "nr_children"),
    feature_keys = feature_keys,
    sample_key = sample_key,
    varm_columns = varm$columns,
    varm_annotations = varm$annotations,
    contrasts = list(
      contrast_name = unique(results$contrast),
      contrast = unname(derive_contrasts(
        readr::read_tsv(example_dea_pair()$annot_file, show_col_types = FALSE)
      ))
    )
  )
  adata <- anndataR::AnnData(
    X = transformed,
    obs = obs,
    var = var,
    layers = list(
      rawData = transformed,
      transformedData = transformed,
      nr_children = nr_children
    ),
    varm = varm$values,
    uns = list(prolfquapp = namespace)
  )
  invisible(rhdf5::H5get_libversion())
  adata$write_h5ad(path, compression = "gzip", mode = "w")
  path
}

#' Build the Small Final MuData Used by Package Vignettes
#'
#' The artifact follows the same H5AD import, typed R6 build, and H5MU storage
#' path as a pipeline run. Kinase enrichment is disabled because the statistics
#' report needs only the DPA, DPU, and CorrectFirst components.
#'
#' @param path Destination H5MU path.
#' @return Normalized `path`, invisibly.
#' @keywords internal
example_ptm_results_h5mu <- function(path = tempfile(fileext = ".h5mu")) {
  output_dir <- tempfile("prophosqua_ptm_results_example_")
  dir.create(output_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  dirs <- example_dea_pair()
  enriched_h5ad <- file.path(output_dir, "enriched.h5ad")
  total_h5ad <- file.path(output_dir, "total.h5ad")
  input_h5mu <- file.path(output_dir, "inputs.h5mu")
  statistics_h5mu <- file.path(output_dir, "statistics.h5mu")

  .write_example_ptm_dea_h5ad(dirs$phospho, enriched_h5ad, site = TRUE)
  .write_example_ptm_dea_h5ad(dirs$protein, total_h5ad, site = FALSE)
  parameters <- list(
    run_kinase = FALSE,
    analyses = list(
      dpa = list(stat_column = "statistic.site"),
      dpu = list(stat_column = "statistic.site"),
      cf = list(stat_column = "statistic.site")
    )
  )
  import_ptm_h5mu(
    enriched_h5ad,
    total_h5ad,
    input_h5mu,
    parameters = parameters
  )
  statistics <- suppressWarnings(compute_ptm_results_h5mu(
    input_h5mu,
    statistics_h5mu
  ))
  result <- PTM_results$new(statistics, enrichments = list())
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  result$write_h5mu(path)
  invisible(normalizePath(path, mustWork = TRUE))
}
