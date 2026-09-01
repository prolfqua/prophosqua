# Synthetic prolfquapp DEA H5AD files derived from the same legacy fixture the
# existing computation tests use. The converter lives only in tests: production
# code consumes the producer-owned H5AD contract and never reconstructs it from
# XLSX or Parquet.

make_anndata_pair <- function() {
  dirs <- example_dea_pair()
  output_dir <- tempfile("prophosqua_anndata_pair")
  dir.create(output_dir)
  paths <- list(
    site = file.path(output_dir, "site.h5ad"),
    protein = file.path(output_dir, "protein.h5ad"),
    annot_file = dirs$annot_file,
    legacy = dirs
  )
  write_dea_fixture_h5ad(dirs$phospho, paths$site, site = TRUE)
  write_dea_fixture_h5ad(dirs$protein, paths$protein, site = FALSE)
  paths
}

anndata_pair_fixture <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      value <<- make_anndata_pair()
    }
    value
  }
})

ptm_result_fixture <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      paths <- anndata_pair_fixture()
      input_hashes <- tools::md5sum(c(paths$site, paths$protein))
      output <- tempfile(fileext = ".h5ad")
      suppressWarnings(compute_ptm_results_h5ad(
        paths$site,
        paths$protein,
        paths$annot_file,
        output
      ))
      value <<- c(
        paths,
        list(output = output, input_hashes = unname(input_hashes))
      )
    }
    value
  }
})

write_dea_fixture_h5ad <- function(dea_dir, path, site) {
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

  transformed <- dea_fixture_matrix(long, sample_key, feature_key, "normalized_abundance")
  transformed <- transformed[rownames(obs), rownames(var), drop = FALSE]
  nr_children <- dea_fixture_matrix(long, sample_key, feature_key, config$nr_children)
  nr_children <- nr_children[rownames(obs), rownames(var), drop = FALSE]
  varm <- dea_fixture_varm(results, var, feature_key)

  namespace <- list(
    artifact_type = "dea_results",
    schema_version = "1.0.0",
    source_software = "prophosqua-test-fixture",
    analysis_configuration = config,
    layer_names = c("raw", "transformed", "nr_children"),
    feature_keys = feature_keys,
    sample_key = sample_key,
    varm_columns = varm$columns,
    varm_annotations = varm$annotations,
    contrasts = list(contrast_name = unique(results$contrast))
  )
  adata <- anndataR::AnnData(
    X = transformed,
    obs = obs,
    var = var,
    layers = list(
      raw = transformed,
      transformed = transformed,
      nr_children = nr_children
    ),
    varm = varm$values,
    uns = list(prolfquapp = namespace)
  )
  invisible(rhdf5::H5get_libversion())
  adata$write_h5ad(path, compression = "gzip", mode = "w")
  path
}

dea_fixture_matrix <- function(long, sample_key, feature_key, value_column) {
  wide <- long |>
    dplyr::select(
      tidyselect::all_of(c(sample_key, feature_key, value_column))
    ) |>
    tidyr::pivot_wider(
      names_from = tidyselect::all_of(feature_key),
      values_from = tidyselect::all_of(value_column)
    ) |>
    as.data.frame()
  rownames(wide) <- wide[[sample_key]]
  as.matrix(wide[, setdiff(names(wide), sample_key), drop = FALSE])
}

dea_fixture_varm <- function(results, var, feature_key) {
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
      "dea__",
      utils::URLencode(contrast, reserved = TRUE, repeated = TRUE)
    )
    table <- results[results$contrast == contrast, , drop = FALSE]
    table <- table[match(var[[feature_key]], table[[feature_key]]), , drop = FALSE]
    values[[key]] <- as.matrix(table[, numeric_columns, drop = FALSE])
    columns[[key]] <- numeric_columns
    annotations[[key]] <- lapply(
      table[, annotation_columns, drop = FALSE],
      unname
    )
  }
  list(values = values, columns = columns, annotations = annotations)
}

rewrite_fixture_h5ad <- function(path, transform) {
  adata <- anndataR::read_h5ad(path)
  adata <- transform(adata)
  output <- tempfile(fileext = ".h5ad")
  adata$write_h5ad(output, compression = "gzip", mode = "w")
  output
}
