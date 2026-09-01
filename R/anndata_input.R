#' Read a Pair of prolfquapp DEA AnnData Files
#'
#' Reads and validates the site-level and total-proteome `AnnData.h5ad` files
#' produced by prolfquapp. The returned records contain ordinary R data frames
#' and prolfqua configuration objects; downstream statistics do not depend on
#' the AnnData implementation.
#'
#' @param site_h5ad Path to the site-level prolfquapp DEA H5AD file.
#' @param protein_h5ad Path to the total-proteome prolfquapp DEA H5AD file.
#' @return A list with validated `site` and `protein` experiment records.
#' @export
read_ptm_anndata_pair <- function(site_h5ad, protein_h5ad) {
  site <- .read_site_dea_h5ad(site_h5ad)
  protein <- .read_protein_dea_h5ad(protein_h5ad)
  protein <- .align_protein_samples(site, protein)
  .validate_shared_design(site, protein)

  structure(
    list(site = site, protein = protein),
    class = c("prophosqua_anndata_pair", "list")
  )
}

.read_site_dea_h5ad <- function(path) {
  experiment <- .read_dea_h5ad(path)
  .validate_site_experiment(experiment)
  experiment$site_info <- .site_info_from_var(experiment$var, path)
  experiment
}

.read_protein_dea_h5ad <- function(path) {
  experiment <- .read_dea_h5ad(path)
  .validate_protein_experiment(experiment)
  experiment
}

.read_dea_h5ad <- function(path) {
  if (!file.exists(path)) {
    stop("AnnData file not found: ", path, call. = FALSE)
  }

  invisible(rhdf5::H5get_libversion())
  adata <- tryCatch(
    anndataR::read_h5ad(path),
    error = function(error) {
      stop(
        "Cannot read prolfquapp AnnData file '",
        path,
        "': ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  namespace <- .prolfquapp_namespace(adata, path)
  .validate_dea_anndata(adata, namespace, path)

  obs <- as.data.frame(adata$obs)
  var <- as.data.frame(adata$var)
  config <- prolfqua::list_to_AnalysisConfiguration(
    namespace$analysis_configuration
  )
  transformed <- as.matrix(adata$layers[["transformed"]])

  structure(
    list(
      source_path = normalizePath(path, mustWork = TRUE),
      schema_version = namespace$schema_version,
      sample_key = namespace$sample_key,
      feature_keys = as.character(namespace$feature_keys),
      obs = obs,
      var = var,
      configuration = config,
      normalized_abundances = .anndata_long_table(
        adata,
        obs,
        var,
        transformed,
        config
      ),
      differential_results = .anndata_dea_results(
        adata,
        namespace,
        var
      )
    ),
    class = c("prophosqua_dea_experiment", "list")
  )
}

.prolfquapp_namespace <- function(adata, path) {
  namespace <- adata$uns[["prolfquapp"]]
  if (!is.list(namespace)) {
    stop(
      "AnnData file has no prolfquapp metadata namespace: ",
      path,
      call. = FALSE
    )
  }
  namespace
}

.validate_dea_anndata <- function(adata, namespace, path) {
  required_metadata <- c(
    "analysis_configuration",
    "artifact_type",
    "feature_keys",
    "sample_key",
    "schema_version",
    "varm_columns"
  )
  missing_metadata <- setdiff(required_metadata, names(namespace))
  if (length(missing_metadata) > 0L) {
    stop(
      "AnnData prolfquapp metadata is missing: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }
  if (!identical(namespace$artifact_type, "dea_results")) {
    stop(
      "Expected a prolfquapp dea_results artifact in ",
      path,
      "; found '",
      namespace$artifact_type,
      "'.",
      call. = FALSE
    )
  }
  if (!identical(namespace$schema_version, "1.0.0")) {
    stop(
      "Unsupported prolfquapp DEA-results schema '",
      namespace$schema_version,
      "' in ",
      path,
      "; supported schema is 1.0.0.",
      call. = FALSE
    )
  }

  obs <- as.data.frame(adata$obs)
  var <- as.data.frame(adata$var)
  .validate_axis(rownames(obs), "sample", path)
  .validate_axis(rownames(var), "feature", path)
  .validate_sample_key(obs, namespace$sample_key, path)

  required_layers <- c("raw", "transformed")
  missing_layers <- setdiff(required_layers, adata$layers_keys())
  if (length(missing_layers) > 0L) {
    stop(
      "AnnData file is missing required layer(s): ",
      paste(missing_layers, collapse = ", "),
      call. = FALSE
    )
  }
  expected_shape <- c(nrow(obs), nrow(var))
  for (layer_name in adata$layers_keys()) {
    layer_shape <- dim(adata$layers[[layer_name]])
    if (!identical(as.integer(layer_shape), as.integer(expected_shape))) {
      stop(
        "AnnData layer '",
        layer_name,
        "' is not aligned to the sample and feature axes.",
        call. = FALSE
      )
    }
  }

  dea_keys <- grep("^dea__", adata$varm_keys(), value = TRUE)
  if (length(dea_keys) == 0L) {
    stop("AnnData file contains no DEA result matrices: ", path, call. = FALSE)
  }
  for (key in dea_keys) {
    .validate_dea_varm(adata, namespace, key, nrow(var))
  }
}

.validate_axis <- function(axis_names, label, path) {
  if (
    is.null(axis_names) ||
      anyNA(axis_names) ||
      any(!nzchar(axis_names)) ||
      anyDuplicated(axis_names)
  ) {
    stop(
      "AnnData ",
      label,
      " names must be present, non-empty, and unique in ",
      path,
      ".",
      call. = FALSE
    )
  }
}

.validate_sample_key <- function(obs, sample_key, path) {
  if (!is.character(sample_key) || length(sample_key) != 1L || !nzchar(sample_key)) {
    stop("AnnData sample_key must be one non-empty column name in ", path, ".", call. = FALSE)
  }
  if (!sample_key %in% names(obs)) {
    stop("AnnData sample_key column is absent from obs: ", sample_key, call. = FALSE)
  }

  sample_ids <- as.character(obs[[sample_key]])
  .validate_axis(sample_ids, "sample-key", path)
  if (!identical(sample_ids, rownames(obs))) {
    stop("AnnData sample_key values must equal obs names in the same order.", call. = FALSE)
  }
}

.validate_dea_varm <- function(adata, namespace, key, n_features) {
  matrix <- as.matrix(adata$varm[[key]])
  columns <- namespace$varm_columns[[key]]
  if (is.null(columns)) {
    stop("AnnData metadata has no column names for varm '", key, "'.", call. = FALSE)
  }
  if (nrow(matrix) != n_features || ncol(matrix) != length(columns)) {
    stop("AnnData varm '", key, "' is not aligned to its declared columns.", call. = FALSE)
  }
}

.validate_site_experiment <- function(experiment) {
  required_keys <- c("protein_Id", "site")
  missing_keys <- setdiff(required_keys, experiment$feature_keys)
  if (length(missing_keys) > 0L) {
    stop(
      "The site AnnData is missing feature role(s): ",
      paste(missing_keys, collapse = ", "),
      ". The site and protein inputs may be swapped.",
      call. = FALSE
    )
  }
  .require_columns(
    experiment$var,
    c(
      "protein_Id",
      "site",
      "posInProtein",
      "modAA",
      "SequenceWindow",
      "description",
      "protein_length"
    ),
    "site AnnData var"
  )
  if (anyNA(experiment$var$protein_Id) || any(!nzchar(experiment$var$protein_Id))) {
    stop("Every site feature must declare a non-empty protein_Id.", call. = FALSE)
  }
}

.validate_protein_experiment <- function(experiment) {
  if (!"protein_Id" %in% experiment$feature_keys || "site" %in% experiment$feature_keys) {
    stop(
      "The protein AnnData must declare protein_Id, but not site, as a feature role. ",
      "The site and protein inputs may be swapped.",
      call. = FALSE
    )
  }
  .require_columns(
    experiment$var,
    c("protein_Id", "description", "protein_length"),
    "protein AnnData var"
  )
}

.require_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

.anndata_long_table <- function(adata, obs, var, transformed, config) {
  duplicate_columns <- intersect(names(obs), names(var))
  if (length(duplicate_columns) > 0L) {
    stop(
      "AnnData obs and var columns overlap: ",
      paste(duplicate_columns, collapse = ", "),
      call. = FALSE
    )
  }

  n_obs <- nrow(obs)
  n_vars <- nrow(var)
  long <- cbind(
    obs[rep(seq_len(n_obs), each = n_vars), , drop = FALSE],
    var[rep(seq_len(n_vars), times = n_obs), , drop = FALSE]
  )
  rownames(long) <- NULL
  long$normalized_abundance <- as.vector(t(transformed))

  if ("nr_children" %in% adata$layers_keys()) {
    nr_children_name <- config$nr_children
    if (
      !is.character(nr_children_name) ||
        length(nr_children_name) != 1L ||
        !nzchar(nr_children_name)
    ) {
      stop("The AnnData nr_children layer has no configured column name.", call. = FALSE)
    }
    long[[nr_children_name]] <- as.vector(
      t(as.matrix(adata$layers[["nr_children"]]))
    )
  }
  long
}

.anndata_dea_results <- function(adata, namespace, var) {
  dea_keys <- grep("^dea__", adata$varm_keys(), value = TRUE)
  tables <- lapply(
    dea_keys,
    function(key) .anndata_dea_table(adata, namespace, var, key)
  )
  results <- dplyr::bind_rows(tables)

  contrast_order <- as.character(namespace$contrasts$contrast_name)
  if ("contrast" %in% names(results) && length(contrast_order) > 0L) {
    results <- results[
      order(match(results$contrast, contrast_order), results$.feature_order),
      ,
      drop = FALSE
    ]
  }
  results$.feature_order <- NULL
  rownames(results) <- NULL
  results
}

.anndata_dea_table <- function(adata, namespace, var, key) {
  matrix <- as.matrix(adata$varm[[key]])
  colnames(matrix) <- as.character(namespace$varm_columns[[key]])
  numeric_results <- as.data.frame(matrix)
  annotations <- namespace$varm_annotations[[key]]
  if (is.null(annotations) || length(annotations) == 0L) {
    annotations <- data.frame(row.names = seq_len(nrow(var)))
  } else {
    annotations <- as.data.frame(annotations, stringsAsFactors = FALSE)
  }
  if (nrow(annotations) != nrow(var)) {
    stop("AnnData annotations for varm '", key, "' are not feature-aligned.", call. = FALSE)
  }

  duplicate_columns <- intersect(
    names(var),
    c(names(numeric_results), names(annotations))
  )
  if (length(duplicate_columns) > 0L) {
    stop(
      "AnnData DEA result repeats feature column(s): ",
      paste(duplicate_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(intersect(names(numeric_results), names(annotations))) > 0L) {
    stop("AnnData DEA numeric and annotation columns overlap for varm '", key, "'.", call. = FALSE)
  }

  data.frame(
    var,
    numeric_results,
    annotations,
    .feature_order = seq_len(nrow(var)),
    check.names = FALSE
  )
}

.site_info_from_var <- function(var, path) {
  columns <- c(
    "site",
    "posInProtein",
    "modAA",
    "SequenceWindow",
    "protein_Id",
    "gene_name",
    "protein_length"
  )
  columns <- intersect(columns, names(var))
  site_info <- unique(var[, columns, drop = FALSE])
  if (anyDuplicated(site_info$site)) {
    stop("PTM site metadata is not unique by site in: ", path, call. = FALSE)
  }
  site_info
}

.align_protein_samples <- function(site, protein) {
  site_ids <- as.character(site$obs[[site$sample_key]])
  protein_ids <- as.character(protein$obs[[protein$sample_key]])
  missing_from_protein <- setdiff(site_ids, protein_ids)
  missing_from_site <- setdiff(protein_ids, site_ids)
  if (length(missing_from_protein) > 0L || length(missing_from_site) > 0L) {
    stop(
      "Site and protein AnnData sample sets differ; missing from protein: ",
      .collapse_or_none(missing_from_protein),
      "; missing from site: ",
      .collapse_or_none(missing_from_site),
      ".",
      call. = FALSE
    )
  }

  order_index <- match(site_ids, protein_ids)
  protein$obs <- protein$obs[order_index, , drop = FALSE]
  protein$normalized_abundances <- protein$normalized_abundances[
    order(
      match(
        protein$normalized_abundances[[protein$sample_key]],
        site_ids
      )
    ),
    ,
    drop = FALSE
  ]
  rownames(protein$normalized_abundances) <- NULL
  protein
}

.validate_shared_design <- function(site, protein) {
  site_factors <- names(site$configuration$factors)
  protein_factors <- names(protein$configuration$factors)
  if (!setequal(site_factors, protein_factors)) {
    stop(
      "Site and protein AnnData files declare different design factors.",
      call. = FALSE
    )
  }

  for (factor_name in site_factors) {
    .require_columns(site$obs, factor_name, "site AnnData obs")
    .require_columns(protein$obs, factor_name, "protein AnnData obs")
    if (
      !identical(
        as.character(site$obs[[factor_name]]),
        as.character(protein$obs[[factor_name]])
      )
    ) {
      stop(
        "Site and protein AnnData disagree on design factor '",
        factor_name,
        "'.",
        call. = FALSE
      )
    }
  }
}

.collapse_or_none <- function(values) {
  if (length(values) == 0L) "none" else paste(values, collapse = ", ")
}
