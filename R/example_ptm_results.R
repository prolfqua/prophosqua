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

.write_example_ptm_dea_h5ad <- function(dea_dir, path, site, annot_file) {
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
        readr::read_tsv(annot_file, show_col_types = FALSE)
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

# Build the richer paired DEA input used only by the statistics vignette.
.example_ptm_results_dea_pair <- function(root) {
  samples <- expand.grid(
    replicate = seq_len(3L),
    G_ = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
  samples$Name <- paste0(samples$G_, samples$replicate)
  samples$raw_file <- paste0(samples$Name, ".raw")
  samples$Control <- ifelse(samples$G_ == "b", "C", "T")

  proteins <- sprintf("P%03d", seq_len(24L))
  site_map <- data.frame(
    protein_Id = rep(proteins, each = 3L),
    posInProtein = rep(c(10L, 20L, 30L), length(proteins)),
    modAA = rep(c("S", "T", "Y"), length(proteins)),
    stringsAsFactors = FALSE
  )
  site_map$site <- paste0(
    site_map$protein_Id,
    "~",
    site_map$modAA,
    site_map$posInProtein
  )

  protein_a <- rep(c(-0.7, -0.3, 0.3, 0.7, 0.1, -0.1), length.out = length(proteins))
  protein_c <- rep(c(0.6, -0.6, 0.2, -0.2, 0.4, -0.4), length.out = length(proteins))
  usage_a <- rep(c(1.4, -1.3, 0.15, -0.15), length.out = nrow(site_map))
  usage_c <- rep(c(-1.25, 0.15, 1.35, -0.1), length.out = nrow(site_map))
  protein_index <- match(site_map$protein_Id, proteins)
  site_a <- protein_a[protein_index] + usage_a
  site_c <- protein_c[protein_index] + usage_c

  amino_acids <- strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1]]
  site_map$SequenceWindow <- vapply(
    seq_len(nrow(site_map)),
    function(index) {
      sequence <- amino_acids[
        ((seq_len(15L) + index * 3L - 2L) %% length(amino_acids)) + 1L
      ]
      sequence[[8L]] <- site_map$modAA[[index]]
      if (usage_a[[index]] > 0.5) {
        sequence[[9L]] <- "P"
      } else if (usage_a[[index]] < -0.5) {
        sequence[c(5L, 6L)] <- "R"
      }
      if (usage_c[[index]] > 0.5) {
        sequence[c(3L, 4L)] <- c("D", "E")
      } else if (usage_c[[index]] < -0.5) {
        sequence[c(11L, 12L)] <- c("L", "V")
      }
      paste0(sequence, collapse = "")
    },
    character(1)
  )
  site_map$description <- paste("protein", site_map$protein_Id)
  site_map$gene_name <- sub("^P", "GENE", site_map$protein_Id)
  site_map$protein_length <- 300L + protein_index * 7L

  group_effect <- function(group, a, c) {
    ifelse(group == "a", a, ifelse(group == "c", c, 0))
  }
  replicate_noise <- c(-0.18, 0.04, 0.14)

  long_protein <- expand.grid(
    Name = samples$Name,
    protein_Id = proteins,
    stringsAsFactors = FALSE
  )
  protein_sample <- match(long_protein$Name, samples$Name)
  protein_feature <- match(long_protein$protein_Id, proteins)
  long_protein$raw_file <- samples$raw_file[protein_sample]
  long_protein$G_ <- samples$G_[protein_sample]
  long_protein$normalized_abundance <-
    .example_base_abundance(long_protein$protein_Id, proteins) +
    group_effect(
      long_protein$G_,
      protein_a[protein_feature],
      protein_c[protein_feature]
    ) +
    replicate_noise[samples$replicate[protein_sample]] *
      (1 + (protein_feature %% 5L) / 20)

  long_site <- expand.grid(
    Name = samples$Name,
    site = site_map$site,
    stringsAsFactors = FALSE
  )
  site_sample <- match(long_site$Name, samples$Name)
  site_feature <- match(long_site$site, site_map$site)
  long_site$protein_Id <- site_map$protein_Id[site_feature]
  long_site$raw_file <- samples$raw_file[site_sample]
  long_site$G_ <- samples$G_[site_sample]
  long_site$normalized_abundance <-
    .example_base_abundance(long_site$site, site_map$site) +
    group_effect(
      long_site$G_,
      site_a[site_feature],
      site_c[site_feature]
    ) +
    c(-0.23, 0.06, 0.17)[samples$replicate[site_sample]] *
      (1 + (site_feature %% 7L) / 25)
  faint <- site_feature <= 12L & samples$replicate[site_sample] == 3L
  long_site$normalized_abundance[faint] <- NA_real_

  contrasts <- c("a_vs_b", "c_vs_b")
  protein_grid <- expand.grid(
    protein_Id = proteins,
    contrast = contrasts,
    stringsAsFactors = FALSE
  )
  protein_feature <- match(protein_grid$protein_Id, proteins)
  protein_grid$diff <- ifelse(
    protein_grid$contrast == "a_vs_b",
    protein_a[protein_feature],
    protein_c[protein_feature]
  )
  protein_grid$std.error <- 0.18
  protein_grid$df <- 6
  protein_grid$std.error.unmoderated <- 0.22
  protein_grid$df.unmoderated <- 6
  protein_grid$statistic <- protein_grid$diff / protein_grid$std.error
  protein_grid$p_value <- 2 * stats::pt(-abs(protein_grid$statistic), protein_grid$df)
  protein_grid$FDR <- stats::ave(
    protein_grid$p_value,
    protein_grid$contrast,
    FUN = stats::p.adjust,
    method = "BH"
  )
  protein_grid$p_value <- NULL
  protein_grid$description <- paste("protein", protein_grid$protein_Id)
  protein_grid$gene_name <- sub("^P", "GENE", protein_grid$protein_Id)
  protein_grid$protein_length <- 300L + protein_feature * 7L
  protein_grid$estimate_type <- "observed"

  site_grid <- expand.grid(
    site = site_map$site,
    contrast = contrasts,
    stringsAsFactors = FALSE
  )
  site_feature <- match(site_grid$site, site_map$site)
  site_grid <- cbind(
    site_grid,
    site_map[site_feature, setdiff(names(site_map), "site"), drop = FALSE]
  )
  site_grid$diff <- ifelse(
    site_grid$contrast == "a_vs_b",
    site_a[site_feature],
    site_c[site_feature]
  )
  site_grid$std.error <- 0.2
  site_grid$df <- 6
  site_grid$std.error.unmoderated <- 0.24
  site_grid$df.unmoderated <- 6
  site_grid$statistic <- site_grid$diff / site_grid$std.error
  site_grid$p_value <- 2 * stats::pt(-abs(site_grid$statistic), site_grid$df)
  site_grid$FDR <- stats::ave(
    site_grid$p_value,
    site_grid$contrast,
    FUN = stats::p.adjust,
    method = "BH"
  )
  site_grid$p_value <- NULL
  site_grid$estimate_type <- "observed"

  paths <- list(
    phospho = file.path(root, "DEA_phospho"),
    protein = file.path(root, "DEA_protein"),
    annot_file = file.path(root, "annotation.tsv")
  )
  .example_dea_dir(
    paths$phospho,
    long_site,
    .example_dea_config(list(protein_Id = "protein_Id", site = "site")),
    site_grid,
    site_map
  )
  .example_dea_dir(
    paths$protein,
    long_protein,
    .example_dea_config(list(protein_Id = "protein_Id")),
    protein_grid
  )
  readr::write_tsv(
    samples[, c("Name", "G_", "Control")] |>
      dplyr::rename(Group = "G_"),
    paths$annot_file
  )
  paths
}

#' Build the Final MuData Used by Package Vignettes
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
  dirs <- .example_ptm_results_dea_pair(file.path(output_dir, "dea"))
  enriched_h5ad <- file.path(output_dir, "enriched.h5ad")
  total_h5ad <- file.path(output_dir, "total.h5ad")
  input_h5mu <- file.path(output_dir, "inputs.h5mu")
  statistics_h5mu <- file.path(output_dir, "statistics.h5mu")

  .write_example_ptm_dea_h5ad(
    dirs$phospho,
    enriched_h5ad,
    site = TRUE,
    annot_file = dirs$annot_file
  )
  .write_example_ptm_dea_h5ad(
    dirs$protein,
    total_h5ad,
    site = FALSE,
    annot_file = dirs$annot_file
  )
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
