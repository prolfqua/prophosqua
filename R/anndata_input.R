# prophosqua reads a DEA artifact only through prolfquapp's DEAResultReader.
.read_dea_pair <- function(phospho_dea_dir, protein_dea_dir, remove_contaminants = FALSE) {
  reader <- function(dea_dir) prolfquapp::DEAResultReader$new(get_dea_file(dea_dir, "AnnData.h5ad"))
  .ptm_pair(reader(phospho_dea_dir), reader(protein_dea_dir), remove_contaminants)
}

.ptm_pair <- function(site_reader, protein_reader, remove_contaminants = FALSE) {
  site <- .without_decoys(.dea_record(site_reader))
  protein <- .without_decoys(.dea_record(protein_reader))
  if (remove_contaminants) {
    site <- .without_contaminants(site)
    protein <- .without_contaminants(protein)
  }
  .validate_site_experiment(site)
  .validate_protein_experiment(protein)
  site$site_info <- .site_info_from_var(site$var)
  .validate_paired_samples(site, protein)
  .validate_shared_design(site, protein)
  list(site = site, protein = protein)
}

# What the PTM analyses take from one DEA artifact, every table keyed by the
# artifact's feature keys. The abundance column is renamed to
# normalized_abundance, the name every PTM computation uses.
.dea_record <- function(reader) {
  metadata <- reader$metadata
  configuration <- reader$lfq_transformed$get_config()$clone(deep = TRUE)
  configuration$set_response("normalized_abundance")
  list(
    sample_key = metadata$sample_key,
    feature_keys = reader$subject_id,
    formula = reader$formula,
    contrasts = reader$contrast_definitions,
    configuration = configuration,
    obs = reader$samples,
    var = reader$annotation,
    normalized_abundances = .dea_abundances(reader$lfq_transformed),
    imputed_abundances = if (!is.null(reader$lfq_imputed)) .dea_abundances(reader$lfq_imputed),
    imputation = reader$imputation,
    differential_results = dplyr::left_join(reader$contrast_table, reader$annotation, by = reader$subject_id)
  )
}

# The DEA keeps contaminants and flags them in its annotation's CON column.
.without_contaminants <- function(record) {
  .require_columns(record$var, "CON", "DEA annotation")
  .without_features(record, record$var$CON)
}

# The DEA fits no decoy but exports their abundances; they are recognized by
# the decoy pattern the DEA recorded, as prolfqua recognizes them.
.without_decoys <- function(record) {
  .without_features(record, prolfqua::is_decoy(record$var$protein_Id, record$configuration$pattern_decoys))
}

# Drop the features of the annotation marked in `drop` from every table.
.without_features <- function(record, drop) {
  dropped <- record$var[drop, record$feature_keys, drop = FALSE]
  record$var <- record$var[!drop, , drop = FALSE]
  for (table in c("normalized_abundances", "imputed_abundances", "imputation", "differential_results")) {
    if (!is.null(record[[table]])) {
      record[[table]] <- dplyr::anti_join(record[[table]], dropped, by = record$feature_keys)
    }
  }
  record
}

.dea_abundances <- function(lfq) {
  dplyr::rename(lfq$data_long(), normalized_abundance = tidyselect::all_of(lfq$response()))
}

.validate_site_experiment <- function(experiment) {
  missing_keys <- setdiff(c("protein_Id", "site"), experiment$feature_keys)
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
    c("protein_Id", "site", "posInProtein", "modAA", "SequenceWindow", "description", "protein_length"),
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
  .require_columns(experiment$var, c("protein_Id", "description", "protein_length"), "protein AnnData var")
}

.require_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

.site_info_from_var <- function(var) {
  columns <- c("site", "posInProtein", "modAA", "SequenceWindow", "protein_Id", "gene_name", "protein_length")
  site_info <- unique(var[, intersect(columns, names(var)), drop = FALSE])
  if (anyDuplicated(site_info$site)) {
    stop("PTM site metadata is not unique by site.", call. = FALSE)
  }
  site_info
}

.validate_paired_samples <- function(site, protein) {
  site_ids <- as.character(site$obs[[site$sample_key]])
  protein_ids <- as.character(protein$obs[[protein$sample_key]])
  if (!setequal(site_ids, protein_ids)) {
    listed <- function(values) if (length(values) == 0L) "none" else paste(values, collapse = ", ")
    stop(
      "Site and protein AnnData sample sets differ; missing from protein: ",
      listed(setdiff(site_ids, protein_ids)),
      "; missing from site: ",
      listed(setdiff(protein_ids, site_ids)),
      ".",
      call. = FALSE
    )
  }
}

.validate_shared_design <- function(site, protein) {
  factors <- names(site$configuration$factors)
  if (!setequal(factors, names(protein$configuration$factors))) {
    stop("Site and protein AnnData files declare different design factors.", call. = FALSE)
  }
  .require_columns(site$obs, factors, "site AnnData obs")
  .require_columns(protein$obs, factors, "protein AnnData obs")
  design <- dplyr::inner_join(
    site$obs[c(site$sample_key, factors)],
    protein$obs[c(protein$sample_key, factors)],
    by = stats::setNames(protein$sample_key, site$sample_key),
    suffix = c(".site", ".protein")
  )
  for (factor_name in factors) {
    if (
      !identical(
        as.character(design[[paste0(factor_name, ".site")]]),
        as.character(design[[paste0(factor_name, ".protein")]])
      )
    ) {
      stop("Site and protein AnnData disagree on design factor '", factor_name, "'.", call. = FALSE)
    }
  }
}
