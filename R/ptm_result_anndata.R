.prophosqua_ptm_schema_version <- "1.0.0"

#' Compute and Write Site-Level PTM Results as AnnData
#'
#' Validates a pair of prolfquapp DEA H5AD files, computes DPA, moderated and
#' unmoderated DPU, and CorrectFirst results, then writes a new site-level H5AD
#' file. The input files are never modified. Site measurements and upstream DEA
#' results are preserved, while PTM statistics are added as feature-aligned
#' `varm` matrices under a versioned `uns["prophosqua"]` contract.
#'
#' @param site_h5ad Path to the site-level prolfquapp DEA H5AD file.
#' @param protein_h5ad Path to the total-proteome prolfquapp DEA H5AD file.
#' @param annot_file Sample annotation defining groups and contrasts.
#' @param output_h5ad Destination for the new site-level PTM result H5AD.
#' @return The normalized output path, invisibly.
#' @export
#' @examples
#' \dontrun{
#' compute_ptm_results_h5ad(
#'   site_h5ad = "phospho/AnnData.h5ad",
#'   protein_h5ad = "total_proteome/AnnData.h5ad",
#'   annot_file = "phospho_dataset.tsv",
#'   output_h5ad = "PTM_results.h5ad"
#' )
#' }
compute_ptm_results_h5ad <- function(
  site_h5ad,
  protein_h5ad,
  annot_file,
  output_h5ad
) {
  .validate_ptm_output_path(output_h5ad, site_h5ad, protein_h5ad)
  if (!file.exists(annot_file)) {
    stop("Annotation file not found: ", annot_file, call. = FALSE)
  }

  pair <- read_ptm_anndata_pair(site_h5ad, protein_h5ad)
  source_hashes <- .ptm_input_hashes(pair)
  annot <- readr::read_tsv(annot_file, show_col_types = FALSE)
  dpa_dpu <- .compute_dpa_dpu_from_pair(pair)
  correct_first <- .compute_cf_dea_from_pair(
    pair,
    annot,
    basename(annot_file)
  )
  .write_ptm_result_h5ad(
    pair,
    dpa_dpu,
    correct_first,
    source_hashes,
    output_h5ad
  )
}

.validate_ptm_output_path <- function(output_h5ad, site_h5ad, protein_h5ad) {
  if (!is.character(output_h5ad) || length(output_h5ad) != 1L || !nzchar(output_h5ad)) {
    stop("output_h5ad must be one non-empty path.", call. = FALSE)
  }
  output_dir <- dirname(output_h5ad)
  if (!dir.exists(output_dir)) {
    stop("PTM result output directory does not exist: ", output_dir, call. = FALSE)
  }

  resolved_output <- file.path(
    normalizePath(output_dir, mustWork = TRUE),
    basename(output_h5ad)
  )
  input_paths <- normalizePath(c(site_h5ad, protein_h5ad), mustWork = FALSE)
  if (resolved_output %in% input_paths) {
    stop("PTM result output must not overwrite either input H5AD file.", call. = FALSE)
  }
}

.write_ptm_result_h5ad <- function(
  pair,
  dpa_dpu,
  correct_first,
  source_hashes,
  output_h5ad
) {
  payload <- .ptm_result_payload(pair, dpa_dpu, correct_first, source_hashes)
  adata <- .ptm_result_anndata(pair$site$source_path, payload)
  temporary <- tempfile(
    pattern = ".PTM-results-",
    tmpdir = dirname(output_h5ad),
    fileext = ".h5ad"
  )
  on.exit(unlink(temporary), add = TRUE)

  invisible(rhdf5::H5get_libversion())
  adata$write_h5ad(temporary, compression = "gzip", mode = "w")
  restored <- anndataR::read_h5ad(temporary)
  .validate_written_ptm_result(restored, adata, payload)
  if (!identical(.ptm_input_hashes(pair), source_hashes)) {
    stop("An input H5AD file changed while the PTM result was being written.", call. = FALSE)
  }
  if (!file.rename(temporary, output_h5ad)) {
    stop("Could not publish validated PTM result H5AD: ", output_h5ad, call. = FALSE)
  }
  invisible(normalizePath(output_h5ad, mustWork = TRUE))
}

.ptm_input_hashes <- function(pair) {
  c(
    site = unname(tools::md5sum(pair$site$source_path)),
    protein = unname(tools::md5sum(pair$protein$source_path))
  )
}

.ptm_result_payload <- function(pair, dpa_dpu, correct_first, source_hashes) {
  site_var <- pair$site$var
  analyses <- list(
    dpa = .ptm_analysis_payload(
      dpa_dpu$combined_site_prot,
      site_var,
      "dpa",
      unique(dpa_dpu$combined_site_prot$contrast)
    ),
    dpu = .ptm_analysis_payload(
      dpa_dpu$combined_test_diff,
      site_var,
      "dpu",
      unique(dpa_dpu$combined_test_diff$contrast)
    ),
    dpu_unmoderated = .ptm_analysis_payload(
      dpa_dpu$combined_test_diff_unmoderated,
      site_var,
      "dpu_unmoderated",
      unique(dpa_dpu$combined_test_diff_unmoderated$contrast)
    ),
    correct_first = .ptm_analysis_payload(
      correct_first$results,
      site_var,
      "correct_first",
      names(correct_first$contrasts)
    )
  )
  combined <- .combine_ptm_analysis_payloads(analyses)
  combined$namespace <- .ptm_result_namespace(
    pair,
    analyses,
    combined,
    source_hashes,
    dpa_dpu,
    correct_first
  )
  combined
}

.ptm_analysis_payload <- function(data, site_var, prefix, contrasts) {
  .require_columns(data, "contrast", paste(prefix, "results"))
  .require_columns(site_var, "site", "site AnnData var")
  contrasts <- unique(as.character(contrasts))
  contrasts <- contrasts[!is.na(contrasts) & nzchar(contrasts)]
  if (length(contrasts) == 0L) {
    stop("No contrasts are available for ", prefix, " results.", call. = FALSE)
  }

  result_site_column <- site_column(data)
  numeric_columns <- names(data)[vapply(
    data,
    function(column) is.numeric(column) || is.logical(column),
    logical(1)
  )]
  numeric_columns <- setdiff(numeric_columns, names(site_var))
  if (length(numeric_columns) == 0L) {
    stop("No numeric statistics are available for ", prefix, " results.", call. = FALSE)
  }
  annotation_columns <- names(data)[
    !vapply(
      data,
      function(column) is.numeric(column) || is.logical(column),
      logical(1)
    )
  ]
  annotation_columns <- setdiff(
    annotation_columns,
    c(names(site_var), "contrast")
  )

  payload <- list(values = list(), columns = list(), annotations = list(), present = list())
  payload$keys <- stats::setNames(character(length(contrasts)), contrasts)
  for (contrast in contrasts) {
    key <- .ptm_varm_key(prefix, contrast)
    aligned <- .align_ptm_result_table(
      data[data$contrast == contrast, , drop = FALSE],
      site_var,
      result_site_column,
      numeric_columns,
      annotation_columns,
      key
    )
    payload$values[[key]] <- aligned$values
    payload$columns[[key]] <- numeric_columns
    payload$annotations[[key]] <- aligned$annotations
    payload$present[[key]] <- aligned$present
    payload$keys[[contrast]] <- key
  }
  payload
}

.ptm_varm_key <- function(prefix, contrast) {
  paste0(
    prefix,
    "__",
    utils::URLencode(contrast, reserved = TRUE, repeated = TRUE)
  )
}

.align_ptm_result_table <- function(
  data,
  site_var,
  result_site_column,
  numeric_columns,
  annotation_columns,
  key
) {
  site_axis <- as.character(site_var$site)
  result_sites <- as.character(data[[result_site_column]])
  aligned_rows <- !is.na(result_sites) & nzchar(result_sites)
  data <- data[aligned_rows, , drop = FALSE]
  result_sites <- result_sites[aligned_rows]

  unknown_sites <- setdiff(result_sites, site_axis)
  if (length(unknown_sites) > 0L) {
    stop(
      "PTM result '",
      key,
      "' contains site(s) absent from the site AnnData axis: ",
      paste(unknown_sites, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(result_sites)) {
    stop("PTM result '", key, "' has duplicate rows for one site.", call. = FALSE)
  }

  positions <- match(result_sites, site_axis)
  values <- matrix(
    NA_real_,
    nrow = nrow(site_var),
    ncol = length(numeric_columns),
    dimnames = list(rownames(site_var), numeric_columns)
  )
  if (length(positions) > 0L) {
    values[positions, ] <- as.matrix(data[, numeric_columns, drop = FALSE])
  }

  annotations <- lapply(
    annotation_columns,
    function(column) {
      .align_ptm_annotation(data[[column]], positions, nrow(site_var))
    }
  )
  names(annotations) <- annotation_columns
  present <- rep(FALSE, nrow(site_var))
  present[positions] <- TRUE
  list(values = values, annotations = annotations, present = present)
}

.align_ptm_annotation <- function(values, positions, size) {
  if (is.factor(values)) {
    values <- as.character(values)
  }
  if (is.logical(values)) {
    aligned <- rep(NA, size)
  } else if (is.integer(values)) {
    aligned <- rep(NA_integer_, size)
  } else if (is.numeric(values)) {
    aligned <- rep(NA_real_, size)
  } else {
    values <- as.character(values)
    aligned <- rep(NA_character_, size)
  }
  aligned[positions] <- values
  unname(aligned)
}

.combine_ptm_analysis_payloads <- function(analyses) {
  combined <- list(values = list(), columns = list(), annotations = list(), present = list())
  for (analysis in analyses) {
    duplicate_keys <- intersect(names(combined$values), names(analysis$values))
    if (length(duplicate_keys) > 0L) {
      stop(
        "PTM result keys collide: ",
        paste(duplicate_keys, collapse = ", "),
        call. = FALSE
      )
    }
    combined$values <- c(combined$values, analysis$values)
    combined$columns <- c(combined$columns, analysis$columns)
    combined$annotations <- c(combined$annotations, analysis$annotations)
    combined$present <- c(combined$present, analysis$present)
  }
  combined
}

.ptm_result_namespace <- function(
  pair,
  analyses,
  combined,
  source_hashes,
  dpa_dpu,
  correct_first
) {
  list(
    artifact_type = "ptm_results",
    schema_version = .prophosqua_ptm_schema_version,
    source_software = "prophosqua",
    site_input = .ptm_input_provenance(pair$site, source_hashes[["site"]]),
    protein_input = .ptm_input_provenance(
      pair$protein,
      source_hashes[["protein"]]
    ),
    result_keys = lapply(analyses, function(analysis) analysis$keys),
    varm_columns = combined$columns,
    varm_annotations = combined$annotations,
    varm_present = combined$present,
    methods = list(
      dpa = list(name = "Differential PTM Abundance"),
      dpu = list(name = "Differential PTM Usage", variance = "moderated"),
      dpu_unmoderated = list(
        name = "Differential PTM Usage",
        variance = "unmoderated"
      ),
      correct_first = list(
        name = "CorrectFirst",
        formula = "ptm_usage ~ G_",
        contrasts = unname(as.character(correct_first$contrasts))
      )
    ),
    result_counts = list(
      dpa = nrow(dpa_dpu$combined_site_prot),
      dpu = nrow(dpa_dpu$combined_test_diff),
      dpu_unmoderated = nrow(dpa_dpu$combined_test_diff_unmoderated),
      correct_first = nrow(correct_first$results)
    )
  )
}

.ptm_input_provenance <- function(experiment, hash) {
  list(
    path = experiment$source_path,
    hash_algorithm = "md5",
    hash = unname(hash),
    schema_version = experiment$schema_version
  )
}

.ptm_result_anndata <- function(site_h5ad, payload) {
  adata <- anndataR::read_h5ad(site_h5ad)
  if (!is.null(adata$uns[["prophosqua"]])) {
    stop("Site AnnData already contains a prophosqua result namespace.", call. = FALSE)
  }
  duplicate_keys <- intersect(adata$varm_keys(), names(payload$values))
  if (length(duplicate_keys) > 0L) {
    stop(
      "Site AnnData already contains PTM result key(s): ",
      paste(duplicate_keys, collapse = ", "),
      call. = FALSE
    )
  }

  for (key in names(payload$values)) {
    adata$varm[[key]] <- payload$values[[key]]
  }
  adata$uns[["prophosqua"]] <- payload$namespace
  adata
}

.validate_written_ptm_result <- function(restored, expected, payload) {
  .validate_ptm_axes(restored, expected)
  .validate_ptm_layers(restored, expected)
  .validate_ptm_varm(restored, expected)
  .validate_ptm_namespace(restored, payload)
}

.validate_ptm_axes <- function(restored, expected) {
  if (!identical(rownames(as.data.frame(restored$obs)), rownames(as.data.frame(expected$obs)))) {
    stop("Written PTM result changed the sample axis.", call. = FALSE)
  }
  if (!identical(rownames(as.data.frame(restored$var)), rownames(as.data.frame(expected$var)))) {
    stop("Written PTM result changed the feature axis.", call. = FALSE)
  }
  .require_equal_ptm_value(as.matrix(restored$X), as.matrix(expected$X), "X")
}

.validate_ptm_layers <- function(restored, expected) {
  if (!setequal(restored$layers_keys(), expected$layers_keys())) {
    stop("Written PTM result changed the input layer keys.", call. = FALSE)
  }
  for (key in expected$layers_keys()) {
    .require_equal_ptm_value(
      as.matrix(restored$layers[[key]]),
      as.matrix(expected$layers[[key]]),
      paste0("layer '", key, "'")
    )
  }
}

.validate_ptm_varm <- function(restored, expected) {
  if (!setequal(restored$varm_keys(), expected$varm_keys())) {
    stop("Written PTM result changed the expected result keys.", call. = FALSE)
  }
  for (key in expected$varm_keys()) {
    .require_equal_ptm_value(
      as.matrix(restored$varm[[key]]),
      as.matrix(expected$varm[[key]]),
      paste0("varm '", key, "'")
    )
  }
}

.validate_ptm_namespace <- function(restored, payload) {
  namespace <- restored$uns[["prophosqua"]]
  if (
    !is.list(namespace) ||
      !identical(namespace$artifact_type, "ptm_results") ||
      !identical(namespace$schema_version, .prophosqua_ptm_schema_version)
  ) {
    stop("Written PTM result has an invalid prophosqua namespace.", call. = FALSE)
  }
  if (!setequal(names(payload$values), names(namespace$varm_columns))) {
    stop("Written PTM result metadata does not cover every PTM result matrix.", call. = FALSE)
  }
}

.require_equal_ptm_value <- function(actual, expected, label) {
  comparison <- all.equal(actual, expected, check.attributes = FALSE)
  if (!isTRUE(comparison)) {
    stop("Written PTM result changed ", label, ": ", comparison[[1]], call. = FALSE)
  }
}
