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
