#' Filter significant PTM sites
#'
#' Filters phosphosite data by FDR and fold-change thresholds and adds a
#' `regulation` column, `"upregulated"` or `"downregulated"` by the sign of the
#' fold change, as [plot_diff_logo()] and [plot_seqlogo_with_diff()] expect.
#'
#' @param data Data frame with PTM results containing FDR and fold-change columns
#' @param fdr_col Name of the FDR column.
#' @param diff_col Name of the log2 fold-change column.
#' @param fdr_threshold Sites with FDR below it are kept.
#' @param fc_threshold Sites with absolute log2 fold change above it are kept.
#' @param require_sequence If TRUE, drop rows whose SequenceWindow is NA or
#'   starts or ends with an underscore.
#' @return Filtered data frame with a `regulation` column.
#' @export
#' @examples
#' data <- data.frame(
#'   contrast = rep("A_vs_B", 6),
#'   SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCSCCCC",
#'                      "DDDSDDDD", "EEESEEEE", "FFFSFFF"),
#'   FDR.site = c(0.01, 0.03, 0.08, 0.02, 0.15, 0.04),
#'   diff.site = c(1.2, -0.8, 0.5, -1.5, 0.3, 0.9)
#' )
#' filter_significant_sites(data)$regulation
filter_significant_sites <- function(
  data,
  fdr_col = "FDR.site",
  diff_col = "diff.site",
  fdr_threshold = 0.05,
  fc_threshold = 0.6,
  require_sequence = FALSE
) {
  .require_columns(data, c(fdr_col, diff_col), "PTM results")
  result <- data |>
    dplyr::filter(
      .data[[fdr_col]] < fdr_threshold,
      abs(.data[[diff_col]]) > fc_threshold
    ) |>
    dplyr::mutate(
      regulation = dplyr::case_when(
        .data[[diff_col]] > 0 ~ "upregulated",
        .data[[diff_col]] < 0 ~ "downregulated",
        TRUE ~ "no_change"
      )
    )
  if (require_sequence) {
    .require_columns(result, "SequenceWindow", "PTM results")
    result <- dplyr::filter(
      result,
      !is.na(.data$SequenceWindow),
      !grepl("^_", .data$SequenceWindow),
      !grepl("_$", .data$SequenceWindow)
    )
  }
  result
}

#' Validate sequence window alignment
#'
#' Keeps only rows where the central residue of the sequence window matches the
#' reported modified amino acid, compared case-insensitively.
#'
#' @param data Data frame with PTM results containing SequenceWindow and modAA columns
#' @param seq_col Name of the sequence window column.
#' @param mod_col Name of the modified amino acid column.
#' @param center_pos Position of the central residue (1-indexed).
#' @return Filtered data frame with only valid sequence windows
#' @export
#' @examples
#' data <- data.frame(
#'   SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCACCC"),
#'   modAA = c("S", "S", "S")
#' )
#' validate_sequence_window(data)
validate_sequence_window <- function(data, seq_col = "SequenceWindow", mod_col = "modAA", center_pos = 8L) {
  .require_columns(data, c(seq_col, mod_col), "PTM results")
  matches <- toupper(substr(data[[seq_col]], center_pos, center_pos)) == toupper(data[[mod_col]])
  data[!is.na(matches) & matches, , drop = FALSE]
}
