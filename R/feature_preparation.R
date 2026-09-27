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
#'   padded at a protein terminus.
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
      !.is_padded_window(.data$SequenceWindow)
    )
  }
  result
}

# The PTM readers cut every window from the FASTA, padding it with X beyond a
# protein terminus.
.is_padded_window <- function(windows) grepl("^X|X$", windows)
