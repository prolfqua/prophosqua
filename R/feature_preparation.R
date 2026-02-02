# Feature preparation helpers for phosphoproteomics analysis

#' Filter significant PTM sites
#'
#' Filters phosphosite data by FDR and fold-change thresholds and adds a
#' regulation column indicating direction of change. This prepares data for
#' downstream analysis like sequence logo visualization with
#' \code{\link{plot_diff_logo}} and \code{\link{plot_seqlogo_with_diff}}.
#'
#' @param data Data frame with PTM results containing FDR and fold-change columns
#' @param fdr_col Character. Name of FDR/adjusted p-value column. Default "FDR.site"
#' @param diff_col Character. Name of log2 fold-change column. Default "diff.site"
#' @param fdr_threshold Numeric. FDR cutoff (sites with FDR < threshold kept).
#'   Default 0.05
#' @param fc_threshold Numeric. Absolute log2 fold-change cutoff (sites with
#'   |FC| > threshold kept). Default 0.6 (approximately 1.5-fold change)
#' @param require_sequence Logical. If TRUE, filter out rows with invalid
#'   SequenceWindow (NA or starting/ending with underscore). Default FALSE
#'
#' @return Filtered data frame with added 'regulation' column containing
#'   "upregulated" or "downregulated" based on sign of fold-change
#'
#' @details
#' The function performs two operations:
#' \enumerate{
#'   \item Filters to significant sites based on FDR and fold-change thresholds
#'   \item Adds a 'regulation' column based on the sign of the fold-change
#' }
#'
#' The resulting data frame is suitable for use with \code{\link{plot_diff_logo}}
#' and \code{\link{plot_seqlogo_with_diff}} which expect a 'regulation' column.
#'
#' @export
#'
#' @examples
#' # Example with mock data
#' data <- data.frame(
#'   contrast = rep("A_vs_B", 6),
#'   SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCSCCCC",
#'                      "DDDSDDDD", "EEESEEEE", "FFFSFFF"),
#'   FDR.site = c(0.01, 0.03, 0.08, 0.02, 0.15, 0.04),
#'   diff.site = c(1.2, -0.8, 0.5, -1.5, 0.3, 0.9),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Filter with default thresholds (FDR < 0.05, |FC| > 0.6)
#' sig_sites <- filter_significant_sites(data)
#' sig_sites$regulation
#'
#' # Use with DPU column names
#' data_dpu <- data.frame(
#'   contrast = "A_vs_B",
#'   SequenceWindow = c("AAASAAAA", "BBBSBBB"),
#'   FDR_I = c(0.01, 0.03),
#'   diff_diff = c(1.2, -0.8)
#' )
#' sig_dpu <- filter_significant_sites(data_dpu, fdr_col = "FDR_I", diff_col = "diff_diff")
#'
#' @seealso \code{\link{plot_diff_logo}}, \code{\link{plot_seqlogo_with_diff}}
filter_significant_sites <- function(data,
                                     fdr_col = "FDR.site",
                                     diff_col = "diff.site",
                                     fdr_threshold = 0.05,
                                     fc_threshold = 0.6,
                                     require_sequence = FALSE) {
  # Validate columns exist
 required_cols <- c(fdr_col, diff_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

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
    if (!"SequenceWindow" %in% colnames(result)) {
      warning("require_sequence = TRUE but 'SequenceWindow' column not found")
    } else {
      result <- result |>
        dplyr::filter(
          !is.na(.data$SequenceWindow),
          !grepl("^_", .data$SequenceWindow),
          !grepl("_$", .data$SequenceWindow)
        )
    }
  }

  return(result)
}


#' Summarize significant sites by group
#'
#' Generates summary statistics for significant sites, counting up- and
#' down-regulated sites by specified grouping variables.
#'
#' @param data Data frame with 'regulation' column (typically output from
#'   \code{\link{filter_significant_sites}})
#' @param group_cols Character vector of column names to group by.
#'   Default \code{c("contrast")}
#'
#' @return A tibble with counts pivoted wide, showing upregulated and
#'   downregulated counts per group
#'
#' @export
#'
#' @examples
#' # Example with mock data
#' data <- data.frame(
#'   contrast = c(rep("A_vs_B", 4), rep("C_vs_D", 3)),
#'   regulation = c("upregulated", "upregulated", "downregulated", "upregulated",
#'                  "downregulated", "downregulated", "upregulated")
#' )
#'
#' summarize_significant_sites(data)
#'
#' @seealso \code{\link{filter_significant_sites}}
summarize_significant_sites <- function(data, group_cols = c("contrast")) {
  if (!"regulation" %in% colnames(data)) {
    stop("Data must contain 'regulation' column. Use filter_significant_sites() first.")
  }

  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "regulation")))) |>
    dplyr::summarize(n = dplyr::n(), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = "regulation",
      values_from = "n",
      values_fill = 0
    )
}
