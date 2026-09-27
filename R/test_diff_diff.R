#' @importFrom stats pt p.adjust
#' @importFrom rlang .data
NULL

# The difference of two differences, its standard error, and Welch-Satterthwaite
# degrees of freedom. With require_positive_df, a row whose degrees of freedom
# are not finite and positive has no test.
.test_diff_diff <- function(
  dataframe_a,
  dataframe_b,
  by,
  diff = "diff",
  std_err = "std.error",
  df = "df",
  suffix_a = ".site",
  suffix_b = ".protein",
  require_positive_df = FALSE
) {
  dataf <- dplyr::inner_join(dataframe_a, dataframe_b, by = by, suffix = c(suffix_a, suffix_b))
  se_a <- dataf[[paste0(std_err, suffix_a)]]
  se_b <- dataf[[paste0(std_err, suffix_b)]]
  df_a <- dataf[[paste0(df, suffix_a)]]
  df_b <- dataf[[paste0(df, suffix_b)]]
  dataf$diff_diff <- dataf[[paste0(diff, suffix_a)]] - dataf[[paste0(diff, suffix_b)]]
  dataf$SE_I <- sqrt(se_a^2 + se_b^2)
  dataf$df_I <- (se_a^2 + se_b^2)^2 / (se_a^4 / df_a + se_b^4 / df_b)
  dataf$tstatistic_I <- dataf$diff_diff / dataf$SE_I
  if (require_positive_df) {
    untestable <- !(is.finite(df_a) & df_a > 0 & is.finite(df_b) & df_b > 0)
    dataf$df_I[untestable] <- NA_real_
    dataf$tstatistic_I[untestable] <- NA_real_
  }
  dataf$pValue_I <- 2 * pt(q = abs(dataf$tstatistic_I), df = dataf$df_I, lower.tail = FALSE)
  # p.adjust leaves NA p-values out of the number of tests.
  dataf |>
    dplyr::group_by(.data$contrast) |>
    dplyr::mutate(FDR_I = p.adjust(.data$pValue_I, method = "BH")) |>
    dplyr::ungroup()
}

#' Compute MSstats-like test statistics for differential PTM usage
#'
#' Joins site-level and protein-level results, computes difference-of-differences,
#' and performs t-tests for differential PTM usage. Only a site whose protein
#' has a result has a usage difference, so the table holds the matched pairs.
#'
#' @param phos_res Data frame with phospho site-level results
#' @param tot_res Data frame with protein-level results
#' @param join_column Character vector of columns to join by
#' @param variant Variance and degrees-of-freedom pair to use. `moderated`
#'   uses the reported moderated pair; `unmoderated` uses the pre-moderation
#'   pair retained in the DEA output.
#' @return Data frame with one row per matched site and protein, including the
#'   diff_diff test statistics
#' @export
test_diff <- function(
  phos_res,
  tot_res,
  join_column = c("protein_Id", "contrast", "description", "protein_length", "nr_tryptic_peptides"),
  variant = c("moderated", "unmoderated")
) {
  variant <- match.arg(variant)
  std_err <- if (variant == "moderated") "std.error" else "std.error.unmoderated"
  df <- if (variant == "moderated") "df" else "df.unmoderated"
  .require_columns(phos_res, c(std_err, df), "DPU site result")
  .require_columns(tot_res, c(std_err, df), "DPU protein result")
  .test_diff_diff(
    phos_res,
    tot_res,
    by = join_column,
    std_err = std_err,
    df = df,
    require_positive_df = variant == "unmoderated"
  )
}
