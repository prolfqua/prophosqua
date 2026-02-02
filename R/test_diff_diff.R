#' @importFrom stats pt p.adjust
#' @importFrom rlang .data
NULL

#' Test if differences of differences are significant (internal)
#'
#' @param dataframe_a First data frame (e.g., site-level results)
#' @param dataframe_b Second data frame (e.g., protein-level results)
#' @param by Columns to join by
#' @param diff Column name for difference values
#' @param std_err Column name for standard error values
#' @param df Column name for degrees of freedom
#' @param suffix_a Suffix for columns from dataframe_a
#' @param suffix_b Suffix for columns from dataframe_b
#' @return Data frame with diff_diff test results
#' @keywords internal
.test_diff_diff <- function(dataframe_a, dataframe_b,
                            by,
                            diff = c("diff"),
                            std_err = c("std.error"),
                            df = c("df"),
                            suffix_a = ".site",
                            suffix_b = ".protein") {
  dataf <- dplyr::inner_join(dataframe_a, dataframe_b, by = by, suffix = c(suffix_a, suffix_b))

  f_se <- function(stde_a, stde_b) {
    sqrt(stde_a^2 + stde_b^2)
  }
  f_df <- function(stde_a, stde_b, df_a, df_b) {
    (stde_a^2 + stde_b^2)^2 / ((stde_a^4 / df_a + stde_b^4 / df_b))
  }

  diff_a <- rlang::sym(paste0(diff, suffix_a))
  diff_b <- rlang::sym(paste0(diff, suffix_b))
  std_error_a <- rlang::sym(paste0(std_err, suffix_a))
  std_error_b <- rlang::sym(paste0(std_err, suffix_b))
  df_a <- rlang::sym(paste0(df, suffix_a))
  df_b <- rlang::sym(paste0(df, suffix_b))

  dataf <- dataf |> dplyr::mutate(
    diff_diff = !!diff_a - !!diff_b,
    SE_I = f_se(!!std_error_a, !!std_error_b),
    df_I = f_df(!!std_error_a, !!std_error_b, !!df_a, !!df_b)
  )

  dataf <- dataf |> dplyr::mutate(tstatistic_I = .data$diff_diff / .data$SE_I)
  dataf <- dataf |> dplyr::mutate(pValue_I = 2 * pt(q = abs(.data$tstatistic_I), df = .data$df_I, lower.tail = FALSE))
  dataf <- dataf |>
    dplyr::group_by(.data$contrast) |>
    dplyr::mutate(FDR_I = p.adjust(.data$pValue_I, method = "BH")) |>
    dplyr::ungroup()
  return(dataf)
}


.reverse_join_column <- function(join_column){
  reverse_join_column <- vector(mode = "character", length(join_column))
  for (i in seq_along(join_column)) {
    reverse_join_column[i] <- if (names(join_column)[i] != "") {  names(join_column)[i]} else { join_column[i]}
    names(reverse_join_column)[i] <- if (names(join_column)[i] != "") { join_column[i]} else {""}
  }
  return(reverse_join_column)
}


#' Compute MSstats-like test statistics for differential PTM usage
#'
#' Joins site-level and protein-level results, computes difference-of-differences,
#' and performs t-tests for differential PTM usage.
#'
#' @param phos_res Data frame with phospho site-level results
#' @param tot_res Data frame with protein-level results
#' @param join_column Character vector of columns to join by
#' @return Data frame with combined results including diff_diff test statistics
#' @export
test_diff <- function(phos_res,
                      tot_res,
                      join_column = c(
                        "protein_Id",
                        "contrast",
                        "description",
                        "protein_length",
                        "nr_tryptic_peptides"
                      )) {
  test_diff <- .test_diff_diff(phos_res, tot_res, by = join_column)
  test_diff$measured_In <- "both"

  removed_from_site <- dplyr::anti_join(phos_res, tot_res, by = join_column)
  removed_from_site$measured_In <- "site"

  removed_from_prot <- dplyr::anti_join(tot_res, phos_res, by = .reverse_join_column(join_column))
  removed_from_prot$measured_In <- "prot"

  common_columns <- setdiff(
    intersect(
      colnames(removed_from_site), colnames(removed_from_prot)
    ),
    c(join_column, "measured_In")
  )
  removed_from_site_renamed <- removed_from_site |>
    dplyr::rename_with(~ paste0(., ".site"), tidyselect::all_of(common_columns))
  removed_from_prot_renamed <- removed_from_prot |>
    dplyr::rename_with(~ paste0(., ".protein"), dplyr::all_of(common_columns))

  combined_test_diff <- dplyr::bind_rows(test_diff, removed_from_site_renamed, removed_from_prot_renamed)
  return(combined_test_diff)
}
