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
#' @param require_positive_df Whether non-finite or non-positive degrees of
#'   freedom make a row untestable.
#' @return Data frame with diff_diff test results
#' @keywords internal
.test_diff_diff <- function(
  dataframe_a,
  dataframe_b,
  by,
  diff = c("diff"),
  std_err = c("std.error"),
  df = c("df"),
  suffix_a = ".site",
  suffix_b = ".protein",
  require_positive_df = FALSE
) {
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

  valid_df <- rep(TRUE, nrow(dataf))
  if (require_positive_df) {
    df_a_values <- dataf[[paste0(df, suffix_a)]]
    df_b_values <- dataf[[paste0(df, suffix_b)]]
    valid_df <- is.finite(df_a_values) & df_a_values > 0 & is.finite(df_b_values) & df_b_values > 0
  }

  dataf <- dataf |>
    dplyr::mutate(
      diff_diff = !!diff_a - !!diff_b,
      SE_I = f_se(!!std_error_a, !!std_error_b),
      df_I = f_df(!!std_error_a, !!std_error_b, !!df_a, !!df_b)
    )
  if (require_positive_df) {
    dataf$df_I[!valid_df] <- NA_real_
  }

  dataf <- dataf |> dplyr::mutate(tstatistic_I = .data$diff_diff / .data$SE_I)
  dataf <- dataf |>
    dplyr::mutate(
      pValue_I = 2 * pt(q = abs(.data$tstatistic_I), df = .data$df_I, lower.tail = FALSE)
    )

  if (require_positive_df) {
    invalid_df <- !valid_df
    dataf$df_I[invalid_df] <- NA_real_
    dataf$tstatistic_I[invalid_df] <- NA_real_
    dataf$pValue_I[invalid_df] <- NA_real_
    dataf <- dataf |>
      dplyr::group_by(.data$contrast) |>
      dplyr::mutate(
        FDR_I = {
          adjusted <- rep(NA_real_, dplyr::n())
          testable <- which(!is.na(.data$pValue_I))
          adjusted[testable] <- p.adjust(.data$pValue_I[testable], method = "BH")
          adjusted
        }
      ) |>
      dplyr::ungroup()
  } else {
    dataf <- dataf |>
      dplyr::group_by(.data$contrast) |>
      dplyr::mutate(FDR_I = p.adjust(.data$pValue_I, method = "BH")) |>
      dplyr::ungroup()
  }
  return(dataf)
}


.reverse_join_column <- function(join_column) {
  reverse_join_column <- vector(mode = "character", length(join_column))
  join_names <- names(join_column)
  if (is.null(join_names)) {
    join_names <- rep("", length(join_column))
  }
  for (i in seq_along(join_column)) {
    reverse_join_column[i] <- if (join_names[i] != "") {
      join_names[i]
    } else {
      join_column[i]
    }
    names(reverse_join_column)[i] <- if (join_names[i] != "") {
      join_column[i]
    } else {
      ""
    }
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
#' @param variant Variance and degrees-of-freedom pair to use. `moderated`
#'   uses the reported moderated pair; `unmoderated` uses the pre-moderation
#'   pair retained in the DEA output.
#' @return Data frame with combined results including diff_diff test statistics
#' @export
test_diff <- function(
  phos_res,
  tot_res,
  join_column = c(
    "protein_Id",
    "contrast",
    "description",
    "protein_length",
    "nr_tryptic_peptides"
  ),
  variant = c("moderated", "unmoderated")
) {
  variant <- match.arg(variant)
  required_schema <- c("std.error", "df", "std.error.unmoderated", "df.unmoderated")
  missing_site <- setdiff(required_schema, colnames(phos_res))
  missing_protein <- setdiff(required_schema, colnames(tot_res))
  if (length(missing_site) > 0 || length(missing_protein) > 0) {
    missing_description <- c(
      if (length(missing_site) > 0) {
        paste0("site result: ", paste(missing_site, collapse = ", "))
      },
      if (length(missing_protein) > 0) {
        paste0("protein result: ", paste(missing_protein, collapse = ", "))
      }
    )
    stop(
      "DPU requires current t/Wald contrast columns; missing ",
      paste(missing_description, collapse = "; "),
      ". Rerun both DEAs with the current prolfqua.",
      call. = FALSE
    )
  }

  std_err <- if (variant == "moderated") "std.error" else "std.error.unmoderated"
  df <- if (variant == "moderated") "df" else "df.unmoderated"
  test_diff <- .test_diff_diff(
    phos_res,
    tot_res,
    by = join_column,
    std_err = std_err,
    df = df,
    require_positive_df = variant == "unmoderated"
  )
  test_diff$measured_In <- "both"

  removed_from_site <- dplyr::anti_join(phos_res, tot_res, by = join_column)
  removed_from_site$measured_In <- rep("site", nrow(removed_from_site))

  removed_from_prot <- dplyr::anti_join(tot_res, phos_res, by = .reverse_join_column(join_column))
  removed_from_prot$measured_In <- rep("prot", nrow(removed_from_prot))

  common_columns <- setdiff(
    intersect(
      colnames(removed_from_site),
      colnames(removed_from_prot)
    ),
    c(join_column, "measured_In")
  )
  removed_from_site_renamed <- removed_from_site |>
    dplyr::rename_with(~ paste0(., ".site"), tidyselect::all_of(common_columns))
  removed_from_prot_renamed <- removed_from_prot |>
    dplyr::rename_with(~ paste0(., ".protein"), dplyr::all_of(common_columns))

  combined_test_diff <- dplyr::bind_rows(
    test_diff,
    removed_from_site_renamed,
    removed_from_prot_renamed
  )
  return(combined_test_diff)
}
