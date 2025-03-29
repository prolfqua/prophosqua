#' Test if differences of differences are significant
#'
.test_diff_diff <- function(dfA, dfB,
                           by,
                           diff = c("diff"),
                           std.err = c("std.error"),
                           df = c("df"),
                           suffixA = ".site",
                           suffixB = ".protein"
){
  dataf <- dplyr::inner_join(dfA, dfB, by = by, suffix = c(suffixA,suffixB))

  f_SE <- function(stdeA, stdeB){
    sqrt(stdeA ^ 2 + stdeB ^ 2 )
  }
  f_df <- function(stdeA, stdeB, dfA, dfB){
    (stdeA ^ 2 + stdeB ^ 2 )^2 / ((stdeA^4/dfA + stdeB^4/dfB ))
  }

  diff.A = rlang::sym(paste0(diff,suffixA))
  diff.B = rlang::sym(paste0(diff,suffixB))
  std.error.A = rlang::sym(paste0(std.err, suffixA))
  std.error.B = rlang::sym(paste0(std.err, suffixB))
  df.A = rlang::sym(paste0(df, suffixA))
  df.B = rlang::sym(paste0(df, suffixB))

  dataf <- dataf |> dplyr::mutate(
    diff_diff = !!diff.A - !!diff.B,
    SE_I = f_SE(!!std.error.A, !!std.error.B),
    df_I = f_df(!!std.error.A,!!std.error.B,!!df.A,!!df.B ))

  dataf <- dataf |> dplyr::mutate(tstatistic_I = diff_diff / SE_I)
  dataf <- dataf |> dplyr::mutate(pValue_I = 2*pt(q = abs(tstatistic_I), df = df_I, lower.tail = FALSE))
  dataf <- dataf |> dplyr::group_by(contrast) |> dplyr::mutate(FDR_I = p.adjust(pValue_I, method = "BH")) |> dplyr::ungroup()
  return(dataf)

}



#' compute MSstats like test statsitics
#' @export
#'
test_diff <- function(phosRes,
                      totRes,
                      join_column = c("protein_Id", "contrast","description", "protein_length", "nr_tryptic_peptides")){
  test_diff <- .test_diff_diff(phosRes,totRes, by = join_column)
  test_diff$measured_In <- "both"

  removed_from_site <- dplyr::anti_join(phosRes, totRes, by = join_column )
  removed_from_site$measured_In <- "site"

  removed_from_prot <- dplyr::anti_join(totRes, phosRes, by = .reverse_join_column(join_column))
  removed_from_prot$measured_In <- "prot"

  common_columns <- setdiff(intersect(colnames(removed_from_site), colnames(removed_from_prot)),c(join_column, "measured_In"))
  removed_from_site_renamed <- removed_from_site |>
    dplyr::rename_with(~ paste0(., ".site"), tidyselect::all_of(common_columns))
  removed_from_prot_renamed <- removed_from_prot |>
    dplyr::rename_with(~ paste0(., ".protein"), dplyr::all_of(common_columns))

  combined_test_diff <- dplyr::bind_rows(test_diff , removed_from_site_renamed , removed_from_prot_renamed)
  return(combined_test_diff)
}


