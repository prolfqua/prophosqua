# MS-stats like normalization for protein change
# https://bioc.ism.ac.jp/packages/3.12/bioc/vignettes/MSstatsTMTPTM/inst/doc/MSstatsTMTPTM.Workflow.html
# look for: how to adjust PTMs

#' Apply MSstatsPTM like site normalization (adjustment) for the protein fold-change and pvalues
#' @param mycombo dataframe from PTMsite and protein prolfqua statistics
#' @export
#'
doMSstatsLikeSiteNormalizationUsingProteinStatsOnComboObject <- function (mycombo)
{
  resultCombo <- data.frame(stringsAsFactors = TRUE)
  for (i in 1:length(unique(mycombo$contrast))) {
    OneC <- mycombo[mycombo$contrast == unique(mycombo$contrast)[i],]
    OneC$MSstatsPTMadj_log2fc <- OneC$diff.x - OneC$diff.y
    OneC$MSstatsPTMadj_s2 <- OneC$std.error.x^2
    OneC$MSstatsPTMadj_s2prot <- OneC$std.error.y^2
    OneC$MSstatsPTMadj_stderr <- sqrt(OneC$MSstatsPTMadj_s2 + OneC$MSstatsPTMadj_s2prot)
    OneC$MSstatsPTMadj_numer <- (OneC$MSstatsPTMadj_s2 + OneC$MSstatsPTMadj_s2prot)^2
    OneC$MSstatsPTMadj_denom <- (OneC$MSstatsPTMadj_s2^2 / OneC$df.x + OneC$MSstatsPTMadj_s2prot^2 / OneC$df.y)
    OneC$MSstatsPTMadj_df <- OneC$MSstatsPTMadj_numer / OneC$MSstatsPTMadj_denom
    OneC$MSstatsPTMadj_tval <- OneC$MSstatsPTMadj_log2fc / OneC$MSstatsPTMadj_stderr
    OneC$MSstatsPTMadj_pVals <- 2 * stats::pt(abs(OneC$MSstatsPTMadj_tval), OneC$MSstatsPTMadj_df, lower.tail = FALSE)
    #adjust pV for multiple testing
    OneC$MSstatsPTMadj_FDR <- p.adjust(OneC$MSstatsPTMadj_pVals, method = "BH")
    resultCombo <- rbind(resultCombo, OneC)
  }
  return(resultCombo)
}

#' Test if differences of differences are significant
#'
#' @export
#'
test_diff_diff <- function(dfA, dfB,
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


