#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_phosphoDEA_FragPipe_TMT <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c("application/_DiffExpQC_Phospho.Rmd",
                  "application/_DiffExpQC_Phospho_V2.Rmd",
                  "application/_Grp2Analysis_Phospho.Rmd",
                  "application/_Grp2Analysis_Phospho_V2.Rmd",
                  "application/_Overview_PhosphoAndIntegration.Rmd",

    if (run_script) {c("application/FP_TMT_phos/FP_TMT_cleanVersion_DEA_enriched_and_total.R",
                       "application/FP_TMT_phos/FP_TMT_cleanVersion_integration.R")}
    )

  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prophosqua")
}


