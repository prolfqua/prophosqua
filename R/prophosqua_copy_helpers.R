#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_phosphoDEA_FragPipe_TMT <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c("/inst/application/_DiffExpQC_Phospho.Rmd",
                  "/inst/application/_Grp2Analysis_Phospho.Rmd",
                  "/inst/application/_Overview_PhosphoAndIntegration.Rmd"
    if (run_script) {c("/inst/application/FP_TMT_phos/FP_TMT_cleanVersion_DEA_enriched_and_total.R",
                       "/inst/application/FP_TMT_phos/FP_TMT_cleanVersion_integration.R"}
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


