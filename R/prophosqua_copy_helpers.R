#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_DEA_FragPipe_DDA <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c("application/_DiffExpQC_Phospho.Rmd",
                  "application/_Grp2Analysis_Phospho.Rmd",
                  "application/_Overview_PhosphoAndIntegration.Rmd"
    if (run_script) {"application/FP_TMT_phos/FP_DDA.R"}
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}


