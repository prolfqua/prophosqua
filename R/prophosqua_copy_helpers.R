#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_phosphoDEA_FragPipe_TMT <- function(workdir = getwd()) {
  runscripts <- c("application/_DiffExpQC_Phospho_V2.Rmd",
                  "application/_Grp2Analysis_Phospho_V2.Rmd")
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prophosqua")
}


