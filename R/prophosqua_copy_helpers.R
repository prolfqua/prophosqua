#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_phospho_Integration <- function(workdir = getwd()) {
  runscripts <- c("application/_Overview_PhosphoAndIntegration.Rmd",
                  "inst/application/bibliography2025.blib"
                  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prophosqua")
}

