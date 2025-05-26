#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_phospho_integration <- function(workdir = getwd()) {
  runscripts <- c(
    "application/_Overview_PhosphoAndIntegration_site.Rmd",
    "application/bibliography2025.bib"
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prophosqua")
}
