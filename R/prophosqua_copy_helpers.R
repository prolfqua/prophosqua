#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_phosphoDEA_FragPipe_TMT <- function(workdir = getwd()) {
  runscripts <- c("application/_DiffExpQC_Phospho_V2.Rmd",
                  "application/_Grp2Analysis_Phospho_V2.Rmd",
                  "application/_Overview_PhosphoAndIntegration.Rmd",
                  "inst/application/bibliography2025.blib"
                  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prophosqua")
}

#' copy Markdown and runscript for Spectronaut Ubi workflow
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_ubiDEA_Spectronaut <- function(workdir = getwd()) {
  runscripts <- c(
    "inst/application/_Overview_Ubi_Integration_WEW.Rmd",
    "inst/application/_Grp2Analysis_Ubi.Rmd",
    "inst/application/_DiffExpQC_Ubi.Rmd",
    "inst/application/bibliography2025.blib"
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prophosqua")
}

#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_phosphoDEA_Spectronaut <- function(workdir = getwd()) {
  runscripts <- c("application/_DiffExpQC_Phospho_V2.Rmd",
                  "application/_Grp2Analysis_Phospho_V2.Rmd",
                  "application/_Overview_PhosphoAndIntegration_site.Rmd",
                  "inst/application/bibliography2025.blib"
  )
  prolfqua::scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prophosqua")
}
