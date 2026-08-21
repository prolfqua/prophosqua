#' copy Markdown and runscript for FragPipe combined_protein.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_phospho_integration <- function(workdir = getwd()) {
  runscripts <- c(
    "application/_Overview_PhosphoAndIntegration_site.Rmd",
    "application/bibliography2025.bib"
  )
  prolfqua::script_copy_helper_vec(runscripts, workdir = workdir, packagename = "prophosqua")
}

#' Copy the PTM Command Wrapper into a Working Directory
#'
#' Places `ptm.sh` from `inst/application/bin` in `workdir`. It is one script
#' taking the step to run as its first argument -- `ptm.sh dpa_dpu`,
#' `ptm.sh render`, `ptm.sh help` -- so a work directory holds a single wrapper
#' however many steps the package grows.
#'
#' The copy carries no analysis logic: it resolves both the list of commands and
#' the script each one runs from the installed package, and so cannot drift from
#' it. It exists so that a person at a prompt runs the pipeline's steps through
#' the same entry point the Snakefile uses.
#'
#' @param workdir directory where to copy the wrapper - default is the current
#'   working directory.
#' @return The path copied, from [prolfqua::script_copy_helper_vec()].
#' @export
#' @examples
#' basename(copy_ptm_shell_script(tempdir()))
copy_ptm_shell_script <- function(workdir = getwd()) {
  copied <- prolfqua::script_copy_helper_vec(
    file.path("application", "bin", "ptm.sh"),
    workdir = workdir,
    packagename = "prophosqua"
  )
  # file.copy preserves the mode, but a package built on a filesystem that does
  # not carry it would leave the wrapper unexecutable.
  Sys.chmod(copied, mode = "0755")
  copied
}
