#' Render a Report Shipped with prophosqua
#'
#' Renders one of the installed Quarto reports, `ptm_statistics.qmd`,
#' `ptm_enrichment.qmd` or the landing page `ptm_index.qmd`, so a project never
#' carries a copy of a template.
#'
#' The template is staged into a private directory under `output_dir`: Quarto
#' names its intermediates after the input document, and two renders of the
#' same template running at once would otherwise overwrite each other's.
#'
#' @param name File name of the report, e.g. `"ptm_statistics.qmd"`.
#' @param output_file File name to write, e.g. `"ptm_statistics.html"`.
#' @param output_dir Directory to write the report to.
#' @param params Named list passed to the report's `params`.
#' @return Invisibly, the path of the rendered file.
#' @export
#' @examples
#' # Renders one of the installed templates; needs the data it asks for.
#' \dontrun{
#' render_ptm_report(
#'   "ptm_statistics.qmd", "ptm_statistics.html", "PTM_results",
#'   params = list(input_h5mu = "PTM_statistics.h5mu", fdr_threshold = 0.25)
#' )
#' }
render_ptm_report <- function(name, output_file, output_dir, params = list()) {
  source <- report_file(name)
  if (!is.null(params$input_h5mu)) {
    params$input_h5mu <- normalizePath(params$input_h5mu, mustWork = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  render_dir <- tempfile(".prophosqua_qmd_", tmpdir = normalizePath(output_dir, mustWork = TRUE))
  dir.create(render_dir)
  on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)
  staged <- file.path(render_dir, name)
  if (!file.copy(source, staged)) {
    stop("Could not stage Quarto report: ", source, call. = FALSE)
  }
  fgczQuartoTemplate::fgcz_render(staged, output_file = output_file, execute_params = params, quiet = FALSE)
  destination <- file.path(output_dir, output_file)
  if (!file.copy(file.path(render_dir, output_file), destination, overwrite = TRUE)) {
    stop("Could not copy Quarto report to: ", destination, call. = FALSE)
  }
  invisible(destination)
}

#' Path of an Installed Report Template
#'
#' The reports are the package's vignettes, and the vignette machinery installs
#' their sources into `doc/`; a package installed without its vignettes built
#' cannot render them.
#'
#' @param name File name, e.g. `"ptm_statistics.qmd"`.
#' @return Full path to the installed template.
#' @export
#' @examples
#' \dontrun{
#' report_file("ptm_statistics.qmd")
#' }
report_file <- function(name) {
  path <- system.file("doc", name, package = "prophosqua")
  if (!nzchar(path)) {
    stop(
      "prophosqua report template not found: ",
      name,
      ". The reports are installed from vignettes/ into doc/, ",
      "so the package has to be installed with its vignettes built (make install).",
      call. = FALSE
    )
  }
  path
}
