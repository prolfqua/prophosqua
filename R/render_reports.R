#' Render a Report Shipped with prophosqua
#'
#' Renders one of the R Markdown documents under `inst/application` from where
#' it is installed. A project therefore never carries a copy of a report
#' template, and cannot end up rendering an edited copy while the package holds
#' a newer one.
#'
#' `knit_root_dir` is the working directory, not the template's directory, so
#' relative paths in `params` mean what they mean to the caller.
#'
#' The output format is the one the report's own YAML declares. Forcing one here
#' would be wrong for the reports that ask for something specific: the index page
#' declares its own theme, and only some reports want bookdown's numbered
#' sections and cross-references.
#'
#' Each render is given its own `intermediates_dir` because knitr names its
#' intermediate files after the input document: two renders of the same template
#' running at once would otherwise overwrite each other's `.knit.md` and both
#' reports would end up with whichever content finished last.
#'
#' @param name File name of the report, e.g. `"Analysis_seqlogo.Rmd"`.
#' @param output_file File name to write, e.g. `"Analysis_seqlogo.html"`.
#' @param output_dir Directory to write the report to.
#' @param params Named list passed to the report's `params`.
#' @param intermediates_dir Directory for knitr's intermediates. Defaults to a
#'   private directory of this R process, removed when it exits.
#' @return Invisibly, the path of the rendered file.
#' @export
#' @examples
#' # Renders one of the installed templates; needs the data it asks for.
#' \dontrun{
#' render_ptm_report(
#'   "Analysis_seqlogo.Rmd", "Analysis_seqlogo.html", "PTM_DPA",
#'   params = list(xlsx_file = "PTM_results.xlsx", sheet = "DPA")
#' )
#' }
render_ptm_report <- function(name, output_file, output_dir, params = list(), intermediates_dir = NULL) {
  rmd_path <- application_file(name)

  if (is.null(intermediates_dir)) {
    intermediates_dir <- file.path(tempdir(), paste0("render_", basename(name)))
  }
  dir.create(intermediates_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rmarkdown::render(
    rmd_path,
    output_file = output_file,
    output_dir = output_dir,
    knit_root_dir = getwd(),
    intermediates_dir = intermediates_dir,
    params = params,
    envir = new.env(parent = globalenv())
  )

  invisible(file.path(output_dir, output_file))
}

#' Render the Phospho and Protein Integration Overview
#'
#' Renders the integration overview from the DPU result object. This report
#' takes R objects rather than file paths as parameters, so it cannot go through
#' [render_ptm_report()]: the objects are built here from the saved result and
#' handed over directly.
#'
#' The template and its bibliography are copied into a private render directory
#' first. Rendering in place would write knitr's intermediates and both copied
#' files into the project root, and would let two concurrent renders overwrite
#' each other's intermediates.
#'
#' @param input_rds The `combined_test_diff.rds` written by the DPU computation.
#' @param output_dir Directory to write `Result_DPU.html` to.
#' @param project_id Project identifier shown in the report header.
#' @param work_unit_id Work unit identifier shown in the report header.
#' @return Invisibly, the path of the rendered file.
#' @export
#' @examples
#' # Needs the DPU result of a pipeline run.
#' \dontrun{
#' render_dpu_overview("PTM_DPU/combined_test_diff.rds", "PTM_DPU")
#' }
render_dpu_overview <- function(input_rds, output_dir, project_id = "PTM_analysis", work_unit_id = "DPU_Integration") {
  if (!file.exists(input_rds)) {
    stop("Input RDS file not found: ", input_rds, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading data from: ", input_rds)
  combined_test_diff <- readRDS(input_rds)

  grp <- prolfquapp::make_DEA_config_R6(
    PROJECTID = project_id,
    ORDERID = "fgcz_project",
    WORKUNITID = work_unit_id
  )

  render_dir <- file.path(tempdir(), "dpu_overview")
  dir.create(render_dir, recursive = TRUE, showWarnings = FALSE)

  template <- "_Overview_PhosphoAndIntegration_site.Rmd"
  for (file in c(template, "bibliography2025.bib")) {
    file.copy(application_file(file), file.path(render_dir, file), overwrite = TRUE)
  }

  message("Rendering DPU overview report...")
  rmarkdown::render(
    file.path(render_dir, template),
    knit_root_dir = getwd(),
    params = list(data = combined_test_diff, grp = grp),
    output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE),
    envir = new.env(parent = globalenv())
  )

  output_file <- file.path(output_dir, "Result_DPU.html")
  file.copy(
    from = file.path(render_dir, sub("\\.Rmd$", ".html", template)),
    to = output_file,
    overwrite = TRUE
  )
  message("DPU overview report saved to: ", output_file)

  invisible(output_file)
}

#' Resolve a File Shipped under inst/application
#'
#' @param name File name.
#' @return Full path to the installed file.
#' @keywords internal
application_file <- function(name) {
  path <- system.file("application", name, package = "prophosqua")
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "prophosqua application file not found: ",
      name,
      ". Reinstall the package.",
      call. = FALSE
    )
  }
  path
}
