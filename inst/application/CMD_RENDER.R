#!/usr/bin/env Rscript
# Render one of the reports shipped under inst/application.
#
# The report parameters are given as trailing key=value arguments:
#
#   CMD_RENDER.R --report Analysis_seqlogo.Rmd \
#       --output_file Analysis_seqlogo.html --output_dir PTM_DPA \
#       xlsx_file=PTM_results.xlsx sheet=DPA fdr=0.25
#
# A value that reads as a number becomes numeric and TRUE/FALSE becomes
# logical, because a report declaring `max_fig: 10` must not be handed the
# string "10". Everything else stays a string.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
})

option_list <- list(
  make_option("--report", type = "character", help = "report file name, e.g. Analysis_seqlogo.Rmd"),
  make_option("--output_file", type = "character", help = "file name to write, e.g. Analysis_seqlogo.html"),
  make_option("--output_dir", type = "character", help = "directory to write the report to"),
  make_option("--intermediates_dir", type = "character", default = NULL, help = "directory for knitr intermediates")
)
parsed <- parse_args(
  OptionParser(option_list = option_list, usage = "%prog [options] key=value ..."),
  positional_arguments = TRUE
)
opt <- parsed$options

for (required in c("report", "output_file", "output_dir")) {
  if (is.null(opt[[required]])) {
    stop("--", required, " is required", call. = FALSE)
  }
}

coerce_param <- function(value) {
  if (value %in% c("TRUE", "FALSE")) {
    return(as.logical(value))
  }
  if (grepl("^-?[0-9]+$", value)) {
    return(as.integer(value))
  }
  if (grepl("^-?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$", value)) {
    return(as.numeric(value))
  }
  value
}

report_params <- list()
for (arg in parsed$args) {
  if (!grepl("=", arg, fixed = TRUE)) {
    stop("report parameter is not key=value: ", arg, call. = FALSE)
  }
  key <- sub("=.*$", "", arg)
  report_params[[key]] <- coerce_param(sub("^[^=]*=", "", arg))
}

prophosqua::render_ptm_report(
  name = opt$report,
  output_file = opt$output_file,
  output_dir = opt$output_dir,
  params = report_params,
  intermediates_dir = opt$intermediates_dir
)
