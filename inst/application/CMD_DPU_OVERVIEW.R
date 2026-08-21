#!/usr/bin/env Rscript
# Render the phospho and protein integration overview from the DPU result.
#
# This report takes R objects rather than file paths as parameters, so it has
# its own entry point instead of going through CMD_RENDER.R.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
})

option_list <- list(
  make_option("--input_rds", type = "character", help = "combined_test_diff.rds written by the DPU computation"),
  make_option("--output_dir", type = "character", help = "directory to write Result_DPU.html to"),
  make_option(
    "--project_id",
    type = "character",
    default = "PTM_analysis",
    help = "project identifier shown in the report header"
  ),
  make_option(
    "--work_unit_id",
    type = "character",
    default = "DPU_Integration",
    help = "work unit identifier shown in the report header"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

for (required in c("input_rds", "output_dir")) {
  if (is.null(opt[[required]])) {
    stop("--", required, " is required", call. = FALSE)
  }
}

prophosqua::render_dpu_overview(
  input_rds = opt$input_rds,
  output_dir = opt$output_dir,
  project_id = opt$project_id,
  work_unit_id = opt$work_unit_id
)
