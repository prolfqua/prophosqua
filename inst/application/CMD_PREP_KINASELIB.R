#!/usr/bin/env Rscript
# Write the kinase-library inputs: sequence windows and ranked files.
#
# One list of unique sequence windows for motif scanning, and one ranked file
# per contrast for motif enrichment.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
})

option_list <- list(
  make_option("--xlsx_file", type = "character", help = "combined PTM_results.xlsx"),
  make_option("--output_dir", type = "character", help = "output directory, normally <analysis>/KinaseLib"),
  make_option("--analysis_type", type = "character", help = "DPA, DPU or CF; names the output files"),
  make_option("--sheet", type = "character", help = "sheet to read"),
  make_option("--stat_column", type = "character", help = "column to rank the sites on")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("xlsx_file", "output_dir", "analysis_type", "sheet", "stat_column")
for (name in required) {
  if (is.null(opt[[name]])) {
    stop("--", name, " is required", call. = FALSE)
  }
}

prophosqua::prep_kinaselib_inputs(
  xlsx_file = opt$xlsx_file,
  output_dir = opt$output_dir,
  analysis_type = opt$analysis_type,
  sheet = opt$sheet,
  stat_column = opt$stat_column
)
