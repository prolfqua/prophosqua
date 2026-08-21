#!/usr/bin/env Rscript
# Combine every analysis into the one workbook downstream reports read.
#
# The DPA, DPU and CorrectFirst results, plus the normalized abundances.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
})

option_list <- list(
  make_option("--dpa_xlsx", type = "character", help = "Result_DPA.xlsx"),
  make_option("--dpu_xlsx", type = "character", help = "Result_DPU.xlsx"),
  make_option("--cf_xlsx", type = "character", help = "CorrectFirst_PTM_usage_results.xlsx"),
  make_option("--protein_parquet", type = "character", help = "normalized protein abundances"),
  make_option("--site_parquet", type = "character", help = "normalized site abundances"),
  make_option("--output_xlsx", type = "character", help = "combined workbook"),
  make_option("--output_rds", type = "character", help = "combined RDS")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c(
  "dpa_xlsx",
  "dpu_xlsx",
  "cf_xlsx",
  "protein_parquet",
  "site_parquet",
  "output_xlsx",
  "output_rds"
)
for (name in required) {
  if (is.null(opt[[name]])) {
    stop("--", name, " is required", call. = FALSE)
  }
}

prophosqua::combine_ptm_results(
  dpa_xlsx = opt$dpa_xlsx,
  dpu_xlsx = opt$dpu_xlsx,
  cf_xlsx = opt$cf_xlsx,
  protein_parquet = opt$protein_parquet,
  site_parquet = opt$site_parquet,
  output_xlsx = opt$output_xlsx,
  output_rds = opt$output_rds
)

message("Done.")
