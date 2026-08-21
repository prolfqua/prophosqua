#!/usr/bin/env Rscript
# Compute PTM-SEA against PTMsigDB.
#
# Data only. The workbook, the result objects and the summary tables the report
# shows are written here, so that fixing a caption costs a render and not a set
# of permutation tests. An enrichment that found nothing still writes both
# files, with zero rows: that is what lets a rule declare them as outputs.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
  library(dplyr)
})

option_list <- list(
  make_option("--xlsx_file", type = "character", help = "combined PTM_results.xlsx"),
  make_option("--sheet", type = "character", help = "sheet to read"),
  make_option("--stat_column", type = "character", default = "statistic.site", help = "per-site statistic to rank on"),
  make_option("--ptmsigdb_file", type = "character", help = "filtered PTMsigDB, .rds or .gmt"),
  make_option("--analysis_type", type = "character", help = "DPA, DPU or CF; names the output files"),
  make_option("--output_dir", type = "character", help = "output directory"),
  make_option("--trim_to", type = "integer", default = 15, help = "flanking width the database uses"),
  make_option("--min_size", type = "integer", default = 10, help = "smallest testable signature"),
  make_option("--max_size", type = "integer", default = 500, help = "largest testable signature"),
  make_option("--n_perm", type = "integer", default = 1000, help = "permutations of the enrichment test")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("xlsx_file", "sheet", "ptmsigdb_file", "analysis_type", "output_dir")
for (name in required) {
  if (is.null(opt[[name]])) {
    stop("--", name, " is required", call. = FALSE)
  }
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

res <- prophosqua::compute_ptmsea(
  xlsx_file = opt$xlsx_file,
  sheet = opt$sheet,
  stat_column = opt$stat_column,
  ptmsigdb_file = opt$ptmsigdb_file,
  trim_to = opt$trim_to,
  min_size = opt$min_size,
  max_size = opt$max_size,
  n_perm = opt$n_perm
)

stem <- file.path(opt$output_dir, paste0("PTMSEA_", opt$analysis_type, "_results"))

export_data <- res$all_clean |> arrange(contrast, pvalue)
export_list <- list(all_clean = export_data)
for (ct in unique(export_data$contrast)) {
  sheet_name <- gsub("[^a-zA-Z0-9_]", "_", substr(ct, 1, 31))
  export_list[[sheet_name]] <- export_data |> filter(contrast == ct)
}
export_list[["significant_FDR10"]] <- export_data |> filter(p.adjust < 0.1)

message("Writing ", paste0(stem, ".xlsx"))
writexl::write_xlsx(export_list, paste0(stem, ".xlsx"))

message("Writing ", paste0(stem, ".rds"))
saveRDS(res, paste0(stem, ".rds"))

message("Done.")
