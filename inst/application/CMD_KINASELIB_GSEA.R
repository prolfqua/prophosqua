#!/usr/bin/env Rscript
# Compute the kinase-library GSEA.
#
# Data only, and always written: an enrichment that found nothing produces a
# workbook with zero rows rather than no workbook, so a rule can declare it.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
  library(dplyr)
})

option_list <- list(
  make_option("--xlsx_file", type = "character", help = "combined PTM_results.xlsx"),
  make_option("--sheet", type = "character", help = "sheet to read"),
  make_option("--term2gene_file", type = "character", help = "kinase-to-site assignment from the motif scan"),
  make_option("--analysis_type", type = "character", help = "DPA, DPU or CF; names the output files"),
  make_option("--output_dir", type = "character", help = "output directory"),
  make_option("--min_size", type = "integer", default = 15, help = "smallest testable kinase set"),
  make_option("--max_size", type = "integer", default = 5000, help = "largest testable kinase set"),
  make_option("--n_perm", type = "integer", default = 1000, help = "permutations of the enrichment test")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("xlsx_file", "sheet", "term2gene_file", "analysis_type", "output_dir")
for (name in required) {
  if (is.null(opt[[name]])) {
    stop("--", name, " is required", call. = FALSE)
  }
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

res <- prophosqua::compute_kinaselib_gsea(
  xlsx_file = opt$xlsx_file,
  sheet = opt$sheet,
  term2gene_file = opt$term2gene_file,
  min_size = opt$min_size,
  max_size = opt$max_size,
  n_perm = opt$n_perm
)

stem <- file.path(opt$output_dir, paste0("KinaseLib_GSEA_", opt$analysis_type))

export_list <- list(
  all_results = res$all_results |> arrange(contrast, FDR),
  significant = res$all_results |> filter(FDR < 0.1) |> arrange(contrast, FDR),
  summary = res$gsea_info
)

message("Writing ", paste0(stem, ".xlsx"))
writexl::write_xlsx(export_list, paste0(stem, ".xlsx"))

message("Writing ", paste0(stem, ".rds"))
saveRDS(res, paste0(stem, ".rds"))

message("Done.")
