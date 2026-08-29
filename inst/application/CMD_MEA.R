#!/usr/bin/env Rscript
# Collect the per-contrast motif enrichment results into one workbook.
#
# Data only. The enrichment itself already happened, one contrast per
# mea_*.csv; this brings them onto shared column names and writes the result
# the report renders from.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
  library(dplyr)
})

option_list <- list(
  make_option("--kinaselib_dir", type = "character", help = "directory holding the mea_*.csv files"),
  make_option("--analysis_type", type = "character", help = "DPA, DPU or CF; names the output files")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (required in c("kinaselib_dir", "analysis_type")) {
  if (is.null(opt[[required]])) {
    stop("--", required, " is required", call. = FALSE)
  }
}

res <- prophosqua::compute_mea(kinaselib_dir = opt$kinaselib_dir)

stem <- file.path(
  opt$kinaselib_dir,
  paste0("MEA_", opt$analysis_type, "_results")
)

export_cols <- c(
  "contrast",
  "kinase",
  "NES",
  "pvalue",
  "FDR",
  "n_leading",
  "set_size"
)
export_list <- list(
  all_results = res$mea_clean |>
    select(all_of(export_cols)) |>
    arrange(contrast, FDR),
  significant = res$mea_clean |>
    filter(FDR < 0.1) |>
    select(all_of(export_cols)) |>
    arrange(contrast, FDR),
  summary = res$summary_df
)

message("Writing ", paste0(stem, ".xlsx"))
writexl::write_xlsx(export_list, paste0(stem, ".xlsx"))

message("Writing ", paste0(stem, ".rds"))
saveRDS(res, paste0(stem, ".rds"))

message("Writing ", paste0(stem, ".json"))
ranks <- prophosqua::read_mea_ranks(opt$kinaselib_dir)
term2gene <- utils::read.csv(
  file.path(opt$kinaselib_dir, "term2gene.csv"),
  stringsAsFactors = FALSE
)
prophosqua::write_gsea_result_json(
  prophosqua::mea_gsea_result_data(res$mea_clean, ranks, term2gene),
  paste0(stem, ".json")
)

message("Done.")
