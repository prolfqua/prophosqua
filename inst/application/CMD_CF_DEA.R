#!/usr/bin/env Rscript
# Compute PTM usage the CorrectFirst way: correct, then model.
#
# Site abundances are corrected by their protein abundance before modelling.
#
# Data only. cf_objects.rds carries the fitted contrasts and the counts the
# report quotes, so a change to the report's prose or captions costs a render
# and not a refit.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
})

option_list <- list(
  make_option("--phospho_dea_dir", type = "character", help = "phospho DEA output directory"),
  make_option("--protein_dea_dir", type = "character", help = "total-proteome DEA output directory"),
  make_option("--annot_file", type = "character", help = "sample annotation defining groups and contrasts"),
  make_option("--output_dir", type = "character", help = "output directory for the CorrectFirst results")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (required in c("phospho_dea_dir", "protein_dea_dir", "annot_file", "output_dir")) {
  if (is.null(opt[[required]])) {
    stop("--", required, " is required", call. = FALSE)
  }
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

res <- prophosqua::compute_cf_dea(
  phospho_dea_dir = opt$phospho_dea_dir,
  protein_dea_dir = opt$protein_dea_dir,
  annot_file = opt$annot_file
)

results_xlsx <- file.path(opt$output_dir, "CorrectFirst_PTM_usage_results.xlsx")

message("Writing ", results_xlsx)
writexl::write_xlsx(list(results = res$results), path = results_xlsx)

message("Writing the corrected intensity matrix")
writexl::write_xlsx(
  res$wide_data,
  file.path(opt$output_dir, "CorrectFirst_intensities_wide.xlsx")
)
writexl::write_xlsx(
  res$wide_annotation,
  file.path(opt$output_dir, "CorrectFirst_intensities_file_annotation.xlsx")
)

# The fitted objects the report plots from, without the two intensity tables
# that are already on disk as workbooks.
message("Writing ", file.path(opt$output_dir, "cf_objects.rds"))
saveRDS(
  res[setdiff(names(res), c("wide_data", "wide_annotation"))],
  file.path(opt$output_dir, "cf_objects.rds")
)

message("Done.")
