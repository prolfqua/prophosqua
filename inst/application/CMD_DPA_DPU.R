#!/usr/bin/env Rscript
# Compute Differential PTM Abundance and Differential PTM Usage.
#
# Data only: the two result workbooks, the DPU object the integration overview
# renders from, and the small object the report quotes its numbers out of. The
# report itself is a separate step.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
})

option_list <- list(
  make_option("--phospho_dea_dir", type = "character", help = "phospho DEA output directory"),
  make_option("--protein_dea_dir", type = "character", help = "total-proteome DEA output directory"),
  make_option("--site_h5ad", type = "character", help = "site-level prolfquapp AnnData.h5ad"),
  make_option("--protein_h5ad", type = "character", help = "total-proteome prolfquapp AnnData.h5ad"),
  make_option("--dpa_dir", type = "character", help = "output directory for the DPA results"),
  make_option("--dpu_dir", type = "character", help = "output directory for the DPU results")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (required in c("dpa_dir", "dpu_dir")) {
  if (is.null(opt[[required]])) {
    stop("--", required, " is required", call. = FALSE)
  }
}

h5ad_any <- !is.null(opt$site_h5ad) || !is.null(opt$protein_h5ad)
h5ad_pair <- !is.null(opt$site_h5ad) && !is.null(opt$protein_h5ad)
dea_dir_any <- !is.null(opt$phospho_dea_dir) || !is.null(opt$protein_dea_dir)
dea_dir_pair <- !is.null(opt$phospho_dea_dir) && !is.null(opt$protein_dea_dir)
valid_h5ad_input <- h5ad_pair && !dea_dir_any
valid_dea_dir_input <- dea_dir_pair && !h5ad_any
if (!valid_h5ad_input && !valid_dea_dir_input) {
  stop(
    "Provide exactly one complete input pair: --site_h5ad and --protein_h5ad, ",
    "or --phospho_dea_dir and --protein_dea_dir.",
    call. = FALSE
  )
}

dir.create(opt$dpa_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$dpu_dir, recursive = TRUE, showWarnings = FALSE)

if (h5ad_pair) {
  res <- prophosqua::compute_dpa_dpu_h5ad(
    site_h5ad = opt$site_h5ad,
    protein_h5ad = opt$protein_h5ad
  )
  phospho_input <- dirname(opt$site_h5ad)
  protein_input <- dirname(opt$protein_h5ad)
} else {
  res <- prophosqua::compute_dpa_dpu(
    phospho_dea_dir = opt$phospho_dea_dir,
    protein_dea_dir = opt$protein_dea_dir
  )
  phospho_input <- opt$phospho_dea_dir
  protein_input <- opt$protein_dea_dir
}

dpa_xlsx <- file.path(opt$dpa_dir, "Result_DPA.xlsx")
dpu_xlsx <- file.path(opt$dpu_dir, "Result_DPU.xlsx")

message("Writing ", dpa_xlsx)
writexl::write_xlsx(
  list(combinedSiteProteinData = res$combined_site_prot),
  path = dpa_xlsx
)

message("Writing ", dpu_xlsx)
writexl::write_xlsx(
  list(combinedStats = res$combined_test_diff),
  path = dpu_xlsx
)

message("Writing ", file.path(opt$dpu_dir, "combined_test_diff.rds"))
saveRDS(res$combined_test_diff, file.path(opt$dpu_dir, "combined_test_diff.rds"))

# Everything the report quotes, so that rendering it needs neither of the two
# large workbooks nor a second pass over the DEA output.
message("Writing ", file.path(opt$dpa_dir, "dpa_dpu_objects.rds"))
report_objects <- list(
  match_rates = res$match_rates,
  n_dpa_rows = nrow(res$combined_site_prot),
  n_dpu_rows = nrow(res$combined_test_diff),
  phospho_dea_dir = phospho_input,
  protein_dea_dir = protein_input,
  dpa_xlsx = dpa_xlsx,
  dpu_xlsx = dpu_xlsx
)
if (h5ad_pair) {
  report_objects$site_h5ad <- opt$site_h5ad
  report_objects$protein_h5ad <- opt$protein_h5ad
}
saveRDS(report_objects, file.path(opt$dpa_dir, "dpa_dpu_objects.rds"))

message("Done.")
