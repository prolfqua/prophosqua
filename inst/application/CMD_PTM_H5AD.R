#!/usr/bin/env Rscript
# Compute DPA, DPU, and CorrectFirst and write one site-level PTM result H5AD.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
})

option_list <- list(
  make_option("--site_h5ad", type = "character", help = "site-level prolfquapp AnnData.h5ad"),
  make_option("--protein_h5ad", type = "character", help = "total-proteome prolfquapp AnnData.h5ad"),
  make_option("--annot_file", type = "character", help = "sample annotation defining groups and contrasts"),
  make_option("--output_h5ad", type = "character", help = "output path for PTM_results.h5ad")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (required in c("site_h5ad", "protein_h5ad", "annot_file", "output_h5ad")) {
  if (is.null(opt[[required]])) {
    stop("--", required, " is required", call. = FALSE)
  }
}

dir.create(dirname(opt$output_h5ad), recursive = TRUE, showWarnings = FALSE)
message("Writing ", opt$output_h5ad)
prophosqua::compute_ptm_results_h5ad(
  site_h5ad = opt$site_h5ad,
  protein_h5ad = opt$protein_h5ad,
  annot_file = opt$annot_file,
  output_h5ad = opt$output_h5ad
)
message("Done.")
