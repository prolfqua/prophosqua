#!/usr/bin/env Rscript
# Compute complete DPA, DPU and CorrectFirst statistics through MuData.
suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--input", type = "character", help = "paired-input H5MU"),
    make_option("--output", type = "character", help = "complete statistics H5MU")
  )
))
stopifnot(!is.null(opt$input), !is.null(opt$output))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
prophosqua::compute_ptm_results_h5mu(opt$input, opt$output)
