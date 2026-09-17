#!/usr/bin/env Rscript
# Write terminal delivery workbooks and RDS files from final MuData.
suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--input", type = "character", help = "final PTM_results H5MU"),
    make_option("--output_dir", type = "character", help = "terminal delivery directory")
  )
))
stopifnot(!is.null(opt$input), !is.null(opt$output_dir))
prophosqua::export_ptm_h5mu(opt$input, opt$output_dir)
