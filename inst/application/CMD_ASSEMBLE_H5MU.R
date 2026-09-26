#!/usr/bin/env Rscript
# Assemble the final MuData from the statistics and the enrichment files beside it.
suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--statistics", type = "character", help = "statistics H5MU"),
    make_option("--output", type = "character", help = "final PTM_results H5MU, in the directory of the analyses")
  )
))
stopifnot(!is.null(opt$statistics), !is.null(opt$output))
prophosqua::assemble_ptm_results(opt$statistics, opt$output)
