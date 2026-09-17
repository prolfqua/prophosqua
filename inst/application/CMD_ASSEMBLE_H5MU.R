#!/usr/bin/env Rscript
# Assemble complete statistics and enabled enrichment MuData artifacts.
suppressPackageStartupMessages(library(optparse))
parsed <- parse_args(
  OptionParser(
    option_list = list(
      make_option("--statistics", type = "character", help = "statistics H5MU"),
      make_option("--output", type = "character", help = "final PTM_results H5MU")
    )
  ),
  positional_arguments = TRUE
)
opt <- parsed$options
stopifnot(!is.null(opt$statistics), !is.null(opt$output))
prophosqua::assemble_ptm_h5mu(opt$statistics, parsed$args, opt$output)
