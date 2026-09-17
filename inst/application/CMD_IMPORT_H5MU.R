#!/usr/bin/env Rscript
# Import paired DEA experiments, settings and reference data into MuData.
suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--enriched", type = "character", help = "enriched DEA H5AD"),
    make_option("--total", type = "character", help = "total DEA H5AD"),
    make_option("--config_json", type = "character", help = "merged pipeline configuration as JSON"),
    make_option(
      "--ptmsigdb",
      type = "character",
      default = NULL,
      help = "existing filtered reference RDS or GMT to import"
    ),
    make_option("--output", type = "character", help = "paired-input H5MU")
  )
))
stopifnot(!is.null(opt$enriched), !is.null(opt$total), !is.null(opt$config_json), !is.null(opt$output))
parameters <- jsonlite::fromJSON(opt$config_json, simplifyVector = FALSE)
resources <- prophosqua:::.import_ptm_resources(parameters, opt$ptmsigdb)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
prophosqua::import_ptm_h5mu(opt$enriched, opt$total, opt$output, resources, parameters)
