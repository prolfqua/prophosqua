#!/usr/bin/env Rscript
# Build the PTM-SEA signature set from PTMsigDB.
#
# Downloads human and mouse, merges them, keeps the requested sub-sources
# and trims the flanking sequences.

suppressPackageStartupMessages({
  library(optparse)
  library(prophosqua)
})

option_list <- list(
  make_option(
    "--output_dir",
    type = "character",
    default = "data/ptmsigdb",
    help = "output directory for the filtered database"
  ),
  make_option(
    "--keep_sources",
    type = "character",
    default = "KINASE-PSP",
    help = "comma separated sub-sources to keep, e.g. KINASE-PSP,PATH-NP"
  ),
  make_option("--trim_to", type = "integer", default = 15, help = "trim flanking sequences to N residues")
)
opt <- parse_args(OptionParser(option_list = option_list))

prophosqua::prepare_ptmsigdb(
  output_dir = opt$output_dir,
  keep_sources = strsplit(opt$keep_sources, ",")[[1]],
  trim_to = opt$trim_to
)

message("Done.")
