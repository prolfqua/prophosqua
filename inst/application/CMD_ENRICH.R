#!/usr/bin/env Rscript
# Compute one enrichment result, or the kinase-library inputs, of one analysis.
suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--statistics", type = "character", help = "shared PTM_statistics.h5mu"),
    make_option("--output", type = "character", help = "target file"),
    make_option("--stage", type = "character", help = "PTMSEA, KinaseInputs or KinaseGSEA"),
    make_option("--analysis", type = "character", help = "DPA, DPU or CF"),
    make_option("--kinase_inputs", type = "character", help = "KinaseInputs file (json.gz)"),
    make_option("--kinase_assignments", type = "character", help = "KinaseAssignments file (json.gz)")
  )
))
stopifnot(!is.null(opt$statistics), !is.null(opt$output), !is.null(opt$stage), !is.null(opt$analysis))
preparation <- list()
if (!is.null(opt$kinase_inputs)) {
  preparation$KinaseInputs <- opt$kinase_inputs
}
if (!is.null(opt$kinase_assignments)) {
  preparation$KinaseAssignments <- opt$kinase_assignments
}
prophosqua::compute_ptm_enrichment(
  opt$statistics,
  opt$output,
  opt$stage,
  opt$analysis,
  preparation
)
