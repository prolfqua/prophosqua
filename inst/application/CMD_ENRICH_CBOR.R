#!/usr/bin/env Rscript
# Compute one enrichment result or preparation as a compact CBOR handoff.
suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--statistics", type = "character", help = "shared PTM_statistics.h5mu"),
    make_option("--output", type = "character", help = "target CBOR artifact"),
    make_option("--stage", type = "character", help = "PTMSEA, KinaseInputs, KinaseGSEA or MEA"),
    make_option("--analysis", type = "character", help = "DPA, DPU or CF"),
    make_option("--kinase_inputs", type = "character", help = "KinaseInputs CBOR"),
    make_option("--kinase_assignments", type = "character", help = "KinaseAssignments CBOR"),
    make_option("--motif_enrichment", type = "character", help = "MotifEnrichment CBOR")
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
if (!is.null(opt$motif_enrichment)) {
  preparation$MotifEnrichment <- opt$motif_enrichment
}
prophosqua::compute_ptm_enrichment_cbor(
  opt$statistics,
  opt$output,
  opt$stage,
  opt$analysis,
  preparation
)
