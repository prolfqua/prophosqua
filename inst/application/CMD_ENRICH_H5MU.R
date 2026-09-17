#!/usr/bin/env Rscript
# Build a complete preparation or enrichment stage using only MuData.
suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--input", type = "character", help = "source stage H5MU"),
    make_option("--output", type = "character", help = "completed target H5MU"),
    make_option("--stage", type = "character", help = "PTMSEA, KinaseInputs, KinaseGSEA or MEA"),
    make_option("--analysis", type = "character", help = "DPA, DPU or CF")
  )
))
stopifnot(!is.null(opt$input), !is.null(opt$output), !is.null(opt$analysis))
types <- list(
  PTMSEA = prophosqua::PTMSEA,
  KinaseInputs = prophosqua::KinaseInputs,
  KinaseGSEA = prophosqua::KinaseGSEA,
  MEA = prophosqua::MEA
)
Type <- types[[opt$stage]]
if (is.null(Type)) {
  stop("Unknown target stage: ", opt$stage)
}
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
prophosqua::read_ptm_h5mu(opt$input)$build(Type, analysis = opt$analysis)$write_h5mu(opt$output)
