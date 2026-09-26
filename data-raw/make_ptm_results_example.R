## Rebuild the example PTM results the vignettes read: the final MuData and,
## beside it, the enrichment files of each analysis.
##
##   Rscript data-raw/make_ptm_results_example.R

devtools::load_all(quiet = TRUE)

root <- file.path("inst", "extdata", "ptm_results_example")
unlink(root, recursive = TRUE)
path <- prophosqua:::example_ptm_results(root)
files <- list.files(root, recursive = TRUE, full.names = TRUE)
message("OK: wrote ", root, ", ", length(files), " files, ", sum(file.info(files)$size), " bytes")
