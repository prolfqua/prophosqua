## Rebuild the final two-contrast MuData artifact used by the statistics vignette.
##
##   Rscript data-raw/make_ptm_results_example.R

devtools::load_all(quiet = TRUE)

output <- file.path("inst", "extdata", "ptm_results_example.h5mu")
prophosqua:::example_ptm_results_h5mu(output)
message("OK: wrote ", output, " (", file.info(output)$size, " bytes)")
