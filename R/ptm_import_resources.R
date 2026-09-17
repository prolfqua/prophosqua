# Import reads external resources once; every downstream stage reads MuData.
.import_ptm_resources <- function(parameters, ptmsigdb_file = NULL) {
  if (!isTRUE(parameters$run_kinase)) {
    return(list())
  }
  if (!is.null(ptmsigdb_file)) {
    return(list(
      ptmsigdb = read_ptmsigdb(ptmsigdb_file),
      ptmsigdb_source = list(path = normalizePath(ptmsigdb_file), md5 = unname(tools::md5sum(ptmsigdb_file)))
    ))
  }
  cache <- tempfile("ptmsigdb-import-")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  sources <- vapply(c("human", "mouse"), function(species) download_ptmsigdb(species, output_dir = cache), character(1))
  pathways <- lapply(sources, fgsea::gmtPathways)
  settings <- parameters$ptmsigdb
  list(
    ptmsigdb = .prepare_ptmsigdb_pathways(pathways$human, pathways$mouse, settings$keep_sources, settings$trim_to),
    ptmsigdb_source = list(files = basename(sources), md5 = unname(tools::md5sum(sources)))
  )
}
