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

read_ptmsigdb <- function(ptmsigdb_file) {
  if (grepl("\\.rds$", ptmsigdb_file)) readRDS(ptmsigdb_file) else fgsea::gmtPathways(ptmsigdb_file)
}

# Human and mouse are merged rather than chosen between: the signatures are
# keyed on flanking sequence, not on organism, so a conserved site contributes
# the same sequence from either, and a signature curated in only one of them
# would otherwise be lost for the other.
.prepare_ptmsigdb_pathways <- function(pathways_human, pathways_mouse, keep_sources, trim_to) {
  all_names <- union(names(pathways_human), names(pathways_mouse))
  merged <- stats::setNames(
    lapply(all_names, function(name) unique(c(pathways_human[[name]], pathways_mouse[[name]]))),
    all_names
  )
  kept <- merged[grepl(paste0("^(", paste(keep_sources, collapse = "|"), ")_"), names(merged))]
  message(
    "PTMsigDB: kept ",
    length(kept),
    " of ",
    length(merged),
    " signatures (",
    paste(keep_sources, collapse = ", "),
    "), trimmed to ",
    trim_to,
    " residues"
  )
  trim_ptmsigdb_pathways(kept, trim_to)
}
