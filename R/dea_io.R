# The file a prolfquapp DEA run writes under its Results_WU_<workunit> folder.
get_dea_file <- function(dea_dir, filename) {
  matches <- Sys.glob(file.path(dea_dir, "Results_WU_*", filename))
  if (length(matches) == 0) {
    stop("No ", filename, " found in: ", dea_dir, call. = FALSE)
  }
  matches[[1]]
}

# The total-proteome DEA names a protein by its FASTA id (sp|P12345|NAME), the
# site DEA by its accession. The mapping must stay one-to-one after decoy
# filtering, or two input identifiers would be joined as one protein.
canonicalize_uniprot_ids <- function(data, id_col = "protein_Id") {
  original <- as.character(data[[id_col]])
  canonical <- vapply(
    strsplit(original, "|", fixed = TRUE),
    function(parts) if (length(parts) >= 2) parts[[2]] else parts[[1]],
    character(1)
  )
  id_map <- unique(data.frame(original = original, canonical = canonical))
  ambiguous <- id_map$canonical[duplicated(id_map$canonical)]
  if (length(ambiguous) > 0) {
    stop("Protein identifier canonicalization is not one-to-one for: ", paste(unique(ambiguous), collapse = ", "))
  }
  data[[id_col]] <- canonical
  data
}
