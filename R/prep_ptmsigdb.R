#' Download, Merge and Filter PTMsigDB
#'
#' Fetches PTMsigDB for human and mouse, merges the two into one signature set,
#' keeps only the requested sub-sources, and trims the flanking sequences to the
#' width the site records of an analysis use. The result is the signature
#' database PTM-SEA scores against.
#'
#' Human and mouse are merged rather than chosen between because the signatures
#' are keyed on flanking sequence, not on organism: a site conserved between the
#' two contributes the same sequence from either, and a signature curated in
#' only one of them would otherwise be lost for the other.
#'
#' Sub-sources present in PTMsigDB v2.0.0 are `KINASE-PSP` (PhosphoSitePlus
#' kinase-substrate, experimental), `KINASE-iKiP` (iKiP, computational),
#' `PATH-NP` / `PATH-WP` / `PATH-BI` (NetPath, WikiPathways, Broad Institute
#' pathways), `PERT-PSP` and `PERT-P100-*` (perturbations), and `DISEASE-PSP`.
#'
#' @param output_dir Directory for the filtered database and the download cache.
#' @param keep_sources Sub-sources to keep, as a character vector.
#' @param trim_to Flanking width to trim the sequences to; one of 11, 13 or 15.
#' @return Invisibly, a list with the `rds` and `gmt` paths written and the
#'   `pathways` themselves.
#' @export
#' @examples
#' # Downloads from PhosphoSitePlus; needs network access.
#' \dontrun{
#' prepare_ptmsigdb("data/ptmsigdb", keep_sources = "KINASE-PSP", trim_to = 15)
#' }
prepare_ptmsigdb <- function(output_dir = "data/ptmsigdb", keep_sources = "KINASE-PSP", trim_to = 15) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea is required to read the PTMsigDB GMT files.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Downloading PTMsigDB for human and mouse...")
  cache_dir <- file.path(output_dir, "raw_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  pathways_human <- fgsea::gmtPathways(
    download_ptmsigdb(species = "human", output_dir = cache_dir)
  )
  pathways_mouse <- fgsea::gmtPathways(
    download_ptmsigdb(species = "mouse", output_dir = cache_dir)
  )
  message(sprintf("Human pathways: %d", length(pathways_human)))
  message(sprintf("Mouse pathways: %d", length(pathways_mouse)))

  message("Merging human and mouse pathways...")
  all_names <- union(names(pathways_human), names(pathways_mouse))
  pathways_merged <- lapply(all_names, function(name) {
    unique(c(pathways_human[[name]], pathways_mouse[[name]]))
  })
  names(pathways_merged) <- all_names
  message(sprintf("Total pathways before filtering: %d", length(pathways_merged)))

  keep_pattern <- paste0("^(", paste(keep_sources, collapse = "|"), ")_")
  pathways_filtered <- pathways_merged[grepl(keep_pattern, names(pathways_merged))]
  message(sprintf(
    "Pathways after filtering (keeping %s): %d",
    paste(keep_sources, collapse = ","),
    length(pathways_filtered)
  ))

  report_ptmsigdb_sources(pathways_merged, pathways_filtered)

  message(sprintf("Trimming flanking sequences to %d-mer...", trim_to))
  pathways_trimmed <- trim_ptmsigdb_pathways(
    pathways_filtered,
    trim_to = as.character(trim_to)
  )

  stem <- sprintf(
    "ptmsigdb_filtered_%s_%dmer",
    paste(keep_sources, collapse = "_"),
    trim_to
  )
  gmt_file <- file.path(output_dir, paste0(stem, ".gmt"))
  rds_file <- file.path(output_dir, paste0(stem, ".rds"))

  write_gmt(pathways_trimmed, gmt_file)
  message(sprintf("Wrote filtered GMT to: %s", gmt_file))

  saveRDS(pathways_trimmed, rds_file)
  message(sprintf("Wrote RDS cache to: %s", rds_file))

  invisible(list(rds = rds_file, gmt = gmt_file, pathways = pathways_trimmed))
}

#' Write a Pathway List as GMT
#'
#' GMT is one line per set: the set name, a description field, then its members,
#' all tab separated. PTMsigDB carries no per-set description, so the field is
#' written as `NA` — the placeholder the format expects rather than an empty
#' column, which some readers would collapse.
#'
#' @param pathways Named list of character vectors.
#' @param file Path to write.
#' @return Invisibly, `file`.
#' @export
#' @examples
#' gmt <- tempfile(fileext = ".gmt")
#' write_gmt(list(KINASE_A = c("SITE1", "SITE2")), gmt)
#' readLines(gmt)
write_gmt <- function(pathways, file) {
  lines <- vapply(
    names(pathways),
    function(name) paste(c(name, "NA", pathways[[name]]), collapse = "\t"),
    character(1)
  )
  writeLines(lines, file)
  invisible(file)
}

#' Report the Sub-Source Breakdown of a PTMsigDB Filter
#'
#' Prints how many signatures each sub-source contributed and how many survived
#' the filter, so a run's log says what was kept without the caller having to
#' inspect the database.
#'
#' @param merged All merged pathways.
#' @param filtered The pathways that survived filtering.
#' @return Invisibly, `NULL`.
#' @keywords internal
report_ptmsigdb_sources <- function(merged, filtered) {
  all_sources <- c(
    "KINASE-PSP",
    "KINASE-iKiP",
    "PATH-NP",
    "PATH-WP",
    "PATH-BI",
    "PERT-PSP",
    "PERT-P100-DIA2",
    "PERT-P100-PRM",
    "PERT-P100-DIA",
    "DISEASE-PSP"
  )
  message("Sub-source breakdown:")
  for (src in all_sources) {
    prefix <- paste0("^", src, "_")
    n <- sum(grepl(prefix, names(merged)))
    if (n == 0) {
      next
    }
    kept <- sum(grepl(prefix, names(filtered)))
    message(sprintf(
      "  %s: %d total, %d kept%s",
      src,
      n,
      kept,
      if (kept > 0) " *" else ""
    ))
  }
  invisible(NULL)
}
