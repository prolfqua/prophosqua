#' Download PTMsigDB signatures
#'
#' Downloads and caches a PTMsigDB GMT file from the Broad Institute ssGSEA2.0
#' repository. PTMsigDB v2.0.0 keys its sites on flanking sequences (15 residues
#' centred on the site), so the signatures hold across vertebrates.
#'
#' @param species "mouse" or "human".
#' @param version PTMsigDB version.
#' @param output_dir Directory for the cached file.
#' @param force_download Re-download even when cached.
#' @return Path to the cached GMT file.
#' @references
#' Krug et al. (2019) A Curated Resource for Phosphosite-specific Signature Analysis.
#' Mol Cell Proteomics. doi:10.1074/mcp.TIR118.000943
#' @keywords internal
download_ptmsigdb <- function(species = "mouse", version = "v2.0.0", output_dir = ".", force_download = FALSE) {
  species <- match.arg(species, c("mouse", "human"))
  gmt_url <- paste0(
    "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/",
    version,
    "/all/ptm.sig.db.all.flanking.",
    species,
    ".",
    version,
    ".gmt"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  gmt_cache <- file.path(output_dir, paste0("ptmsigdb_flanking_", species, "_", version, ".gmt"))
  if (!file.exists(gmt_cache) || force_download) {
    message("Downloading PTMsigDB ", version, " signatures for ", species, "...")
    utils::download.file(gmt_url, gmt_cache, quiet = TRUE, mode = "w")
  }
  gmt_cache
}

# Trims a flanking sequence to trim_to residues around its centre; 15 is the
# untrimmed PTMsigDB width.
trim_flanking_seq <- function(seq, trim_to = 15L) {
  if (trim_to == 15L || is.na(seq) || nchar(seq) < trim_to) {
    return(seq)
  }
  trim_each <- (nchar(seq) - trim_to) / 2
  substr(seq, floor(trim_each) + 1, nchar(seq) - ceiling(trim_each))
}

# One ranked list per contrast, named by the upper-case sequence window and
# sorted for GSEA. A window ranked twice keeps its first statistic.
.rank_windows <- function(data, stat_column, trim_to = 15L) {
  .require_columns(data, c("SequenceWindow", stat_column, "contrast"), "ranked sites")
  trim_to <- as.integer(trim_to)
  contrasts <- unique(data$contrast)
  ranks <- lapply(contrasts, function(contrast) {
    rows <- data[data$contrast == contrast, , drop = FALSE]
    windows <- vapply(
      toupper(trimws(as.character(rows$SequenceWindow))),
      trim_flanking_seq,
      character(1),
      trim_to = trim_to,
      USE.NAMES = FALSE
    )
    ranks <- stats::setNames(as.numeric(rows[[stat_column]]), windows)
    ranks <- ranks[!is.na(ranks) & !is.na(names(ranks)) & names(ranks) != ""]
    sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
  })
  stats::setNames(ranks, contrasts)
}

# PTMsigDB names a site by its flanking sequence with a "-p" suffix and marks
# the direction of the regulation with ";u" or ";d". PTM-SEA ranks sequence
# windows and ignores the direction, so both suffixes are dropped.
.ptmsigdb_windows <- function(pathways) {
  lapply(pathways, function(sites) unique(sub("-p(;[ud])?$", "", sites)))
}

# Pre-ranked GSEA of each contrast against PTMsigDB; a contrast sharing fewer
# than ten sites with the database is skipped with a warning.
run_ptmsea <- function(ranks_list, pathways, min_size = 3, max_size = 500, n_perm = 1000, pvalueCutoff = 0.1) {
  pathways <- .ptmsigdb_windows(pathways)
  term2gene <- ptmsigdb_to_term2gene(pathways)
  database_sites <- unique(unlist(pathways))
  results <- lapply(names(ranks_list), function(contrast_name) {
    ranks <- ranks_list[[contrast_name]]
    overlap <- length(intersect(names(ranks), database_sites))
    if (overlap < 10) {
      warning("Contrast '", contrast_name, "': Low overlap with PTMsigDB (", overlap, " sites). Skipping.")
      return(NULL)
    }
    message("Running PTM-SEA for '", contrast_name, "' (", length(ranks), " sites, ", overlap, " overlap)")
    clusterProfiler::GSEA(
      geneList = sort(ranks, decreasing = TRUE),
      TERM2GENE = term2gene,
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = "BH",
      verbose = FALSE,
      by = "fgsea",
      nPermSimple = n_perm
    )
  })
  names(results) <- names(ranks_list)
  Filter(Negate(is.null), results)
}

# Trims every PTMsigDB site to trim_to residues, keeping its "-p" and direction
# suffixes, so the database matches windows trimmed the same way.
trim_ptmsigdb_pathways <- function(pathways, trim_to = 15L) {
  trim_to <- as.integer(trim_to)
  if (trim_to == 15L) {
    return(pathways)
  }
  lapply(pathways, function(sites) {
    suffix <- ifelse(grepl(";[ud]$", sites), sub(".*(-p;[ud])$", "\\1", sites), "-p")
    windows <- vapply(
      sub("-p(;[ud])?$", "", sites),
      trim_flanking_seq,
      character(1),
      trim_to = trim_to,
      USE.NAMES = FALSE
    )
    unique(paste0(windows, suffix))
  })
}

ptmsigdb_to_term2gene <- function(pathways) {
  data.frame(
    term = rep(names(pathways), lengths(pathways)),
    gene = unlist(pathways, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}
