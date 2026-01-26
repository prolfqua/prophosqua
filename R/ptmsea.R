# PTM-SEA (Post-Translational Modification Signature Enrichment Analysis) Functions
#
# Minimal functions for running PTM-SEA analysis using PTMsigDB signatures and fgsea.


#' Download PTMsigDB signatures
#'
#' Downloads and caches PTMsigDB GMT files from the Broad Institute ssGSEA2.0 repository.
#' PTMsigDB v2.0.0 uses flanking sequences (15 amino acids centered on the phosphosite),
#' making the analysis species-invariant within vertebrates.
#'
#' @param species Character. Species for signatures ("mouse" or "human"). Default "mouse".
#' @param version Character. PTMsigDB version. Default "v2.0.0".
#' @param output_dir Character. Directory for cached files. Default current directory.
#' @param force_download Logical. Force re-download even if cached. Default FALSE.
#'
#' @return Character path to the cached GMT file.
#' @export
#'
#' @references
#' Krug et al. (2019) A Curated Resource for Phosphosite-specific Signature Analysis.
#' Mol Cell Proteomics. doi:10.1074/mcp.TIR118.000943
#'
#' @examples
#' \dontrun{
#' gmt_path <- download_ptmsigdb(species = "mouse", output_dir = "ptmsea_cache")
#' pathways <- fgsea::gmtPathways(gmt_path)
#' }
download_ptmsigdb <- function(species = "mouse",
                              version = "v2.0.0",
                              output_dir = ".",
                              force_download = FALSE) {
  valid_species <- c("mouse", "human")
  if (!species %in% valid_species) {
    stop("Species must be one of: ", paste(valid_species, collapse = ", "))
  }

  base_url <- "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb"
  gmt_url <- paste0(base_url, "/", version, "/all/ptm.sig.db.all.flanking.", species, ".", version, ".gmt")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  gmt_cache <- file.path(output_dir, paste0("ptmsigdb_flanking_", species, "_", version, ".gmt"))

  if (!file.exists(gmt_cache) || force_download) {
    message("Downloading PTMsigDB ", version, " signatures for ", species, "...")
    tryCatch(
      {
        utils::download.file(gmt_url, gmt_cache, quiet = TRUE, mode = "w")
        message("Downloaded: ", gmt_cache)
      },
      error = function(e) {
        stop("Failed to download PTMsigDB: ", e$message, "\nURL: ", gmt_url)
      }
    )
  } else {
    message("Using cached PTMsigDB: ", gmt_cache)
  }

  return(gmt_cache)
}


#' Trim flanking sequence to specified width
#'
#' Trims a flanking sequence centered on the phosphosite (position 8 in 15-mer).
#' Removes equal amino acids from N and C termini.
#'
#' @param seq Character. Flanking sequence (15-mer expected).
#' @param trim_to Integer. Target width: 15 (no trim), 13, or 11. Default 15.
#'
#' @return Character. Trimmed sequence.
#' @keywords internal
trim_flanking_seq <- function(seq, trim_to = 15L) {
  if (trim_to == 15L) {
    return(seq)
  }
  if (is.na(seq) || seq == "") {
    return(seq)
  }
  n <- nchar(seq)
  if (n < trim_to) {
    return(seq)
  }
  trim_each <- (n - trim_to) / 2
  start <- floor(trim_each) + 1
  end <- n - ceiling(trim_each)
  substr(seq, start, end)
}


#' Prepare data for PTM-SEA analysis
#'
#' Converts phosphosite data to a list of named vectors suitable for GSEA.
#' Each contrast produces one named vector where names are site IDs in PTMsigDB
#' format (SEQUENCE-p) and values are t-statistics (or other ranking metric).
#'
#' @param data data.frame containing phosphosite data with statistics per contrast
#' @param stat_column Character. Column name for ranking statistic. Use "statistic.site"
#'   for DPA or "tstatistic_I" for DPU.
#' @param seq_window_col Character. Column name for sequence windows. Default "SequenceWindow".
#' @param contrast_col Character. Column name for contrasts. Default "contrast".
#' @param trim_to Character. Trim flanking sequences to this width. Options: "11" (default),
#'   "13", or "15" (no trimming). Trimming increases overlap with PTMsigDB.
#'
#' @return List with two elements:
#'   \itemize{
#'     \item ranks: Named list of named numeric vectors (one per contrast)
#'     \item dropped: Named list of dropped duplicates (one per contrast with duplicates)
#'   }
#' @export
#'
#' @examples
#' @examples
#' # Prepare mock data
#' data <- data.frame(
#'   contrast = c("A", "A", "B", "B"),
#'   SequenceWindow = c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB", "AAAAAAASAAAAAAA", "CCCCCCCSCCCCCCC"),
#'   statistic.site = c(2.5, 1.5, -2.0, -1.0),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Prepare DPA data using t-statistic for ranking
#' prep_dpa <- ptmsea_data_prep(data, stat_column = "statistic.site")
#' prep_dpa$ranks # ranked lists per contrast
#' prep_dpa$dropped # duplicates that were dropped
#'
#' # Trim to 13-mer to potentially increase PTMsigDB overlap
#' prep_dpa_13 <- ptmsea_data_prep(data, stat_column = "statistic.site", trim_to = "13")
ptmsea_data_prep <- function(data,
                             stat_column,
                             seq_window_col = "SequenceWindow",
                             contrast_col = "contrast",
                             trim_to = c("11", "13", "15")) {
  trim_to <- match.arg(trim_to)
  trim_width <- as.integer(trim_to)

  # Validate input columns
  required_cols <- c(seq_window_col, stat_column, contrast_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  contrasts <- unique(data[[contrast_col]])
  ranks_list <- list()
  dropped_list <- list()

  for (contrast_name in contrasts) {
    contrast_data <- data[data[[contrast_col]] == contrast_name, ]

    # Create site IDs in PTMsigDB format: uppercase sequence + "-p"
    seq_windows <- toupper(trimws(as.character(contrast_data[[seq_window_col]])))

    # Trim sequences if requested
    if (trim_width < 15L) {
      seq_windows <- sapply(seq_windows, trim_flanking_seq, trim_to = trim_width, USE.NAMES = FALSE)
    }

    site_ids <- paste0(seq_windows, "-p")
    statistics <- as.numeric(contrast_data[[stat_column]])

    # Create named vector
    ranks <- statistics
    names(ranks) <- site_ids

    # Remove NA entries
    valid_idx <- !is.na(ranks) & !is.na(names(ranks))
    ranks <- ranks[valid_idx]

    # Handle duplicated flanking sequences - keep first, store dropped
    if (any(duplicated(names(ranks)))) {
      dup_idx <- duplicated(names(ranks))
      dropped_list[[contrast_name]] <- ranks[dup_idx]
      ranks <- ranks[!dup_idx]
    }

    # Sort by statistic (descending) for GSEA
    ranks <- sort(ranks, decreasing = TRUE)

    ranks_list[[contrast_name]] <- ranks
  }

  return(list(ranks = ranks_list, dropped = dropped_list))
}


#' Run PTM-SEA analysis
#'
#' Applies fgsea to each contrast's ranked site list against PTMsigDB pathways.
#'
#' @param ranks_list Named list of named numeric vectors from ptmsea_data_prep().
#' @param pathways Named list of pathway definitions (e.g., from fgsea::gmtPathways()).
#' @param min_size Integer. Minimum pathway size. Default 3.
#' @param max_size Integer. Maximum pathway size. Default 500.
#' @param n_perm Integer. Number of permutations. Default 1000.
#'
#' @return Named list of fgsea result data.frames, one per contrast.
#' @export
#'
#' @examples
#' # Mock ranked lists
#' # Must have > 10 overlap with pathways to run
#' seqs <- paste0("SEQ", 1:15, "-p")
#' rank_vec <- setNames(rnorm(15), seqs)
#' ranks_list <- list(contrast1 = sort(rank_vec, decreasing = TRUE))
#'
#' # Mock pathways
#' # Create pathways that overlap with the sequences
#' pathways <- list(
#'   PathwayA = seqs[1:10],
#'   PathwayB = seqs[6:15]
#' )
#'
#' # Run PTM-SEA (adjusting min_size for small mock data)
#' results <- run_ptmsea(ranks_list, pathways, min_size = 1, pvalueCutoff = 1.0)
# Internal function to split pathways by direction
split_ptmsigdb_pathways <- function(pathways) {
  result <- list()
  for (name in names(pathways)) {
    genes <- pathways[[name]]
    up_genes <- genes[grepl(";u$", genes)]
    down_genes <- genes[grepl(";d$", genes)]
    up_genes <- gsub(";u$", "", up_genes)
    down_genes <- gsub(";d$", "", down_genes)
    if (length(up_genes) > 0) result[[paste0(name, "_up")]] <- unique(up_genes)
    if (length(down_genes) > 0) result[[paste0(name, "_down")]] <- unique(down_genes)
  }
  return(result)
}


#' Run directional PTM-SEA analysis
#'
#' Splits PTMsigDB pathways by direction (up/down) and runs fgsea twice
#' per contrast (ascending and descending rank order).
#'
#' @param ranks_list Named list of named numeric vectors from ptmsea_data_prep().
#' @param pathways Named list of pathway definitions (e.g., from fgsea::gmtPathways()).
#' @param min_size Integer. Minimum pathway size. Default 3.
#' @param max_size Integer. Maximum pathway size. Default 500.
#' @param n_perm Integer. Number of permutations. Default 1000.
#' @param pvalueCutoff Numeric. Adjusted p-value cutoff. Default 0.1.
#'
#' @return Named list of clusterProfiler gseaResult objects, named <contrast>_<ascending/descending>.
#' @export
run_ptmsea_up_down <- function(ranks_list,
                               pathways,
                               min_size = 3,
                               max_size = 500,
                               n_perm = 1000,
                               pvalueCutoff = 0.1) {
  # Split pathways by direction
  pathways_split <- split_ptmsigdb_pathways(pathways)

  # Convert to TERM2GENE format
  term2gene <- ptmsigdb_to_term2gene(pathways_split)

  results <- list()

  for (contrast_name in names(ranks_list)) {
    ranks <- ranks_list[[contrast_name]]

    # Check overlap
    all_pathway_genes <- unique(unlist(pathways_split))
    overlap <- length(intersect(names(ranks), all_pathway_genes))

    if (overlap < 10) {
      warning(
        "Contrast '", contrast_name, "': Low overlap with PTMsigDB (",
        overlap, " sites). Skipping."
      )
      next
    }

    message(
      "Running PTM-SEA for '", contrast_name, "' (",
      length(ranks), " sites, ", overlap, " overlap)"
    )

    # Descending
    ranks_desc <- sort(ranks, decreasing = TRUE)
    res_desc <- clusterProfiler::GSEA(
      geneList = ranks_desc,
      TERM2GENE = term2gene,
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = "BH",
      verbose = FALSE,
      by = "fgsea",
      nPermSimple = n_perm
    )
    results[[paste0(contrast_name, "_descending")]] <- res_desc

    # Ascending (negate values, sort descending - same as ascending on original)
    ranks_asc <- sort(-ranks, decreasing = TRUE)
    res_asc <- clusterProfiler::GSEA(
      geneList = ranks_asc,
      TERM2GENE = term2gene,
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = "BH",
      verbose = FALSE,
      by = "fgsea",
      nPermSimple = n_perm
    )
    results[[paste0(contrast_name, "_ascending")]] <- res_asc
  }

  return(results)
}


#' Run PTM-SEA without direction
#'
#' Runs GSEA on PTMsigDB pathways, ignoring direction suffixes (;u/;d).
#'
#' @param ranks_list Named list of ranked vectors from ptmsea_data_prep().
#' @param pathways List of pathways from fgsea::gmtPathways().
#' @param min_size Integer. Minimum pathway size. Default 3.
#' @param max_size Integer. Maximum pathway size. Default 500.
#' @param n_perm Integer. Number of permutations. Default 1000.
#' @param pvalueCutoff Numeric. Adjusted p-value cutoff. Default 0.1.
#'
#' @return Named list of clusterProfiler gseaResult objects (one per contrast).
#' @export
run_ptmsea <- function(ranks_list,
                       pathways,
                       min_size = 3,
                       max_size = 500,
                       n_perm = 1000,
                       pvalueCutoff = 0.1) {
  # Strip direction suffixes (;u or ;d) from PTMsigDB site IDs
  pathways <- lapply(pathways, function(genes) {
    unique(gsub(";[ud]$", "", genes))
  })

  # Convert to TERM2GENE format for clusterProfiler
  term2gene <- ptmsigdb_to_term2gene(pathways)

  results <- lapply(names(ranks_list), function(contrast_name) {
    ranks <- ranks_list[[contrast_name]]

    # Check overlap
    all_pathway_genes <- unique(unlist(pathways))
    overlap <- length(intersect(names(ranks), all_pathway_genes))

    if (overlap < 10) {
      warning(
        "Contrast '", contrast_name, "': Low overlap with PTMsigDB (",
        overlap, " sites). Skipping."
      )
      return(NULL)
    }

    message(
      "Running PTM-SEA for '", contrast_name, "' (",
      length(ranks), " sites, ", overlap, " overlap)"
    )

    # Ensure sorted descending for GSEA
    ranks <- sort(ranks, decreasing = TRUE)

    res <- clusterProfiler::GSEA(
      geneList = ranks,
      TERM2GENE = term2gene,
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = "BH",
      verbose = FALSE,
      by = "fgsea",
      nPermSimple = n_perm
    )

    return(res)
  })

  names(results) <- names(ranks_list)
  # Remove NULL entries (skipped contrasts)
  results <- results[!sapply(results, is.null)]

  return(results)
}


#' Prepare data for PTM-SEA ORA analysis
#'
#' Prepares enriched and background site lists for over-representation analysis.
#' When fc_col is provided, significant sites are split into up and down sets.
#' Only sites that intersect with PTMsigDB are included.
#'
#' @param data data.frame with flanking sequences, scores, and optionally fold changes
#' @param ptmsigdb_sites Character vector. Site IDs from PTMsigDB (stripped, without ;u/;d).
#'   Use \code{unique(gsub(";[ud]$", "", unlist(pathways)))}.
#' @param score_col Character. Column for significance score (e.g., "FDR_I").
#' @param fc_col Character or NULL. Column for fold change. If provided, returns
#'   separate enriched_up and enriched_down. Default NULL.
#' @param seq_window_col Character. Column for flanking sequences. Default "SequenceWindow".
#' @param threshold Numeric. Significance threshold for score_col. Default 0.05.
#' @param trim_to Character. Trim width: "11" (default), "13", or "15".
#'
#' @return List with enriched (or enriched_up/enriched_down if fc_col) and background.
#' @export
#' @importFrom dplyr mutate filter distinct pull
ptmsea_ora_prep <- function(data,
                            ptmsigdb_sites,
                            score_col,
                            fc_col = NULL,
                            seq_window_col = "SequenceWindow",
                            threshold = 0.05,
                            trim_to = c("11", "13", "15")) {
  trim_to <- match.arg(trim_to)
  trim_width <- as.integer(trim_to)

  # Validate required columns
  required_cols <- c(seq_window_col, score_col, fc_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Prepare site IDs
  data <- data |>
    dplyr::mutate(
      .site_id = paste0(sapply(toupper(trimws(.data[[seq_window_col]])),
        trim_flanking_seq,
        trim_to = trim_width, USE.NAMES = FALSE
      ), "-p")
    ) |>
    dplyr::filter(!is.na(.data[[score_col]]))

  # Background: all sites intersecting with PTMsigDB
  background <- data |>
    dplyr::distinct(.site_id) |>
    dplyr::pull(.site_id) |>
    intersect(ptmsigdb_sites)

  if (!is.null(fc_col)) {
    data <- data |> dplyr::filter(!is.na(.data[[fc_col]]))

    enriched_up <- data |>
      dplyr::filter(.data[[score_col]] < threshold, .data[[fc_col]] > 0) |>
      dplyr::distinct(.site_id) |>
      dplyr::pull(.site_id) |>
      intersect(ptmsigdb_sites)

    enriched_down <- data |>
      dplyr::filter(.data[[score_col]] < threshold, .data[[fc_col]] < 0) |>
      dplyr::distinct(.site_id) |>
      dplyr::pull(.site_id) |>
      intersect(ptmsigdb_sites)

    return(list(enriched_up = enriched_up, enriched_down = enriched_down, background = background))
  } else {
    enriched <- data |>
      dplyr::filter(.data[[score_col]] < threshold) |>
      dplyr::distinct(.site_id) |>
      dplyr::pull(.site_id) |>
      intersect(ptmsigdb_sites)

    return(list(enriched = enriched, background = background))
  }
}


#' Trim PTMsigDB pathways to shorter flanking sequences
#'
#' Trims all site IDs in PTMsigDB pathways to specified width. Use this when
#' trimming your data with ptmsea_data_prep(trim_to=...) to ensure both sides match.
#'
#' @param pathways Named list of pathways from fgsea::gmtPathways().
#' @param trim_to Character. Target width: "11" (default), "13", or "15" (no trim).
#'
#' @return Named list of pathways with trimmed site IDs.
#' @export
#'
#' @examples
#' # Mock pathways
#' pathways <- list(
#'   PathwayA = c("AAAAAAASAAAAAAA-p", "BBBBBBBSBBBBBBB-p"),
#'   PathwayB = c("CCCCCCCSCCCCCCC-p")
#' )
#'
#' # Trim to 11-mer
#' pathways_11 <- trim_ptmsigdb_pathways(pathways, trim_to = "11")
trim_ptmsigdb_pathways <- function(pathways, trim_to = c("11", "13", "15")) {
  trim_to <- match.arg(trim_to)
  trim_width <- as.integer(trim_to)

  if (trim_width == 15L) {
    return(pathways)
  }

  lapply(pathways, function(genes) {
    # Site IDs are like "ABCDEFGHIJKLMNO-p" or "ABCDEFGHIJKLMNO-p;u"
    sapply(genes, function(g) {
      has_dir <- grepl(";[ud]$", g)
      dir_suffix <- ifelse(has_dir, sub(".*(-p;[ud])$", "\\1", g), "-p")
      seq <- sub("-p(;[ud])?$", "", g)
      paste0(trim_flanking_seq(seq, trim_width), dir_suffix)
    }, USE.NAMES = FALSE) |> unique()
  })
}


#' Convert PTMsigDB pathways to TERM2GENE format
#'
#' Converts pathway list from fgsea::gmtPathways() to data.frame format
#' required by clusterProfiler::enricher.
#'
#' @param pathways Named list of pathways (from fgsea::gmtPathways)
#'
#' @return data.frame with columns: term, gene
#' @export
ptmsigdb_to_term2gene <- function(pathways) {
  term2gene <- do.call(rbind, lapply(names(pathways), function(term) {
    data.frame(term = term, gene = pathways[[term]], stringsAsFactors = FALSE)
  }))
  return(term2gene)
}
