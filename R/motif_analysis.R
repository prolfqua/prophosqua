#' Get default list of kinase motifs
#'
#' Returns a data.frame of common kinase motifs defined by regular expressions.
#' Patterns are designed to match 15-mer sequences centered at the phosphosite.
#'
#' @note
#' This list is a compiled "starter set" of well-known consensus motifs from literature
#' (e.g. classical definitions for PKA, PKC, CK2, CDK, etc.). It does NOT contain
#' the entire proprietary PhosphoSitePlus database. Users with access to curated
#' motif databases should load them using `read_kinase_motifs_csv` or pass them
#' as `extra_motifs`.
#'
#' @param extra_motifs Optional data.frame to append to the default list.
#'                     Must have columns: kinase_group, motif_name, pattern, description.
#'
#' @return data.frame with columns: kinase_group, motif_name, pattern, description
#' @export
#' @examples
#' motifs <- get_kinase_motifs()
#' head(motifs)
get_kinase_motifs <- function(extra_motifs = NULL) {
  motifs <- dplyr::bind_rows(
    # --- Basophilic (Arg/Lys rich) ---
    data.frame(
      kinase_group = "Basophilic",
      motif_name = "PKA_AKT_consensus",
      pattern = "R.R..[ST].....", # Core: R.R..[S] (L=5, R=0). Pad R 5.
      description = "R-x-R-x-x-S/T (PKA, AKT, etc. - strong basophilic)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Basophilic",
      motif_name = "PKA_Classic",
      pattern = "R.[ST]..", # Core: R.[S] (L=2, R=0). Pad R 2.
      description = "R-x-S/T (PKA classic)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Basophilic",
      motif_name = "PKC_General",
      pattern = "..[ST].[RK]", # Core: [S].[RK] (L=0, R=2). Pad L 2.
      description = "S/T-x-R/K (PKC general)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Basophilic",
      motif_name = "Akt_RRxS",
      pattern = "RR.[ST]...", # Core: RR.[S] (L=3, R=0). Pad R 3.
      description = "R-R-x-S/T (Akt preference)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Basophilic",
      motif_name = "CAMK2",
      pattern = "R..[ST]...", # Core: R..[S] (L=3, R=0). Pad R 3.
      description = "R-x-x-S/T (CaMKII)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Basophilic",
      motif_name = "Substrate_RK_rich",
      pattern = "[RK][RK].[ST]...", # Core: [RK][RK].[S] (L=3, R=0). Pad R 3.
      description = "Basic residues upstream (R/K-R/K-x-S/T)",
      stringsAsFactors = FALSE
    ),

    # --- Proline-directed (MAPK, CDK, GSK3) ---
    data.frame(
      kinase_group = "Proline-directed",
      motif_name = "Proline_General",
      pattern = ".[ST]P", # Core: [S]P (L=0, R=1). Pad L 1.
      description = "S/T-P (General Proline-directed)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Proline-directed",
      motif_name = "CDK_Consensus",
      pattern = "...[ST]P.[KR]", # Core: [S]P.[KR] (L=0, R=3). Pad L 3.
      description = "S/T-P-x-K/R (CDK1/2/5)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Proline-directed",
      motif_name = "MAPK_Pxh",
      pattern = "P.[ST]P.", # Core: P.[S]P (L=2, R=1). Pad R 1 (makes 2).
      description = "P-x-S/T-P (MAPK)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Proline-directed",
      motif_name = "GSK3_Primed",
      pattern = "....[ST]...[ST]", # Core: [S]...[S] (Target Left). L=0, R=4. Pad L 4.
      description = "S/T-x-x-x-S/T (GSK3 - requires priming pSer at +4)",
      stringsAsFactors = FALSE
    ),

    # --- Acidic (CK2, CK1) ---
    data.frame(
      kinase_group = "Acidic",
      motif_name = "CK2",
      pattern = "...[ST]..E", # Core: [S]..E (L=0, R=3). Pad L 3.
      description = "S/T-x-x-E (CK2 classic)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Acidic",
      motif_name = "CK2_Asp",
      pattern = "...[ST]..D", # Same.
      description = "S/T-x-x-D (CK2 alternative)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Acidic",
      motif_name = "CK1",
      pattern = "...[ST]..[ST]", # User request: Left S is center. (L=0, R=3 in core? No, ...[S]..[S] balances R=3 with L=3? No, core [S]..[S] is R=3. Pad L=3.)
      description = "pS/pT-x-x-S/T (CK1 - requires priming)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Acidic",
      motif_name = "Acidic_General",
      pattern = "..[ST][DE][DE]", # Core: [S][D][D] (L=0, R=2). Pad L 2.
      description = "S/T-D/E-D/E (Acidic patch)",
      stringsAsFactors = FALSE
    ),


    # --- DNA Damage / Other ---
    data.frame(
      kinase_group = "DNA_Damage",
      motif_name = "ATM_ATR",
      pattern = ".[ST]Q", # Core: [S]Q (L=0, R=1). Pad L 1.
      description = "S/T-Q (ATM/ATR)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      kinase_group = "Tyrosine",
      motif_name = "Tyrosine_General",
      pattern = "..Y..", # Core: Y.. (L=0?? No, General is usually context-free? Yxx? If Yxx, L=0, R=2. Pad L 2.)
      description = "Y-x-x (General Tyr - nonspecific)",
      stringsAsFactors = FALSE
    )
  )

  if (!is.null(extra_motifs)) {
    # Validate structure
    req_cols <- c("kinase_group", "motif_name", "pattern", "description")
    if (all(req_cols %in% colnames(extra_motifs))) {
      motifs <- dplyr::bind_rows(motifs, extra_motifs)
    } else {
      warning("extra_motifs must have columns: ", paste(req_cols, collapse = ", "), ". Ignoring.")
    }
  }

  return(motifs)
}


#' Scan sequences for kinase motifs
#'
#' Scans a vector of sequences against a database of regex motifs.
#' Generates a TERM2GENE data.frame for GSEA/Enrichment analysis.
#'
#' @param sequences Character vector of sequences (e.g., 15-mers).
#' @param motif_db data.frame of motifs (must have 'motif_name' and 'pattern' columns).
#'                 Defaults to `get_kinase_motifs()`.
#' @param center_pos Integer. Index of the central phosphosite in the sequence.
#'                   If NULL, matches anywhere. If provided, the regex match MUST
#'                   include/overlap this position. Default is 8 (for 15-mers).
#'
#' @return data.frame with columns 'term' (motif name) and 'gene' (sequence ID).
#'         Sequence ID is the original sequence with "-p" appended (std PTMSEA format).
#' @export
#' @examples
#' # Create some mock sequences (15-mers)
#' # PKA_AKT_consensus pattern: "R.R..[ST]....."
#' # To center S at 8, the pattern R.R..S must match indices 3-13.
#' # We'll create sequences with R at 3 and 5, and S at 8.
#'
#' seqs <- c(
#'   "GGRRRGSEVVVAAAA", # Matches PKA (R at 3,5, S at 8)
#'   "LARRRASVAQLTTAA", # Matches PKA
#'   "AARARAASVAAAAAA", # Matches PKA
#'   "RRRRAASVAAAAAAA", # Matches PKA and others
#'   "GGGGGGSGGGGGGGG" # Non-match control
#' )
#'
#' motif_db <- get_kinase_motifs()
#' res <- scan_motifs(seqs, motif_db)
#' head(res)
scan_motifs <- function(sequences, motif_db = get_kinase_motifs(), center_pos = 8L) {
  sequences <- unique(toupper(trimws(as.character(sequences))))

  results_list <- list()

  for (i in seq_len(nrow(motif_db))) {
    motif_name <- motif_db$motif_name[i]
    pattern <- motif_db$pattern[i]

    # simple grep first to find candidates
    matches <- grep(pattern, sequences, value = TRUE)

    if (length(matches) > 0) {
      if (!is.null(center_pos)) {
        # Check if the match is centered at center_pos
        # The patterns are now defined to be odd-length and centered on the phosphosite.
        # So we can calculate the geometric center of the match in the sequence.

        locs <- regexpr(pattern, matches)
        starts <- as.integer(locs)
        lengths <- attr(locs, "match.length")

        # Calculate geometric center index of the match within the sequence
        # (lengths - 1) / 2 is the offset from start to the center base
        match_centers <- starts + (lengths - 1) / 2

        # Valid if the match center aligns with the sequence center_pos
        valid_idx <- (match_centers == center_pos)

        matches <- matches[valid_idx]
      }

      if (length(matches) > 0) {
        site_ids <- paste0(matches, "-p")
        results_list[[motif_name]] <- data.frame(
          term = motif_name,
          gene = site_ids,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results_list) == 0) {
    return(data.frame(term = character(), gene = character()))
  }

  term2gene <- dplyr::bind_rows(results_list)
  return(term2gene)
}


#' Run Motif-based GSEA
#'
#' @param ranks_list Named list of named numeric vectors (stats per contrast).
#' @param motif_term2gene data.frame from `scan_motifs`.
#' @param min_size Minimum set size.
#' @param max_size Maximum set size.
#' @param n_perm Number of permutations.
#' @param pvalueCutoff P-value cutoff.
#'
#' @return List of GSEA results.
#' @export
#' @examples
#' # Mock data: 100 sequences, first 10 weighted to include PKA motifs
#' seqs <- paste0("SEQ", 1:100, "-p")
#' ranks <- rnorm(100)
#' # Make top ranked sequences PKA-like (MotifA)
#' ranks[1:5] <- c(3.0, 2.5, 2.0, 1.9, 1.8)
#' names(ranks) <- seqs
#' ranks <- sort(ranks, decreasing = TRUE)
#' rankList <- list(contrast1 = ranks)
#'
#' # Define MotifA as PKA-like sequences (matches SEQ1-SEQ5)
#' motif_term2gene <- data.frame(
#'   term = rep("MotifA_PKA", 5),
#'   gene = seqs[1:5],
#'   stringsAsFactors = FALSE
#' )
#'
#' # Run GSEA
#' res <- run_motif_gsea(rankList, motif_term2gene, min_size = 3, pvalueCutoff = 1.0)
run_motif_gsea <- function(ranks_list,
                           motif_term2gene,
                           min_size = 3,
                           max_size = 5000,
                           n_perm = 1000,
                           pvalueCutoff = 0.1) {
  if (nrow(motif_term2gene) == 0) {
    warning("motif_term2gene is empty. Returning empty results.")
    return(list())
  }

  results <- list()

  for (contrast_name in names(ranks_list)) {
    ranks <- ranks_list[[contrast_name]]

    # Sort ranks descending (required for GSEA)
    ranks <- sort(ranks, decreasing = TRUE)

    # Check overlap
    overlap_genes <- intersect(names(ranks), motif_term2gene$gene)
    if (length(overlap_genes) < min_size) {
      warning(paste("Contrast", contrast_name, ": insufficient overlap with motif genes. Skipping."))
      next
    }

    # Attempt GSEA
    tryCatch(
      {
        res <- clusterProfiler::GSEA(
          geneList = ranks,
          TERM2GENE = motif_term2gene,
          minGSSize = min_size,
          maxGSSize = max_size,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = "BH",
          verbose = FALSE,
          nPermSimple = n_perm
        )
        results[[contrast_name]] <- res
      },
      error = function(e) {
        warning(paste("GSEA failed for contrast", contrast_name, ":", e$message))
      }
    )
  }

  return(results)
}

#' Read kinase motifs from CSV
#'
#' Helper function to load a custom motif database from a file.
#'
#' @param file Path to CSV file. Expected columns: "kinase_group", "motif_name", "pattern", "description".
#' @return data.frame suitable for `scan_motifs` or `get_kinase_motifs(extra_motifs=...)`
#' @export
read_kinase_motifs_csv <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  df <- utils::read.csv(file, stringsAsFactors = FALSE)

  req_cols <- c("motif_name", "pattern")
  missing <- setdiff(req_cols, colnames(df))
  if (length(missing) > 0) {
    stop("CSV missing required columns: ", paste(missing, collapse = ", "))
  }

  # Fill optional columns if missing
  if (!"kinase_group" %in% colnames(df)) df$kinase_group <- "Custom"
  if (!"description" %in% colnames(df)) df$description <- ""

  return(df)
}

#' Derive motifs from PTMsigDB GMT file
#'
#' Analyzes a GMT file (containing lists of aligned substrate sequences for each kinase)
#' and generates a "consensus" regular expression for each kinase based on residue frequency.
#'
#' @param gmt_file Path to PTMsigDB GMT file (e.g., "ptm.sig.db.all.flanking...gmt").
#' @param freq_threshold Numeric (0-1). Minimum frequency for a residue (or group) to be included in the regex.
#'                       Default 0.3 (30 percent).
#' @param min_seqs Integer. Minimum number of sequences required to derive a motif for a kinase. Default 10.
#' @param trim_to Integer. Width of the analyzed window (must be odd). Default 15.
#'
#' @return data.frame with columns: kinase_group, motif_name, pattern, description
#' @export
derive_motifs_from_gmt <- function(gmt_file, freq_threshold = 0.3, min_seqs = 10, trim_to = 15) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required to read GMT files.")
  }

  # Read GMT: List of KinaseName -> Vector of "Sequence-p" IDs
  pathways <- fgsea::gmtPathways(gmt_file)

  results <- list()

  center <- ceiling(trim_to / 2)

  for (kinase_name in names(pathways)) {
    site_ids <- pathways[[kinase_name]]

    # Clean IDs to get sequences: "SEQUENCE-p" -> "SEQUENCE"
    # Also remove direction suffixes ;u/;d if present
    seqs <- gsub("-p(;[ud])?$", "", site_ids)
    seqs <- unique(seqs)

    # Filter by length (must be >= trim_to) and min_seqs
    valid_seqs <- seqs[nchar(seqs) >= trim_to]
    if (length(valid_seqs) < min_seqs) next

    # Trim to window
    # Assuming sequences are centered.
    # If standard PTMsigDB, they are 15-mers.
    trimmed_seqs <- sapply(valid_seqs, function(s) {
      n <- nchar(s)
      start <- floor((n - trim_to) / 2) + 1
      substr(s, start, start + trim_to - 1)
    })

    # Build PFM (Position Frequency Matrix)
    # Rows: Positions (1..15), Cols: AA
    aa_alphabet <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    n_seqs <- length(trimmed_seqs)

    pattern_parts <- character(trim_to)

    for (pos in 1:trim_to) {
      residues <- substr(trimmed_seqs, pos, pos)
      counts <- table(factor(residues, levels = aa_alphabet))
      freqs <- counts / n_seqs

      # Determine consensus for this position
      # 1. Is there a single dominant AA?
      max_aa <- names(which.max(freqs))
      max_freq <- max(freqs)

      if (pos == center) {
        # Central residue: usually S/T or Y.
        # let's allow data, but usually it's [ST] or [Y].
        if (sum(freqs[c("S", "T")]) > 0.8) {
          pattern_parts[pos] <- "([ST])"
        } else if (freqs["Y"] > 0.8) {
          pattern_parts[pos] <- "(Y)"
        } else {
          # Fallback
          pattern_parts[pos] <- "([STY])"
        }
        next
      }

      if (max_freq >= freq_threshold) {
        pattern_parts[pos] <- max_aa
      } else {
        # 2. Check for functional groups if single failed
        # Basic: R, K
        if (sum(freqs[c("R", "K")]) >= freq_threshold) {
          pattern_parts[pos] <- "[RK]"
        }
        # Acidic: D, E
        else if (sum(freqs[c("D", "E")]) >= freq_threshold) {
          pattern_parts[pos] <- "[DE]"
        }
        # Hydrophobic (rough set): L, V, I, F, M
        else if (sum(freqs[c("L", "V", "I", "F", "M")]) >= freq_threshold) {
          pattern_parts[pos] <- "[LVIFM]"
        } else {
          pattern_parts[pos] <- "."
        }
      }
    }

    # Construct Regex
    regex <- paste0(pattern_parts, collapse = "")

    # Only keep if it's not effectively empty
    # Count non-dots (excluding center)
    complexity <- sum(pattern_parts[-center] != ".")

    if (complexity >= 1) {
      # Clean naming: strip "KINASE_TARGET" prefixes often found in PTMsigDB
      clean_name <- gsub("_TARGET$", "", kinase_name)

      results[[clean_name]] <- data.frame(
        kinase_group = "Derived",
        motif_name = clean_name,
        pattern = regex,
        description = paste0("Derived from ", n_seqs, " sites (", round(max_freq * 100), "% consensus)"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    return(dplyr::bind_rows(results))
  } else {
    return(data.frame())
  }
}
