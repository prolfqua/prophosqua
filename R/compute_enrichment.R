# The enrichment analyses of one standardized PTM result table. Each keeps
# every tested term; the reports apply their own FDR threshold. A result is its
# gseaResult objects, one per contrast, and the tables derived from them, so a
# result read back from its protsea document is the one that was computed.

.compute_ptmsea_tables <- function(data, pathways, stat_column, trim_to, min_size, max_size, n_perm) {
  .ptmsea_result(run_ptmsea(
    .rank_windows(data, stat_column, trim_to),
    pathways,
    min_size,
    max_size,
    n_perm,
    pvalueCutoff = 1
  ))
}

.ptmsea_result <- function(results) {
  all_clean <- extract_gsea_results(results) |>
    dplyr::mutate(
      pathway = .data$ID,
      pathway_short = substr(gsub("^(KINASE|PERT|PATH|DISEASE)-PSP_", "", .data$ID), 1, 40)
    )
  list(results = results, all_clean = all_clean)
}

.compute_kinase_tables <- function(data, term2gene, stat_column, min_size, max_size, n_perm) {
  # Kinase Library writes the phosphorylated residue in lower case, e.g.
  # "PETITIRsGPPSPLP"; the ranked windows are upper case.
  term2gene <- data.frame(term = term2gene$term, gene = toupper(term2gene$gene))
  .kinasegsea_result(lapply(.rank_windows(data, stat_column), function(ranks) {
    clusterProfiler::GSEA(
      geneList = ranks,
      TERM2GENE = term2gene,
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = 1,
      nPermSimple = n_perm,
      verbose = FALSE
    )
  }))
}

.kinasegsea_result <- function(gsea_results) {
  all_results <- dplyr::bind_rows(lapply(names(gsea_results), function(contrast) {
    dplyr::as_tibble(gsea_results[[contrast]]@result) |>
      dplyr::transmute(
        contrast = contrast,
        kinase = .data$ID,
        .data$NES,
        .data$pvalue,
        FDR = .data$p.adjust,
        .data$setSize
      )
  }))
  gsea_info <- dplyr::tibble(
    Contrast = names(gsea_results),
    `Significant Kinases (FDR < 0.25)` = vapply(
      gsea_results,
      function(result) sum(result@result$p.adjust < 0.25, na.rm = TRUE),
      integer(1)
    )
  )
  list(gsea_results = gsea_results, all_results = all_results, gsea_info = gsea_info)
}

# The kinase-library tool computes the MEA and writes it as a protsea document;
# its core enrichment holds the leading substrates.
.mea_result <- function(results) {
  mea_table <- dplyr::bind_rows(lapply(names(results), function(contrast) {
    table <- results[[contrast]]@result
    leading <- strsplit(as.character(table$core_enrichment), "/", fixed = TRUE)
    dplyr::tibble(
      contrast = contrast,
      kinase = as.character(table$ID),
      ES = as.numeric(table$enrichmentScore),
      NES = as.numeric(table$NES),
      pvalue = as.numeric(table$pvalue),
      FDR = as.numeric(table$p.adjust),
      n_leading = lengths(leading),
      set_size = as.integer(table$setSize),
      Leading.substrates = vapply(leading, paste, character(1), collapse = ";")
    )
  }))
  mea_clean <- prepare_enrichment_data(mea_table, "FDR", 0.1)
  summary_df <- mea_clean |>
    dplyr::group_by(.data$contrast) |>
    dplyr::summarize(
      total_kinases = dplyr::n(),
      sig_up = sum(.data$FDR < 0.1 & .data$NES > 0, na.rm = TRUE),
      sig_down = sum(.data$FDR < 0.1 & .data$NES < 0, na.rm = TRUE),
      .groups = "drop"
    )
  list(results = results, mea_clean = mea_clean, summary_df = summary_df)
}

# A motif scan needs a full, uninterrupted window: windows padded with
# underscores at a protein terminus, or shorter than seven residues, carry too
# little context and are dropped rather than scanned.
filter_sequence_windows <- function(data) {
  data |>
    dplyr::filter(
      !is.na(.data$SequenceWindow),
      .data$SequenceWindow != "",
      !grepl("^_", .data$SequenceWindow),
      !grepl("_$", .data$SequenceWindow),
      nchar(.data$SequenceWindow) >= 7
    ) |>
    dplyr::mutate(SequenceWindow = toupper(.data$SequenceWindow))
}

# Motif enrichment walks one ranked list, so a window may appear only once.
# Where sites share a window the most extreme statistic is kept, the one that
# moves the enrichment score most, rather than an average in which two opposing
# sites would cancel.
rank_sites_for_mea <- function(data, stat_column, contrast_name) {
  data |>
    dplyr::filter(.data$contrast == contrast_name, !is.na(.data[[stat_column]])) |>
    dplyr::mutate(statistic.site = .data[[stat_column]]) |>
    dplyr::select("SequenceWindow", "statistic.site") |>
    dplyr::group_by(.data$SequenceWindow) |>
    dplyr::slice(which.max(abs(.data$statistic.site))) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$statistic.site))
}
