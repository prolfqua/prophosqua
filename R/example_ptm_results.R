# A deterministic final PTM MuData artifact for package vignettes.

# The richer paired DEA input the statistics vignette reads: three groups, `b`
# the control.
.example_ptm_results_dea_pair <- function(root) {
  samples <- expand.grid(replicate = seq_len(3L), G_ = c("a", "b", "c"), stringsAsFactors = FALSE)
  samples$Name <- paste0(samples$G_, samples$replicate)

  proteins <- sprintf("P%03d", seq_len(24L))
  site_index <- seq_len(length(proteins) * 3L)
  site_map <- data.frame(
    protein_Id = rep(proteins, each = 3L),
    posInProtein = rep(c(10L, 20L, 30L), length(proteins)),
    modAA = c("S", "T", "Y", "S")[(site_index - 1L) %% 4L + 1L],
    stringsAsFactors = FALSE
  )
  site_map$site <- paste0(
    site_map$protein_Id,
    "~",
    site_map$modAA,
    site_map$posInProtein
  )

  protein_a <- rep(c(-0.7, -0.3, 0.3, 0.7, 0.1, -0.1), length.out = length(proteins))
  protein_c <- rep(c(0.6, -0.6, 0.2, -0.2, 0.4, -0.4), length.out = length(proteins))
  usage_a <- rep(c(1.4, -1.3, 0.15, -0.15), length.out = nrow(site_map))
  usage_c <- rep(c(-1.25, 0.15, 1.35, -0.1), length.out = nrow(site_map))
  protein_index <- match(site_map$protein_Id, proteins)
  site_a <- protein_a[protein_index] + usage_a
  site_c <- protein_c[protein_index] + usage_c

  background_residues <- strsplit(
    "AAACDEFGGHIIKKLLMNPPQRRSSVVWY",
    "",
    fixed = TRUE
  )[[1]]
  site_map$SequenceWindow <- vapply(
    seq_len(nrow(site_map)),
    function(index) {
      background_group <- (index - 1L) %/% 4L + 1L
      sequence <- vapply(
        seq_len(15L),
        function(position) {
          hash <- digest::digest2int(
            paste(background_group, position, sep = ":")
          )
          background_residues[
            as.double(hash) %% length(background_residues) + 1L
          ]
        },
        character(1)
      )
      sequence[[8L]] <- site_map$modAA[[index]]
      if (usage_a[[index]] > 0.5) {
        sequence[[9L]] <- "P"
      } else if (usage_a[[index]] < -0.5) {
        sequence[c(5L, 6L)] <- "R"
      }
      if (usage_c[[index]] > 0.5) {
        sequence[c(3L, 4L)] <- c("D", "E")
      } else if (usage_c[[index]] < -0.5) {
        sequence[c(11L, 12L)] <- c("L", "V")
      }
      paste0(sequence, collapse = "")
    },
    character(1)
  )
  site_map$description <- paste("protein", site_map$protein_Id)
  site_map$gene_name <- sub("^P", "GENE", site_map$protein_Id)
  site_map$protein_length <- 300L + protein_index * 7L

  group_effect <- function(group, a, c) {
    ifelse(group == "a", a, ifelse(group == "c", c, 0))
  }
  replicate_noise <- c(-0.18, 0.04, 0.14)

  long_protein <- expand.grid(
    Name = samples$Name,
    protein_Id = proteins,
    stringsAsFactors = FALSE
  )
  protein_sample <- match(long_protein$Name, samples$Name)
  protein_feature <- match(long_protein$protein_Id, proteins)
  long_protein$G_ <- samples$G_[protein_sample]
  long_protein$normalized_abundance <-
    .example_base_abundance(long_protein$protein_Id, proteins) +
    group_effect(
      long_protein$G_,
      protein_a[protein_feature],
      protein_c[protein_feature]
    ) +
    replicate_noise[samples$replicate[protein_sample]] *
      (1 + (protein_feature %% 5L) / 20)

  long_site <- expand.grid(
    Name = samples$Name,
    site = site_map$site,
    stringsAsFactors = FALSE
  )
  site_sample <- match(long_site$Name, samples$Name)
  site_feature <- match(long_site$site, site_map$site)
  long_site$protein_Id <- site_map$protein_Id[site_feature]
  long_site$G_ <- samples$G_[site_sample]
  long_site$normalized_abundance <-
    .example_base_abundance(long_site$site, site_map$site) +
    group_effect(
      long_site$G_,
      site_a[site_feature],
      site_c[site_feature]
    ) +
    c(-0.23, 0.06, 0.17)[samples$replicate[site_sample]] *
      (1 + (site_feature %% 7L) / 25)
  faint <- site_feature <= 12L & samples$replicate[site_sample] == 3L
  long_site$normalized_abundance[faint] <- NA_real_

  protein_annotation <- data.frame(
    protein_Id = proteins,
    description = paste("protein", proteins),
    gene_name = sub("^P", "GENE", proteins),
    protein_length = 300L + seq_along(proteins) * 7L,
    stringsAsFactors = FALSE
  )
  contrasts <- c(a_vs_b = "G_a - G_b", c_vs_b = "G_c - G_b")
  paths <- list(phospho = file.path(root, "DEA_phospho"), protein = file.path(root, "DEA_protein"))
  .example_dea_dir(paths$phospho, long_site, list(protein_Id = "protein_Id", site = "site"), site_map, contrasts)
  .example_dea_dir(paths$protein, long_protein, list(protein_Id = "protein_Id"), protein_annotation, contrasts)
  paths
}

.example_enrichment_rank_tables <- function(data) {
  ranked <- data |>
    dplyr::group_by(.data$contrast, .data$SequenceWindow) |>
    dplyr::summarize(
      statistic.site = mean(.data$statistic.site),
      .groups = "drop"
    )
  lapply(split(ranked, ranked$contrast), function(table) {
    table |>
      dplyr::select("SequenceWindow", "statistic.site") |>
      dplyr::arrange(dplyr::desc(.data$statistic.site), .data$SequenceWindow)
  })
}

.example_enrichment_sets <- function(sequences, method) {
  sets <- list(
    proline_directed = sequences[substr(sequences, 9L, 9L) == "P"],
    basophilic = sequences[substr(sequences, 5L, 6L) == "RR"],
    acidophilic = sequences[substr(sequences, 3L, 4L) == "DE"],
    hydrophobic = sequences[substr(sequences, 11L, 12L) == "LV"],
    serine = sequences[substr(sequences, 8L, 8L) == "S"],
    tyrosine = sequences[substr(sequences, 8L, 8L) == "Y"]
  )
  labels <- switch(
    method,
    PTMSEA = c(
      "KINASE-PSP_CDK2",
      "KINASE-PSP_PKACA",
      "KINASE-PSP_CSNK2A1",
      "PERT-PSP_GROWTH_FACTOR",
      "PATHWAY-PSP_SERINE_SIGNALING",
      "DISEASE-PSP_TYROSINE_SIGNALING"
    ),
    KinaseGSEA = c("CDK2", "PKACA", "CSNK2A1", "MAPK1", "PKC", "SRC")
  )
  stats::setNames(sets, labels)
}

.example_term2gene <- function(sets) {
  data.frame(
    term = rep(names(sets), lengths(sets)),
    gene = unlist(sets, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

.example_gsea_results <- function(rank_tables, term2gene, seed) {
  results <- vector("list", length(rank_tables))
  names(results) <- names(rank_tables)
  for (i in seq_along(rank_tables)) {
    table <- rank_tables[[i]]
    ranks <- stats::setNames(table$statistic.site, table$SequenceWindow)
    set.seed(seed + i)
    results[[i]] <- suppressWarnings(clusterProfiler::GSEA(
      sort(ranks, decreasing = TRUE),
      TERM2GENE = term2gene,
      minGSSize = 5L,
      maxGSSize = 100L,
      pvalueCutoff = 1,
      verbose = FALSE,
      seed = TRUE
    ))
  }
  results
}

.example_gsea_table <- function(results, item = "ID") {
  table <- purrr::imap_dfr(results, function(result, contrast) {
    as.data.frame(result) |>
      dplyr::select(
        "ID",
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust",
        "rank",
        "core_enrichment"
      ) |>
      dplyr::mutate(contrast = contrast, .before = 1L)
  })
  if (!identical(item, "ID")) {
    names(table)[names(table) == "ID"] <- item
  }
  table
}

.example_mea_table <- function(results) {
  purrr::imap_dfr(results, function(result, contrast) {
    table <- as.data.frame(result)
    leading <- strsplit(table$core_enrichment, "/", fixed = TRUE)
    data.frame(
      contrast = contrast,
      kinase = table$ID,
      ES = table$enrichmentScore,
      NES = table$NES,
      pvalue = table$pvalue,
      FDR = table$p.adjust,
      n_leading = lengths(leading),
      set_size = table$setSize,
      Leading.substrates = vapply(
        leading,
        paste,
        collapse = ";",
        FUN.VALUE = character(1)
      ),
      stringsAsFactors = FALSE
    )
  })
}

.example_ptm_enrichment_branches <- function(statistics, analysis, table, seed) {
  rank_tables <- .example_enrichment_rank_tables(table)
  sequences <- unique(table$SequenceWindow)

  ptm_results <- .example_gsea_results(
    rank_tables,
    .example_term2gene(.example_enrichment_sets(sequences, "PTMSEA")),
    seed
  )
  ptmsea <- PTMSEA$new(statistics, analysis, list(results = ptm_results, all_clean = .example_gsea_table(ptm_results)))

  kinase_term2gene <- .example_term2gene(.example_enrichment_sets(sequences, "KinaseGSEA"))
  kinase_inputs <- KinaseInputs$new(
    statistics,
    analysis,
    list(seqwindows = data.frame(SequenceWindow = sequences), ranks = rank_tables)
  )
  assignments <- KinaseAssignments$new(kinase_inputs, analysis, list(term2gene = kinase_term2gene))
  kinase_results <- .example_gsea_results(rank_tables, kinase_term2gene, seed + 100L)
  kinase_table <- .example_gsea_table(kinase_results)
  kinase <- KinaseGSEA$new(
    assignments,
    analysis,
    list(
      gsea_results = kinase_results,
      all_results = kinase_table,
      gsea_info = dplyr::count(kinase_table, .data$contrast, name = "terms")
    )
  )

  mea_clean <- .example_mea_table(kinase_results)
  mea_json <- protsea::gsea_result_json_text(protsea::gsea_result_data(
    kinase_results,
    category = "MEA",
    method = "gseapy"
  ))
  motif <- MotifEnrichment$new(assignments, analysis, list(mea_results = mea_clean, gsea_json = mea_json))
  mea_summary <- mea_clean |>
    dplyr::group_by(.data$contrast) |>
    dplyr::summarize(
      total_kinases = dplyr::n(),
      sig_up = sum(.data$FDR < 0.1 & .data$NES > 0),
      sig_down = sum(.data$FDR < 0.1 & .data$NES < 0),
      .groups = "drop"
    )
  mea <- MEA$new(motif, analysis, list(mea_clean = mea_clean, summary_df = mea_summary))
  list(ptmsea, kinase, mea)
}

.example_ptm_enrichments <- function(statistics) {
  tables <- statistics$get_tables()
  analyses <- c("DPA", "DPU", "CF")
  unlist(
    lapply(seq_along(analyses), function(i) {
      .example_ptm_enrichment_branches(statistics, analyses[[i]], tables[[analyses[[i]]]], seed = 4100L + i * 1000L)
    }),
    recursive = FALSE
  )
}

#' Build the Final MuData Used by Package Vignettes
#'
#' The artifact follows the same H5AD import and H5MU storage path as a
#' pipeline run; its enrichments are deterministic GSEA runs on motif sets of
#' the example sequence windows, so no kinase library or PTMsigDB is needed.
#'
#' @param path Destination H5MU path.
#' @return Normalized `path`, invisibly.
#' @keywords internal
example_ptm_results_h5mu <- function(path = tempfile(fileext = ".h5mu")) {
  output_dir <- tempfile("prophosqua_ptm_results_example_")
  dir.create(output_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  dirs <- .example_ptm_results_dea_pair(file.path(output_dir, "dea"))
  input_h5mu <- file.path(output_dir, "inputs.h5mu")
  # The enrichments are built here, not computed, so no reference data is imported.
  parameters <- list(
    run_kinase = TRUE,
    analyses = list(
      dpa = list(stat_column = "statistic.site"),
      dpu = list(stat_column = "statistic.site"),
      cf = list(stat_column = "statistic.site")
    )
  )
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(get_dea_file(dirs$phospho, "AnnData.h5ad")),
    anndataR::read_h5ad(get_dea_file(dirs$protein, "AnnData.h5ad")),
    parameters = parameters
  )
  inputs$write_h5mu(input_h5mu)
  statistics <- suppressWarnings(compute_ptm_results_h5mu(input_h5mu, file.path(output_dir, "statistics.h5mu")))
  result <- PTM_results$new(statistics, enrichments = .example_ptm_enrichments(statistics))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  result$write_h5mu(path)
  invisible(normalizePath(path, mustWork = TRUE))
}
