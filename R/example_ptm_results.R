# Deterministic final PTM results for package vignettes, laid out as a
# pipeline run lays them out.

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

.example_ptm_enrichment_stages <- function(statistics, analysis, table, seed) {
  rank_tables <- .example_enrichment_rank_tables(table)
  sequences <- unique(table$SequenceWindow)

  ptm_results <- .example_gsea_results(
    rank_tables,
    .example_term2gene(.example_enrichment_sets(sequences, "PTMSEA")),
    seed
  )
  ptmsea <- PTMSEA$new(statistics, analysis, .ptmsea_result(ptm_results))

  kinase_term2gene <- .example_term2gene(.example_enrichment_sets(sequences, "KinaseGSEA"))
  kinase_inputs <- KinaseInputs$new(
    statistics,
    analysis,
    list(seqwindows = data.frame(SequenceWindow = sequences), ranks = rank_tables)
  )
  assignments <- KinaseAssignments$new(kinase_inputs, analysis, list(term2gene = kinase_term2gene))
  kinase_results <- .example_gsea_results(rank_tables, kinase_term2gene, seed + 100L)
  kinase <- KinaseGSEA$new(assignments, analysis, .kinasegsea_result(kinase_results))
  # The kinase-library tool computes the MEA of a pipeline run; the example
  # uses the kinase GSEA of the same sets in its place.
  mea <- MEA$new(assignments, analysis, .mea_result(kinase_results))
  list(ptmsea, kinase_inputs, assignments, kinase, mea)
}

#' Build the Final PTM Results Used by Package Vignettes
#'
#' Written as a pipeline run writes them: the final MuData and, beside it, the
#' enrichment files of each analysis. The data follow the same H5AD import and
#' H5MU storage path as a pipeline run; the enrichments are deterministic GSEA
#' runs on motif sets of the example sequence windows, so no kinase library or
#' PTMsigDB is needed.
#'
#' @param root Destination directory.
#' @return The final MuData file, `PTM_results.h5mu` in `root`, invisibly.
#' @keywords internal
example_ptm_results <- function(root = tempfile("prophosqua_ptm_results_")) {
  work <- tempfile("prophosqua_ptm_results_example_")
  dir.create(work, recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  dirs <- .example_ptm_results_dea_pair(file.path(work, "dea"))
  input_h5mu <- file.path(work, "inputs.h5mu")
  statistics_h5mu <- file.path(work, "statistics.h5mu")
  # The enrichments are built here, not computed, so no reference data is imported.
  parameters <- list(
    run_kinase = TRUE,
    analyses = list(
      dpa = list(subdir = "PTM_DPA", stat_column = "statistic.site"),
      dpu = list(subdir = "PTM_DPU", stat_column = "statistic.site"),
      cf = list(subdir = "PTM_CF_DPU", stat_column = "statistic.site")
    )
  )
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(get_dea_file(dirs$phospho, "AnnData.h5ad")),
    anndataR::read_h5ad(get_dea_file(dirs$protein, "AnnData.h5ad")),
    parameters = parameters
  )
  inputs$write_h5mu(input_h5mu)
  statistics <- suppressWarnings(compute_ptm_results_h5mu(input_h5mu, statistics_h5mu))
  statistics_hash <- .ptm_file_sha256(statistics_h5mu)
  files <- .ptm_enrichment_files(parameters, root)
  tables <- statistics$get_tables()
  analyses <- c("DPA", "DPU", "CF")
  for (i in seq_along(analyses)) {
    stages <- .example_ptm_enrichment_stages(statistics, analyses[[i]], tables[[analyses[[i]]]], 4100L + i * 1000L)
    for (stage in stages) {
      .write_ptm_stage_file(stage, files[[.ptm_key(class(stage)[1L], analyses[[i]])]], statistics_hash)
    }
  }
  path <- file.path(root, "PTM_results.h5mu")
  PTM_results$new(statistics, files, statistics_hash)$write_h5mu(path)
  invisible(normalizePath(path, mustWork = TRUE))
}
