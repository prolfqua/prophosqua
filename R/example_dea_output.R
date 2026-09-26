# A synthetic pair of prolfquapp DEA output directories, made by prolfquapp's
# own DEA on synthetic abundances, so the examples and tests read the artifact
# exactly as a pipeline run does: a `Results_WU_*` directory with `AnnData.h5ad`.

# Two groups of three replicates, `b` the control.
.example_dea_samples <- function() {
  samples <- expand.grid(replicate = 1:3, G_ = c("a", "b"), stringsAsFactors = FALSE)
  samples$Name <- paste0(samples$G_, samples$replicate)
  samples
}

.example_dea_config <- function(hierarchy) {
  config <- prolfqua::AnalysisConfiguration$new()
  config$hierarchy <- hierarchy
  config$hierarchy_depth <- length(hierarchy)
  config$factors <- list(G_ = "G_")
  config$factor_depth <- 1
  config$isotope_label <- "isotopeLabel"
  config$ident_q_value <- "qValue"
  config$nr_children <- "nr_children"
  config$sample_name <- "Name"
  config$file_name <- "Name"
  config$set_response("intensity")
  config
}

#' Write One Example DEA Output Directory
#'
#' Runs prolfquapp's DEA, the `lm_impute` model, on example abundances and
#' writes its artifact, `Results_WU_example/AnnData.h5ad`. The abundances are
#' log2 values; prolfquapp's log2 transform of their exponent gives them back.
#'
#' @param dea_dir Directory to create.
#' @param long Long-format log2 abundances in `normalized_abundance`.
#' @param hierarchy Named list of hierarchy keys.
#' @param row_annot Feature annotation keyed by the hierarchy keys, with a
#'   `description` column.
#' @param contrasts Named contrast expressions.
#' @return `dea_dir`.
#' @keywords internal
.example_dea_dir <- function(dea_dir, long, hierarchy, row_annot, contrasts) {
  long$intensity <- 2^long$normalized_abundance
  long$normalized_abundance <- NULL
  long$isotopeLabel <- "light"
  long$qValue <- 0
  long$nr_children <- 1L
  config <- .example_dea_config(hierarchy)
  lfqdata <- prolfqua::LFQData$new(suppressMessages(prolfqua::setup_analysis(long, config)), config)
  row_annot$nrPeptides <- 1L
  annotation <- prolfquapp::ProteinAnnotation$new(lfqdata, row_annot = row_annot, description = "description")
  dea_config <- prolfquapp::make_DEA_config_R6(PATH = tempfile("prophosqua-example-dea-"), model = "lm_impute")
  dea_config$processing_options$transform <- "log2"
  prep <- prolfquapp::ProteinDataPrep$new(lfqdata, annotation, dea_config)
  prep$cont_decoy_summary()
  # The example is already at feature level: nothing to aggregate.
  suppressWarnings(prep$aggregate())
  prep$transform_data()
  dea <- prep$build_deanalyse(contrasts)
  dea$build_default()
  dea$get_annotated_contrasts()
  dea$filter_contrasts()
  se <- prolfquapp::DEAReportGenerator$new(dea, dea_config, name = "")$make_SummarizedExperiment()
  results_dir <- file.path(dea_dir, "Results_WU_example")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  prolfquapp::write_summarized_experiment_h5ad(se, file.path(results_dir, "AnnData.h5ad"))
  dea_dir
}

# Spreads features over a realistic log2 dynamic range, deterministically by
# position, so that some are faint enough to go undetected.
.example_base_abundance <- function(ids, levels) {
  seq(14, 26, length.out = length(levels))[match(ids, levels)]
}

# A matched pair of example DEA output directories, built once per session
# because the modelling is the expensive part.
example_dea_pair <- function() {
  root <- file.path(tempdir(), "prophosqua_example_dea")
  paths <- list(phospho = file.path(root, "DEA_phospho"), protein = file.path(root, "DEA_protein"))
  if (file.exists(file.path(paths$protein, "Results_WU_example", "AnnData.h5ad"))) {
    return(paths)
  }

  samples <- .example_dea_samples()
  proteins <- paste0("P", seq_len(8))
  sites <- paste0(rep(proteins, each = 2), "~S", c(10, 20))
  group_of <- function(name) samples$G_[match(name, samples$Name)]
  # Deterministic values rather than rnorm(): an example must not depend on, or
  # disturb, the session's RNG state. Low-abundance sites go missing in some
  # samples, because that is what an imputing model fits its dropout on.
  wobble <- function(n) rep_len(c(-0.20, 0.05, 0.15, -0.10, 0.25, -0.15), n)

  long_protein <- expand.grid(Name = samples$Name, protein_Id = proteins, stringsAsFactors = FALSE)
  long_protein$G_ <- group_of(long_protein$Name)
  long_protein$normalized_abundance <- .example_base_abundance(long_protein$protein_Id, proteins) +
    wobble(nrow(long_protein))

  long_site <- expand.grid(Name = samples$Name, site = sites, stringsAsFactors = FALSE)
  long_site$protein_Id <- sub("~.*", "", long_site$site)
  long_site$G_ <- group_of(long_site$Name)
  # A shift in one group only, so that the contrasts have something to find.
  long_site$normalized_abundance <- .example_base_abundance(long_site$site, sites) +
    wobble(nrow(long_site)) +
    ifelse(long_site$G_ == "a", 0.8, 0)
  faint <- long_site$normalized_abundance < stats::quantile(long_site$normalized_abundance, 0.3, na.rm = TRUE)
  long_site$normalized_abundance[faint & seq_len(nrow(long_site)) %% 3 == 0] <- NA_real_

  protein_annotation <- data.frame(
    protein_Id = proteins,
    description = paste("protein", proteins),
    gene_name = sub("^P", "GENE", proteins),
    protein_length = 300L,
    stringsAsFactors = FALSE
  )
  # A PTM reader keys its row annotation on protein and site; the site rows
  # carry the protein columns beside the site ones.
  site_annotation <- data.frame(
    protein_Id = sub("~.*", "", sites),
    site = sites,
    posInProtein = 10L,
    modAA = "S",
    SequenceWindow = "AAAAAAASAAAAAAA",
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(protein_annotation, by = "protein_Id")

  contrasts <- c(a_vs_b = "G_a - G_b")
  .example_dea_dir(paths$phospho, long_site, list(protein_Id = "protein_Id", site = "site"), site_annotation, contrasts)
  .example_dea_dir(paths$protein, long_protein, list(protein_Id = "protein_Id"), protein_annotation, contrasts)
  paths
}
