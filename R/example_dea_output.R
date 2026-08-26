# A synthetic pair of prolfquapp DEA output directories.
#
# The DPA/DPU and CorrectFirst reports read a finished DEA run: a phospho one at
# site level and a total-proteome one at protein level. To render either report
# without a pipeline run, the input has to exist in the shape prolfquapp writes
# it -- a `Results_WU_*` directory holding the DE workbook, the normalized
# abundances as parquet, and the analysis configuration as yaml -- so these
# builders produce exactly that, and the example then calls the same compute
# function the pipeline calls.

#' Samples of the Example DEA Run
#'
#' Two groups of three replicates, `b` the control.
#'
#' @return data.frame with `Name`, `raw_file`, `G_` and `Control`.
#' @keywords internal
.example_dea_samples <- function() {
  samples <- expand.grid(
    replicate = 1:3,
    G_ = c("a", "b"),
    stringsAsFactors = FALSE
  )
  samples$Name <- paste0(samples$G_, samples$replicate)
  samples$raw_file <- paste0(samples$Name, ".raw")
  samples$Control <- ifelse(samples$G_ == "b", "C", "T")
  samples
}

#' Analysis Configuration of the Example DEA Run
#'
#' `sample_name` and `file_name` must name different columns: `setup_analysis()`
#' completes the table by joining on both, and one column serving as both ends
#' up duplicated as `Name...1` and `Name...2`.
#'
#' @param hierarchy Named list of hierarchy keys, e.g.
#'   `list(protein_Id = "protein_Id", site = "site")`.
#' @return A [prolfqua::AnalysisConfiguration].
#' @keywords internal
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
  config$file_name <- "raw_file"
  config$set_response("normalized_abundance")
  config
}

#' Write One Example DEA Output Directory
#'
#' @param dea_dir Directory to create.
#' @param long Long-format abundances.
#' @param config Its [prolfqua::AnalysisConfiguration].
#' @param diff_exp The `diff_exp_analysis` sheet.
#' @param normalized The `normalized_abundances` sheet, or NULL.
#' @return `dea_dir`.
#' @keywords internal
.example_dea_dir <- function(dea_dir, long, config, diff_exp, normalized = NULL) {
  # Columns prolfquapp's output carries and setup_analysis() otherwise invents
  # with a warning; naming them keeps the example's shape honest and the test
  # output readable.
  long$isotopeLabel <- "light"
  long$qValue <- 0
  long$nr_children <- 1L

  results_dir <- file.path(dea_dir, "Results_WU_example")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  sheets <- list(diff_exp_analysis = diff_exp)
  if (!is.null(normalized)) {
    sheets$normalized_abundances <- normalized
  }
  writexl::write_xlsx(sheets, file.path(results_dir, "DE_example.xlsx"))

  arrow::write_parquet(
    prolfqua::setup_analysis(long, config),
    file.path(results_dir, "lfqdata_normalized.parquet")
  )
  yaml::write_yaml(
    prolfqua::R6_extract_values(config),
    file.path(results_dir, "lfqdata.yaml")
  )
  dea_dir
}

#' Base Abundance of an Example Feature
#'
#' Spreads features over a realistic log2 dynamic range, deterministically by
#' position, so that some are faint enough to go undetected.
#'
#' @param ids The feature ids of each row.
#' @param levels All feature ids, in order.
#' @return Numeric vector of base abundances.
#' @keywords internal
.example_base_abundance <- function(ids, levels) {
  base <- seq(14, 26, length.out = length(levels))
  base[match(ids, levels)]
}

#' A Matched Pair of Example DEA Output Directories
#'
#' Built once per session, because both reports that read it may be rendered in
#' the same session and the modelling is the expensive part.
#'
#' @return List with `phospho`, `protein` and `annot_file`.
#' @keywords internal
example_dea_pair <- function() {
  root <- file.path(tempdir(), "prophosqua_example_dea")
  marker <- file.path(root, "annotation.tsv")
  paths <- list(
    phospho = file.path(root, "DEA_phospho"),
    protein = file.path(root, "DEA_protein"),
    annot_file = marker
  )
  if (file.exists(marker)) {
    return(paths)
  }

  samples <- .example_dea_samples()
  proteins <- paste0("P", seq_len(8))
  sites <- paste0(rep(proteins, each = 2), "~S", c(10, 20))
  group_of <- function(name) samples$G_[match(name, samples$Name)]

  # Deterministic values rather than rnorm(): an example must not depend on, or
  # disturb, the session's RNG state. The abundances span a realistic dynamic
  # range, and low-abundance sites go missing in some samples, because that is
  # what an imputing aggregator fits its dropout model on -- given a complete
  # matrix of near-identical values it has nothing to fit.
  wobble <- function(n) rep_len(c(-0.20, 0.05, 0.15, -0.10, 0.25, -0.15), n)

  long_protein <- expand.grid(
    Name = samples$Name,
    protein_Id = proteins,
    stringsAsFactors = FALSE
  )
  long_protein$raw_file <- paste0(long_protein$Name, ".raw")
  long_protein$G_ <- group_of(long_protein$Name)
  long_protein$normalized_abundance <- .example_base_abundance(
    long_protein$protein_Id,
    proteins
  ) +
    wobble(nrow(long_protein))

  long_site <- expand.grid(
    Name = samples$Name,
    site = sites,
    stringsAsFactors = FALSE
  )
  long_site$protein_Id <- sub("~.*", "", long_site$site)
  long_site$raw_file <- paste0(long_site$Name, ".raw")
  long_site$G_ <- group_of(long_site$Name)
  # A shift in one group only, so that the contrasts have something to find.
  long_site$normalized_abundance <- .example_base_abundance(long_site$site, sites) +
    wobble(nrow(long_site)) +
    ifelse(long_site$G_ == "a", 0.8, 0)
  # Missing where the signal is weakest, one sample in three, so that detection
  # depends on abundance the way it does in a real run.
  faint <- long_site$normalized_abundance <
    stats::quantile(
      long_site$normalized_abundance,
      0.3,
      na.rm = TRUE
    )
  long_site$normalized_abundance[faint & seq_len(nrow(long_site)) %% 3 == 0] <- NA_real_

  protein_de <- data.frame(
    protein_Id = proteins,
    contrast = "a_vs_b",
    description = paste("protein", proteins),
    gene_name = sub("^P", "GENE", proteins),
    protein_length = 300L,
    diff = seq(-1, 1, length.out = length(proteins)),
    std.error = 0.2,
    df = 6,
    statistic = seq(-3, 3, length.out = length(proteins)),
    FDR = 0.01,
    estimate_type = "observed",
    stringsAsFactors = FALSE
  )
  site_de <- protein_de[rep(seq_len(nrow(protein_de)), each = 2), ]
  site_de$site <- sites
  # A PTM reader keys its row annotation on protein and site, so the site
  # annotation reaches every annotated sheet, diff_exp_analysis included; the
  # compute functions read it from there.
  site_de$posInProtein <- 10L
  site_de$modAA <- "S"
  site_de$SequenceWindow <- "AAAAAAASAAAAAAA"
  site_de$diff <- site_de$diff + 0.5

  site_annotation <- data.frame(
    site = sites,
    posInProtein = 10L,
    modAA = "S",
    SequenceWindow = "AAAAAAASAAAAAAA",
    protein_Id = sub("~.*", "", sites),
    gene_name = sub("^P", "GENE", sub("~.*", "", sites)),
    protein_length = 300L,
    stringsAsFactors = FALSE
  )

  .example_dea_dir(
    paths$phospho,
    long_site,
    .example_dea_config(list(protein_Id = "protein_Id", site = "site")),
    site_de,
    site_annotation
  )
  .example_dea_dir(
    paths$protein,
    long_protein,
    .example_dea_config(list(protein_Id = "protein_Id")),
    protein_de
  )
  readr::write_tsv(
    samples[, c("Name", "G_", "Control")] |> dplyr::rename(Group = "G_"),
    marker
  )
  paths
}
