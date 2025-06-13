# integration_MSStats_improved.R
# Integration of phospho-DIA data analysis
# Author: Functional Genomics Center Zurich
# Description: This script performs integration analysis of phospho-DIA data,
# combining phospho-peptide and total proteome measurements.
# This is an improved version with better code organization and DRY principles.

# Load required libraries
library(tidyverse)
library(prolfqua)
library(prophosqua)
library(dplyr)
library(stringr)
library(seqinr)

# Print version information
message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")

# Configuration parameters
# ----------------------
# Project identifiers
wu_id <- "CustomPhosphoAnalysis"
fgcz_project <- "p38194"
oid_fgcz <- "fgcz_project"
fracti <- "Integration"

# Date and directory settings
datetoday <- format(Sys.Date(), "%Y%m%d")
descri <- wu_id
res_dir <- paste0("V3_", fgcz_project, "_", datetoday, "_", fracti, "_", descri)
path <- "."

# Directory paths
tot_dir <- "DEA_20250516_O38194_WU38194_complete_proteome_vsn/"
ptm_dir <- "DEA_20250516_O38194_WU38194_enriched_vsn/"

# File patterns
tot_pattern <- "DE_.*_complete_proteome\\.xlsx$"
ptm_pattern <- "DE_.*_enriched\\.xlsx$"

# Analysis parameters
fdr_threshold <- 0.01
join_column <- c("fasta.id" = "protein_Id", "contrast", "description", "protein_length")
suffix_a <- ".site"
suffix_b <- ".protein"

# Required columns for validation
required_cols <- c("protein_Id", "protein_length", "contrast")

# Excluded contrasts
excluded_contrasts <- c("ask13FCvsWTFC", "H6FCvsH1FC")

# Helper functions
# ---------------
find_file_by_pattern <- function(dir, pattern, description) {
  file <- list.files(
    path = dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )[1]

  if (is.na(file)) {
    stop("Could not find ", description, " file matching pattern '", pattern, "' in ", dir)
  }

  file
}

filter_contaminants <- function(data, protein_id_col = "protein_Id") {
  data |>
    dplyr::filter(!grepl("FGCZCont", !!sym(protein_id_col))) |>
    dplyr::filter(!grepl("contam_", !!sym(protein_id_col))) |>
    dplyr::filter(!grepl("^rev_", !!sym(protein_id_col)))
}

load_and_preprocess_data <- function(file_path, sheet_name = "diff_exp_analysis") {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  data <- readxl::read_xlsx(path = file_path, sheet = sheet_name)

  # Validate required columns
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Filter contaminants and excluded contrasts
  data |>
    filter_contaminants() |>
    dplyr::filter(!.data$contrast %in% excluded_contrasts)
}

process_phospho_sites <- function(data) {
  data |>
    tidyr::separate(.data$site, c("proteinID", "SiteNlocation", "sequence"),
      remove = FALSE, sep = "_|~"
    ) |>
    mutate(
      posInProtein = as.numeric(gsub("[S|T|Y](\\d+)", "\\1", .data$SiteNlocation)),
      modAA = gsub(pattern = "([S|T|Y])\\d+", replacement = "\\1", .data$SiteNlocation),
      startModSite = .data$posInProtein - 1,
      endModSite = .data$posInProtein + 1,
      AllLocalized = TRUE
    )
}

create_n_to_c_plots <- function(data, plot_type = c("initial", "integrated")) {
  plot_type <- match.arg(plot_type)
  plot_fn <- if (plot_type == "initial") n_to_c_plot else n_to_c_plot_integrated

  # Select required columns based on plot type
  required_cols <- if (plot_type == "initial") {
    c(
      "protein_Id", "contrast", "protein_length", "site", "diff.protein",
      "diff.site", "FDR.site", "posInProtein", "startModSite", "endModSite",
      "AllLocalized", "modAA", "modelName.site"
    )
  } else {
    c(
      "protein_Id", "contrast", "protein_length", "site", "diff_diff",
      "FDR_I", "posInProtein", "modAA"
    )
  }

  # Prepare data for plotting
  plot_data <- data |>
    dplyr::select(all_of(required_cols)) |>
    mutate(site = gsub("~.*", "", site)) |>
    group_by(.data$protein_Id, .data$contrast, .data$protein_length) |>
    tidyr::nest()
  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))
  # Create plots
  for (i in seq_len(nrow(plot_data))) {
    plot_data$plot[[i]] <- plot_fn(
      plot_data$data[[i]],
      plot_data$protein_Id[[i]],
      plot_data$protein_length[[i]],
      plot_data$contrast[[i]]
    )
  }

  # Group by contrast and create PDFs
  plot_data <- ungroup(plot_data)
  plot_data_wr <- dplyr::select(plot_data, .data$contrast, .data$plot, .data$data)
  contrast_groups <- plot_data_wr |>
    group_by(.data$contrast) |>
    group_split()

  # Create PDFs for each contrast
  for (group in contrast_groups) {
    contrast_name <- unique(group$contrast)
    pdf_fn <- paste0(
      "Significant_Sites_",
      if (plot_type == "integrated") "Normalized" else "Proteins_",
      contrast_name,
      "_NtoCplots.pdf"
    )

    pdf(file.path(res_dir, pdf_fn))
    for (i in seq_len(nrow(group))) {
      print(group$plot[[i]])
      grid::grid.newpage()

      # Prepare table data
      table_data <- if (plot_type == "initial") {
        group$data[[i]] |>
          select(-all_of(c("startModSite", "endModSite", "AllLocalized")))
      } else {
        group$data[[i]]
      }

      table_grob <- gridExtra::tableGrob(
        table_data,
        theme = gridExtra::ttheme_default(base_size = 6)
      )
      grid::grid.draw(table_grob)
    }
    dev.off()
  }
}

# Main execution
# -------------
# Find input files
tot_xlsx <- find_file_by_pattern(
  dir = tot_dir,
  pattern = tot_pattern,
  description = "total proteome"
)

phospho_xlsx <- find_file_by_pattern(
  dir = ptm_dir,
  pattern = ptm_pattern,
  description = "phospho"
)

# Load and preprocess data
tot_res <- load_and_preprocess_data(tot_xlsx)
phospho_res <- load_and_preprocess_data(phospho_xlsx)

# Process phospho sites
phospho_res2 <- process_phospho_sites(phospho_res)

# Join total and phospho analysis
combined_site_prot <- dplyr::left_join(
  phospho_res2,
  tot_res,
  by = join_column,
  suffix = c(suffix_a, suffix_b)
)

# Create results directory
if (!dir.exists(res_dir)) {
  dir.create(res_dir, recursive = TRUE)
}

# Save combined results
excel_result_list <- list(combinedStats = combined_site_prot)
writexl::write_xlsx(
  excel_result_list,
  path = file.path(res_dir, "Result_phosphoAndTotal.xlsx")
)

# Create N-to-C plots
debug(create_n_to_c_plots)
create_n_to_c_plots(combined_site_prot, "initial")

# Integration analysis
combined_test_diff <- prophosqua::test_diff(phospho_res2, tot_res, join_column = join_column)

# Create DEA configuration
drumm <- prolfquapp::make_DEA_config_R6(
  PROJECTID = fgcz_project,
  ORDERID = oid_fgcz,
  WORKUNITID = descri
)

# Copy integration files
prophosqua::copy_phospho_integration()

# Render HTML report
rmarkdown::render("_Overview_PhosphoAndIntegration_site.Rmd",
  params = list(
    data = combined_test_diff,
    grp = drumm,
    phosres = phospho_res
  ),
  output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE)
)

# Copy HTML report to results directory
file.copy(
  from = "_Overview_PhosphoAndIntegration_site.html",
  to = file.path(res_dir, "Result_phosphoAndIntegration.html")
)

# Save integration results
excel_result_list <- list(combinedStats = combined_test_diff)
writexl::write_xlsx(
  excel_result_list,
  path = file.path(res_dir, "Result_phosphoAndTotalIntegration.xlsx")
)

# Create integrated N-to-C plots
create_n_to_c_plots(combined_test_diff, "integrated")
