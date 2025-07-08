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

wu_id <- "CustomPhosphoAnalysis"
fgcz_project <- "PTM_example"
oid_fgcz <- "fgcz_project"
ptm_feature <- "singlesite"

tot_dir <- "/Users/witoldwolski/Dropbox/DataAnalysis/PTM/PTM_example/data_total/DEA_20250630_WUtotal_proteome_vsn/"
tot_pattern <- "DE_.*total_proteome\\.xlsx$"

if (ptm_feature == "multisite") {
  ptm_dir <- "/Users/witoldwolski/Dropbox/DataAnalysis/PTM/PTM_example/data_ptm/DEA_20250630_WUmultisite_PTM_vsn/"
  ptm_pattern <- "DE_.*multisite_PTM\\.xlsx$"
} else if (ptm_feature == "singlesite") {
  ptm_dir <- "/Users/witoldwolski/Dropbox/DataAnalysis/PTM/PTM_example/data_ptm/DEA_20250701_WUsinglesite_PTM_vsn/"
  ptm_pattern <- "DE_.*singlesite_PTM\\.xlsx$"
}

# Date and directory settings
datetoday <- format(Sys.Date(), "%Y%m%d")
descri <- wu_id
res_dir <- paste0("user_", fgcz_project, "_", datetoday, "_", ptm_feature, "_", descri)
path <- "."


# File patterns

# Analysis parameters
fdr_threshold <- 0.01
join_column <- c("fasta.id" = "protein_Id", "contrast", "description", "protein_length")
suffix_a <- ".site"
suffix_b <- ".protein"
required_cols <- c("protein_Id", "protein_length", "contrast")

# Load helper functions
source("helper_functions.r")

# Main execution
# -------------
# Find input files
tot_xlsx <- find_file_by_pattern(
  dir = tot_dir,
  pattern = tot_pattern,
  description = "total_proteome"
)
tot_res <- load_and_preprocess_data(tot_xlsx, required_cols)
tot_res <- filter_contaminants(tot_res)

phospho_xlsx <- find_file_by_pattern(
  dir = ptm_dir,
  pattern = ptm_pattern,
  description = "phospho"
)
phospho_res <- load_and_preprocess_data(phospho_xlsx, required_cols)
phospho_res <- filter_contaminants(phospho_res)

# Join total and phospho analysis
combined_site_prot <- dplyr::left_join(
  phospho_res,
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
  path = file.path(res_dir, "Result_phosphoAndTotal_expression.xlsx")
)

if (ptm_feature == "multisite") {
  combined_site_prot_long <- explode_multisites(combined_site_prot)
  plot_data <- n_to_c_expression(combined_site_prot_long, "Late_vs_Uninfect_at_KO", 0.01)
} else if (ptm_feature == "singlesite") {
  plot_data <- n_to_c_expression(combined_site_prot, "Late_vs_Uninfect_at_KO", 0.01)
}
# Create plots
pdf(file.path(res_dir, "Site_differential_Expression.pdf"))
for (i in seq_len(100)) {
  print(plot_data$plot[[i]])
}
dev.off()


# Integration analysis
combined_test_diff <- prophosqua::test_diff(phospho_res, tot_res, join_column = join_column)
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

if (ptm_feature == "multisite") {
  combined_test_diff_long <- explode_multisites(combined_test_diff)
  plot_data <- n_to_c_usage(combined_test_diff_long, "KO_vs_WT_at_Uninfect", 0.01)
} else if (ptm_feature == "singlesite") {
  plot_data <- n_to_c_usage(combined_test_diff, "KO_vs_WT_at_Uninfect", 0.01)
}

pdf(file.path(res_dir, "Site_differential_OccupancyChange.pdf"))
for (i in seq_len(100)) {
  print(plot_data$plot[[i]])
}
dev.off()
