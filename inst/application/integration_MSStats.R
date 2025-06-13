# integration_MSStats.R
# Integration of phospho-DIA data analysis
# Author: Functional Genomics Center Zurich
# Description: This script performs integration analysis of phospho-DIA data,
# combining phospho-peptide and total proteome measurements.
rm(list=ls())
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
fracti <- "Integration_filtered"

# Date and directory settings
datetoday <- format(Sys.Date(), "%Y%m%d")
descri <- wu_id
res_dir <- paste0("V3_", fgcz_project, "_", datetoday, "_", fracti, "_", descri)

if (!dir.exists(res_dir)) {
  dir.create(res_dir, recursive = TRUE)
}


path <- "."

# Directory paths
tot_dir <- "DEA_20250528_O38194_WU38194_complete_proteome_filtered_vsn/"
ptm_dir <- "DEA_20250528_O38194_WU38194_enriched_filtered_vsn/"

# File patterns
tot_pattern <- "DE_.*\\.xlsx"
ptm_pattern <- "DE_.*\\.xlsx$"


# Analysis parameters
fdr_threshold <- 0.01
join_column <- c("fasta.id" = "protein_Id", "contrast", "description", "protein_length")
suffix_a <- ".site"
suffix_b <- ".protein"

source("helpers.R")

tot_res <- load_diff_data(tot_dir, pattern = tot_pattern, description = "total proteome")
phospho_res <- load_diff_data(ptm_dir, pattern = ptm_pattern, description = "ptm peptides")


# Process phospho sites
# -------------------
phospho_res2 <- phospho_res |>
  separate(site, c("proteinID", "SiteNlocation", "sequence"),
           remove = FALSE, sep = "_|~"
  )
phospho_res2$site <- gsub("~.*", "", phospho_res2$site)

phospho_res2 <- phospho_res2 |> mutate(
  posInProtein = as.numeric(gsub("[S|T|Y](\\d+)", "\\1", SiteNlocation)),
  modAA = gsub(pattern = "([S|T|Y])\\d+", replacement = "\\1", SiteNlocation),
  #startModSite = posInProtein - 1,
  #endModSite = posInProtein + 1,
  #AllLocalized = TRUE
)

combined_site_prot <- dplyr::left_join(
  phospho_res2,
  tot_res,
  by = join_column,
  suffix = c(suffix_a, suffix_b)
)


# Save combined results
# -------------------
excel_result_list <- list(combinedStats = combined_site_prot)
writexl::write_xlsx(
  excel_result_list,
  path = file.path(res_dir, "Result_phosphoAndTotal.xlsx")
)



make_n_to_c_plots(combined_site_prot)

#### Integration data

# do diff-diff and additional propagation of errors as suggested by MSstatsPTM
combined_test_diff <- prophosqua::test_diff(phospho_res2, tot_res, join_column = join_column)

drumm <- prolfquapp::make_DEA_config_R6(
  PROJECTID = fgcz_project,
  ORDERID = oid_fgcz,
  WORKUNITID = descri
)

prophosqua::copy_phospho_integration()

# render html
rmarkdown::render("_Overview_PhosphoAndIntegration_site.Rmd",
                  params = list(
                    data = combined_test_diff,
                    grp = drumm,
                    phosres = phospho_res
                  ),
                  output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE)
)
file.copy(
  from = "_Overview_PhosphoAndIntegration_site.html",
  to = file.path(res_dir, "Result_phosphoAndIntegration.html")
)



# write to excel
excel_result_list <- list()
excel_result_list$combinedStats <- combined_test_diff
nrow(combined_test_diff)
writexl::write_xlsx(excel_result_list, path = file.path(res_dir, "Result_phosphoAndTotalIntegration.xlsx"))
# Function to determine the significance for plotting NtoC
# only proteins where in at least one contrast the protein is significant (fdrThreshold) are plotted
fdr_threshold <- 0.01
# build up candidate matrix for plotting
candidate_mat <- combined_test_diff[!is.na(combined_test_diff$FDR_I), ]
cand <- candidate_mat[candidate_mat$FDR_I < fdr_threshold, ]

# proteins with sites regulated in any of the contrasts.
sig_protein_hits <- unique(cand$protein_Id)

# look at one of the contrasts
combo_mat <- candidate_mat
combo_mat <- combo_mat |> dplyr::filter(protein_Id %in% sig_protein_hits)
# how many true are ok -> position parsed properly -> this is irrelevant for site centric approach here
table(!is.na(combo_mat$posInProtein))
colnames(combo_mat)

make_n_to_c_plots_integrated(combo_mat)
