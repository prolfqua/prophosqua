#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#

# libraries necessairy
library(prolfquapp)
library(prophosqua)
library(tidyverse)
library(ggseqlogo)


# from prophosqua chapter and witolds AnalysePTM.Rmd
# for console:
# cd /Users/jonasgrossmann/FGCZ/Interdisziplinar/o39871_HelanHeLa/Hela/prophosqua
# R --vanilla -e "library(prolfquapp); copy_shell_script(workdir ='.')"

# input for total is msstats.csv
# ./prolfqua_dataset.sh --indir ../o39871_FP231_Hela_input/ --dataset p39871_dataset_Hela_Total.xlsx --software prolfquapp.MSSTATS

# read back in and generate contrasts in R
library(openxlsx)
annot <- read.xlsx("p39871_dataset.xlsx", sheet = "Sheet1")
knitr::kable(annot)

# ./prolfqua_yaml.sh -y p39871_config_Hela_total.yaml

# run DEA for the total!
# only 2 Groups -> using CONTROL column to define ratio!
# ./prolfqua_dea.sh -i ../o39871_FP231_Hela_input/ -d p39871_dataset.xlsx -y p39871_config.yaml -w _Hela_TOTAL -s prolfquapp.MSSTATS

# go for PTM
# for PTM we read *_combined_STY.tsv file -> this is not yet filtered for BestLocalizationSite Probability
# copy and rename this file first and then filter and overwrite the original file.

ptmDir <- "../o39871_FP231_Hela_enriched/"
ptmFile <- dir(path = ptmDir, pattern = "^combined_site_STY_.+\\.tsv",
                     recursive = TRUE, full.names = TRUE)
ptmFile
# copy and save wit new name at the same place bak.*
file.copy(from = ptmFile, to = file.path(dirname(ptmFile), "bak.Original_noFilter_combined_site_STY_79.9663.tsv"), overwrite = TRUE)

# Do filtering on best_localization_probability and overwrite the original file
library(readr)
LocalizationThreshold <- 0.75
ptmData <- read_tsv(ptmFile)
colnames(ptmData)
filteredPtmData <- ptmData %>%
    filter(`Best Localization Probability` >= LocalizationThreshold)
# overwrite original file
write_tsv(filteredPtmData, ptmFile)


# dataset again first
# ./prolfqua_dataset.sh -i ../o39871_FP231_Hela_enriched/ -d p39871_dataset_PTM.xlsx -s prolfquappPTMreaders.FP_combined_STY

# ./prolfqua_yaml.sh -y p39871_config_PTM.yaml

# then we run DEA for the phospho!
# ./prolfqua_dea.sh -i ../o39871_FP231_Hela_enriched -d p39871_dataset_PTM.xlsx -y p39871_config_PTM.yaml -w _Hela_BLPfilt_PTM -s prolfquappPTMreaders.FP_combined_STY


# let's rock it with prophosqua
library(prolfquapp)
library(prophosqua)
library(openxlsx)
library(rlang)
library(dplyr)

ptm_feature <- "singlesite"
wu_id <- "Hela_Phospho_wBLPfiltering_Hela_pTJ41_vs_Hela_pTJ40"
fgcz_project <- "p39819"
oid_fgcz <- "o39819"
datetoday <- format(Sys.Date(), "%Y%m%d") #
project_dir <- "."


descri <- wu_id
res_dir <- file.path(
    project_dir,
    paste0(fgcz_project, "_", datetoday, "_", ptm_feature, "_", descri)
)
if (!dir.exists(res_dir)) {
    dir.create(res_dir, recursive = TRUE)
}

res_dir

#
tot_file <-  "DEA_20251030_WU_Hela_TOTAL_vsn/Results_WU__Hela_TOTAL/DE_WU_Hela_TOTAL.xlsx"
ptm_file <- "DEA_20251030_WU_Hela_BLPfilt_PTM_vsn/Results_WU__Hela_BLPfilt_PTM/DE_WU_Hela_BLPfilt_PTM.xlsx"

stopifnot(file.exists(tot_file))
stopifnot(file.exists(ptm_file))

required_cols <- c("protein_Id", "protein_length", "contrast")
tot_res <- load_and_preprocess_data(tot_file, required_cols)
tot_res <- filter_contaminants(tot_res)
phospho_res <- load_and_preprocess_data(ptm_file, required_cols)
phospho_res <- filter_contaminants(phospho_res)

join_column <- c(
    "fasta.id" = "protein_Id",
    "contrast",
    "description",
    "protein_length"
)
suffix_a <- ".site"
suffix_b <- ".protein"
# combine
combined_site_prot <- dplyr::left_join(
    phospho_res,
    tot_res,
    by = join_column,
    suffix = c(suffix_a, suffix_b)
)

combined_test_diff <- prophosqua::test_diff(phospho_res, tot_res, join_column = join_column)

unique(combined_test_diff$contrast)

# some issues fixed with fixing the functions
source("myRhelpers.R")
library(ggplot2)

#rm(plot_data)

# calculating and directly plotting n-to-c-usage plots for all contrasts
for (i in 1:length(unique(combined_test_diff$contrast))) {
    print(unique(combined_test_diff$contrast)[i])
    # plotting
    if (ptm_feature== "multisite") {
        combined_test_diff_long <- explode_multisites(combined_test_diff)
        plot_data <- n_to_c_usage(combined_test_diff_long, unique(combined_test_diff$contrast)[i], FDR_threshold = 0.1)
    } else if (ptm_feature== "singlesite") {
        plot_data <- n_to_c_usage(combined_test_diff, unique(combined_test_diff$contrast)[i], FDR_threshold = 0.1)
    }
    # N-to-C plot
    pdfFN <- paste("Site_differential_UsageChange_", unique(combined_test_diff$contrast)[i], ".pdf", sep = "")
    pdf(file.path(res_dir, pdfFN))
    for (jj in seq_len(length(plot_data$protein_Id))) {
        #print(i)
        p <- n_to_c_plot_integrated(poi_matrix_min = plot_data$data[[jj]],
                                    protein_name = plot_data$protein_Id[jj],
                                    prot_length = plot_data$protein_length[jj],
                                    contrast = plot_data$contrast[[jj]], thr_a = 0.01, thr_b = 0.1)
        print(p)

    }
    dev.off()
}

# # plotting
# if (ptm_feature== "multisite") {
#     combined_test_diff_long <- explode_multisites(combined_test_diff)
#     plot_data <- n_to_c_usage(combined_test_diff_long, "DMSO_bFGF_vs_DMSO", FDR_threshold = 0.1)
# } else if (ptm_feature== "singlesite") {
#     plot_data <- n_to_c_usage(combined_test_diff, "DMSO_bFGF_vs_DMSO", FDR_threshold = 0.1)
# }


# plot_data$data <- lapply(plot_data$data, function(df) {
#   df$posInProtein <- suppressWarnings(as.numeric(df$posInProtein))
#   if (any(is.na(df$posInProtein))) {
#     warning("Some posInProtein values could not be converted to numeric in one or more data frames.")
#   }
#   df
# })


# N-to-C plot
# pdfFN <- paste("Site_differential_UsageChange_", unique(combined_test_diff$contrast)[i], ".pdf", sep = "")
# pdf(file.path(res_dir, "Site_differential_UsageChange_DMSO_bFGF_vs_DMSO.pdf"))
# for (i in seq_len(length(plot_data$protein_Id))) {
#     #print(i)
#   p <- n_to_c_plot_integrated(poi_matrix_min = plot_data$data[[i]],
#                          protein_name = plot_data$protein_Id[i],
#                          prot_length = plot_data$protein_length[i],
#                          contrast = "Constricted_vs_Dilated", thr_a = 0.01, thr_b = 0.1)
#   print(p)
#
# }
# dev.off()



# Same also for n-to-c expression



# # plotting
# if (ptm_feature== "multisite") {
#   combined_test_diff_long <- explode_multisites(combined_test_diff)
#   plot_data_expression <- n_to_c_expression(combined_test_diff_long, "DMSO_bFGF_vs_DMSO", FDR_threshold = 0.1)
# } else if (ptm_feature== "singlesite") {
#   plot_data_expression <- n_to_c_expression(combined_test_diff, "DMSO_bFGF_vs_DMSO", FDR_threshold = 0.1)
# }
#
# # N-to-C plot
# pdf(file.path(res_dir, "Site_differential_Expression_DMSO_bFGF_vs_DMSO.pdf"))
# for (i in seq_len(length(plot_data_expression$protein_Id))) {
#   #print(i)
#   p <- n_to_c_plot(poi_matrix_min = plot_data_expression$data[[i]],
#                               protein_name = plot_data_expression$protein_Id[i],
#                               prot_length = plot_data_expression$protein_length[i],
#                               contrast = "DMSO_bFGF_vs_DMSO", thr_a = 0.05, thr_b = 0.2)
#   print(p)
#
# }
# dev.off()
#
for (i in 1:length(unique(combined_test_diff$contrast))) {
    print(unique(combined_test_diff$contrast)[i])
    # plotting
    if (ptm_feature== "multisite") {
        combined_test_diff_long <- explode_multisites(combined_test_diff)
        plot_data <- n_to_c_expression(combined_test_diff_long, unique(combined_test_diff$contrast)[i], FDR_threshold = 0.1)
    } else if (ptm_feature== "singlesite") {
        plot_data <- n_to_c_expression(combined_test_diff, unique(combined_test_diff$contrast)[i], FDR_threshold = 0.1)
    }
    # N-to-C plot
    pdfFN <- paste("Site_differential_Expression_", unique(combined_test_diff$contrast)[i], ".pdf", sep = "")
    pdf(file.path(res_dir, pdfFN))
    for (jj in seq_len(length(plot_data$protein_Id))) {
        #print(i)
        p <- n_to_c_plot(poi_matrix_min = plot_data$data[[jj]],
                                    protein_name = plot_data$protein_Id[jj],
                                    prot_length = plot_data$protein_length[jj],
                                    contrast = plot_data$contrast[[jj]], thr_a = 0.01, thr_b = 0.1)
        print(p)

    }
    dev.off()
}




# rendering reports
drumm <- prolfquapp::make_DEA_config_R6(
  PROJECTID = fgcz_project,
  ORDERID = oid_fgcz,
  WORKUNITID = descri
)
# Copy integration files
prophosqua::copy_phospho_integration()

# rendering
rmarkdown::render(
  "_Overview_PhosphoAndIntegration_site.Rmd",
  params = list(
    #data = combined_test_diff_long,
      data = combined_test_diff,
      grp = drumm
  ),
  output_format = bookdown::html_document2(
    toc = TRUE,
    toc_float = TRUE
  ),
  envir = new.env(parent = globalenv())
)
# Copy HTML report to results directory
file.copy(
  from = "_Overview_PhosphoAndIntegration_site.html",
  to = file.path(
    res_dir,
    "Result_phosphoAndTotalIntegration.html"
  )
)


# Do sequence windows n logos
col2sel <- c("Protein", "modAA", "posInProtein")
unique_prot_pep_seq <- combined_test_diff |> select(all_of(col2sel)) |> distinct()

rev_pattern = "_rev"
fasta_file <- "../o39871_FP231_Hela_input/2025-08-14-decoys-reviewed-contam-UP000005640-spikein.fasta"
fasta <- prolfquapp::get_annot_from_fasta(fasta_file, pattern_decoys = "^rev_", include_seq = TRUE)
unique_prot_pep_seq <- dplyr::inner_join(unique_prot_pep_seq, fasta, by = c(Protein = "fasta.id"))

seq_window <- prophosqua::get_sequence_windows(unique_prot_pep_seq, pos_in_protein = "posInProtein",flank_size = 7)
colnames(seq_window)
seq_window <- seq_window |> dplyr::select(c(col2sel, SequenceWindow = "sequence_window"))
combined_test_diff <- dplyr::inner_join(combined_test_diff, seq_window, by = col2sel )

# write out to excel
excel_result_list <- list(
  combinedStats = combined_test_diff,
  combinedSiteProteinData = combined_site_prot)

writexl::write_xlsx(
  excel_result_list,
  path = file.path(
    res_dir,
    "Result_phosphoAndTotalIntegration.xlsx"
  )
)



library(ggseqlogo)
# Filter for significantly regulated sites (FDR_I < 0.05)
log_2FC_threshold <- 0.58

combined_test_diff$contrast |> unique()
significant_sites <- combined_test_diff |>
  dplyr::filter(
    FDR_I < 0.05,
    BLP > 0.75,
    abs(diff_diff) > log_2FC_threshold,
    !is.na(posInProtein),
    !grepl("^_", SequenceWindow),
    !grepl("_$", SequenceWindow)
  ) |>
  dplyr::mutate(
    regulation = case_when(
      diff_diff > 0 ~ "upregulated",
      diff_diff < 0 ~ "downregulated",
      TRUE ~ "no_change"
    )
  )

with(significant_sites, table(contrast, regulation, modAA)) |>
  knitr::kable(caption = "number of significantly regulated sites with abs(log2FC)>0.58, BLP>0.75, FDR_I < 0.05") |>
  kableExtra::kable_styling(font_size = 9, latex_options = "scale_down")

# Get 7th character from each sequence window
st <- nrow(significant_sites)
significant_sites$seventh_chars <- toupper(substr(significant_sites$SequenceWindow, 8, 8))
significant_sites <- significant_sites |> filter(seventh_chars == modAA)
end <- nrow(significant_sites)
if(st != end){
  warning("sequence window missalignments : ", st - end)
}


# support
library(kableExtra)
library(webshot2)
library(cowplot)
library(png)
mytbl <- with(significant_sites, table(contrast, regulation, modAA)) |>
  knitr::kable(caption = "regulated sites with abs(log2FC)>0.58, BLP>0.75, fdr smaller 0.05") |>
  kableExtra::kable_styling(font_size = 9, latex_options = "scale_down")

tbl_path <- file.path(tempdir(), "sig_table.png")
save_kable(mytbl, tbl_path, zoom = 2)



# plot to pdf
(pdfName <- paste0("Phospho_Site_Logos_by_AA_and_Regulation_", wu_id, ".pdf"))
pdf(file.path(res_dir, pdfName), width = 8, height = 3)

# support
# Plot the table image
grid::grid.newpage()
grid::grid.raster(png::readPNG(tbl_path))

# the plots
seq_list <- significant_sites |> filter(modAA == "S") |>
  mutate(.grp = paste(contrast, regulation,  sep = "_")) |>
  group_by(.grp) |>
  summarize(
    seqs = list(toupper(SequenceWindow)),
    .groups = "drop"
  ) |>
  with(setNames(seqs, .grp))

ggseqlogo(
  seq_list,
  ncol = 2,
  seq_type = "aa",
  method = "probability"
)

seq_list <- significant_sites |> filter(modAA == "T") |>
  mutate(.grp = paste(contrast, regulation,  sep = "_")) |>
  group_by(.grp) |>
  summarize(
    seqs = list(toupper(SequenceWindow)),
    .groups = "drop"
  ) |>
  with(setNames(seqs, .grp))

ggseqlogo(
  seq_list,
  ncol = 2,
  seq_type = "aa",
  method = "probability"
)
seq_list <- significant_sites |> filter(modAA == "Y") |>
  mutate(.grp = paste(contrast, regulation,  sep = "_")) |>
  group_by(.grp) |>
  summarize(
    seqs = list(toupper(SequenceWindow)),
    .groups = "drop"
  ) |>
  with(setNames(seqs, .grp))
ggseqlogo(
  seq_list,
  ncol = 2,
  seq_type = "aa",
  method = "probability"
)

dev.off()
















