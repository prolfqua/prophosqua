#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

#remotes::install_github('wolski/prolfquapp', dependencies = TRUE)
#remotes::install_github('fgcz/prolfqua', dependencies = TRUE)

rm(list=ls())

# libs
library(tidyverse)
library(prolfqua)
library(prolfquapp)
library(proubiqua)
library(dplyr)
library(stringr)
library(openxlsx)
library(seqinr)

# integration of ubi-DIA
message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")


# parameters and thresholds
# variables
WUID <- "Ubi"
fgczProject <- "p37127"
OIDfgcz <- "o37127"
fracti <- "Integration"
datetoday <- format(Sys.Date(), "%Y%m%d")
descri <- WUID # "OneCondition"
(resDir <- paste0(fgczProject, "_",OIDfgcz,"_",datetoday,"_",fracti,"_", descri))

path <- "."


# find path manually
# read back in results
(totXlsx  <-  paste0("DEA_p37127_specifiedContrasts_totalProtein/DEA_20241220_PIp37127_Oo37127_WUspecifiedContrasts_robscale/Results_WU_specifiedContrasts/DE_WUspecifiedContrasts.xlsx"))
(ubiXlsx <- paste0("DEA_20241220_p37127_specifiedContrasts_ubi/Results_WU_specifiedContrasts/DE_WUspecifiedContrasts.xlsx"))


totRes <- readxl::read_xlsx(path = totXlsx, sheet = "diff_exp_analysis")
totRes$protein_length |> is.na() |> mean() # check how many proteins have no length @ WeW how can this happen? these are all the fw proteins.. only rev proteins are properly parsed?
totRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id) # these are the protein IDs without length! # should be only rev_ or proteins that are not found in this fasta

# if smaller case rev_ are used and not fixed in yaml file before running DEA this here helps to still get the integration done
# totRes$protein_Id <- sapply(strsplit((totRes$protein_Id), split = "\\|"), function(x)x[2])

rev_pattern = "^rev_"
totRes <- totRes |> dplyr::filter(!grepl("FGCZCont", protein_Id))
totRes <- totRes |> dplyr::filter(!grepl("contam_", protein_Id))
totRes <- totRes |> dplyr::filter(!grepl(rev_pattern, protein_Id)) # here we take out all revs!


ubiRes <- readxl::read_xlsx(path = ubiXlsx, sheet = "diff_exp_analysis")
ubiRes$protein_length |> is.na() |> mean()
ubiRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id) # these are the protein IDs without length! # should be only rev_

# drop contaminants and rev sequences
ubiRes <- ubiRes |> dplyr::filter(!grepl("FGCZCont", protein_Id))
ubiRes <- ubiRes |> dplyr::filter(!grepl("cont_", protein_Id))
ubiRes <- ubiRes |> dplyr::filter(!grepl(rev_pattern, protein_Id))

fasta_file = "../o36946_FP_TMTi_enriched/2024-12-02-decoys-contam-UP000000589.fasta"

# thats much faster than before thanks @ witold
# source helper functions
source("/Users/jonasgrossmann/GitHub/prophosqua/R/FP_phosphoHelperFunctions_v3_202310.R")
#ubiRes$AllLocalized <- TRUE
#seq_window <- get_sequence_windows(ubiRes, fasta_file, rev_pattern) # this is very much tailored towards phospho analysis

join_column <- c("fasta.id" = "fasta.id", "protein_Id","contrast", "description", "protein_length")

combined_test_diff <- test_diff(ubiRes, totRes, join_column = join_column)

# write out html
resDir
dir.create(resDir)
drumm <- prolfquapp::make_DEA_config_R6(
  PROJECTID = fgczProject,
  ORDERID = OIDfgcz,
  WORKUNITID = descri)

# render html
rmarkdown::render("_Overview_Ubi_Integration_WEW.Rmd", params = list(data = combined_test_diff, grp = drumm, phosres = ubiRes), output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE))
file.copy(from = "_Overview_Ubi_Integration_WEW.html", to = file.path(resDir, "Result_UbiAndIntegration.html"))

# write to excel
excelResultList <- list()
excelResultList$combinedStats <- combined_test_diff
#excelResultList$seq_window <- seq_window
writexl::write_xlsx(excelResultList, path = file.path(resDir,"Result_UbiAndTotalIntegration.xlsx"))
#

#  N-to-C plotting is not yet working
# # Function to determine the significance for plotting NtoC
# fdrThreshold = 0.00001
#
#
# candidateMat <- combined_test_diff[!is.na(combined_test_diff$FDR.site), ]
# cand <- candidateMat[candidateMat$FDR.site < fdrThreshold ,]
# nrow(cand)
#
# # proteins with sites regulated in any of the contrasts.
# mySigProteinHits <- unique(cand$protein_Id)
# length(mySigProteinHits)
#
# # look at one of the contrasts
# comboMat <- candidateMat
# comboMat <- comboMat |> dplyr::filter(protein_Id %in% mySigProteinHits)
#
# # parse the position in protein
# comboMat$posInProtein <- as.numeric(gsub(x = gsub(x = sapply(strsplit((comboMat$site), split = "~"), function(x)x[2]), pattern = "\\(K", replacement = ""), pattern = "\\)", replacement = ""))
# table(!is.na(comboMat$posInProtein))
#
# # select columns necessary for plotting
# comboMat_min <- dplyr::select(
#   comboMat,
#   c("protein_Id",
#     "contrast",
#     "protein_length",
#     "site",
#     "diff.protein",
#     "diff.site",
#     "FDR.site",
#     "posInProtein",
#     model_site = "modelName.site"
#   ))
#
#
# comboMat_min <- comboMat_min |> dplyr::group_by(protein_Id, contrast, protein_length) |> tidyr::nest()
# comboMat_min$plot <- vector(mode = "list", nrow(comboMat_min))
#
# # how many ubi proteins at the beginning?
# length(unique(ubiRes$protein_Id))
#
#
# # check how many pages are plotted
# length(unique(comboMat_min$protein_Id))
#
# # fill all slots with plots
# for (i in 1:nrow(comboMat_min)) {
#   print(i)
#   comboMat_min$plot[[i]] <- prophosqua::N_to_C_plot(comboMat_min$data[[i]],
#                                                     comboMat_min$protein_Id[[i]],
#                                                     comboMat_min$protein_length[[i]],
#                                                     comboMat_min$contrast[[i]])
# }
#
#
# # Do the plotting only for each contrast individually
# for (j in 1:length(unique(comboMat_min$contrast))) {
#   print(j)
#   oneC_comboMat <- comboMat_min[comboMat_min$contrast == unique(comboMat_min$contrast)[j],]
#   pdfFN <- paste0("SignificantProteins_",unique(comboMat_min$contrast)[j],"_NtoCplots.pdf")
#   pdf(file.path(resDir, pdfFN))
#   for (i in 1:nrow(oneC_comboMat)) {
#     print(oneC_comboMat$plot[[i]])
#     grid::grid.newpage()
#     table <- oneC_comboMat$data[[i]]
#     table <- table |> select(-all_of(c( "startModSite", "endModSite", "AllLocalized")))
#     table_grob <- gridExtra::tableGrob(table, theme = gridExtra::ttheme_default(base_size=6))
#     grid::grid.draw(table_grob)
#   }
#   dev.off()
# }
#
#
# # all in one pdf
# # pdf(file.path(resDir, "NtoCplots2.pdf"))
# # for (i in 1:nrow(comboMat_min)) {
# #   print(comboMat_min$plot[[i]])
# #   grid::grid.newpage()
# #   table <- comboMat_min$data[[i]]
# #   table <- table |> select(-all_of(c( "startModSite", "endModSite", "AllLocalized")))
# #   table_grob <- gridExtra::tableGrob(table, theme = gridExtra::ttheme_default(base_size=6))
# #   grid::grid.draw(table_grob)
# # }
# # dev.off()
#
#
#
#
# #
# #rFN <- paste0(",")
# save.image(paste0(resDir, ".RData"))
#


