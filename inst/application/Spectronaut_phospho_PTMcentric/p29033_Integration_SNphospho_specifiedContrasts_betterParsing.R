#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#

#remotes::install_github('wolski/prolfquapp', dependencies = TRUE)
#remotes::install_github('fgcz/prolfqua', dependencies = TRUE)

rm(list=ls())

# libs
library(tidyverse)
library(prolfqua)
library(prolfquapp)
library(prophosqua)
library(dplyr)
library(stringr)
library(openxlsx)
library(seqinr)

# integration of phospho-DIA
message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")


# parameters and thresholds
# variables
WUID <- "CustomPhosphoAnalysis"
fgczProject <- "p29033"
OIDfgcz <- "p29033"
fracti <- "Integration"
datetoday <- format(Sys.Date(), "%Y%m%d")
descri <- WUID # "OneCondition"
(resDir <- paste0(fgczProject,"_",datetoday,"_",fracti,"_", descri))

path <- "."


# find path manually
# read back in results
(totXlsx  <-  paste0("DEA_p29033_20250304_TotalProteome_betterParsed/DEA_20250304_PIp29033_ObetterParsed_WUTotalProteome_robscale/Results_WU_TotalProteome/DE_WUTotalProteome.xlsx"))
(phosphoXlsx <- paste0("DEA_20250304_p29033_enriched_betterParsed/Results_WU_enriched/DE_WUenriched.xlsx"))


totRes <- readxl::read_xlsx(path = totXlsx, sheet = "diff_exp_analysis")
totRes$protein_length |> is.na() |> mean() # check how many proteins have no length ->  only rev proteins should have not length
totRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id) # these are the protein IDs without length! # should be only rev_ or proteins that are not found in this fasta

# if smaller case rev_ are used and not fixed in yaml file before running DEA this here helps to still get the integration done
# totRes$protein_Id <- sapply(strsplit((totRes$protein_Id), split = "\\|"), function(x)x[2])

rev_pattern = "^rev_"
totRes <- totRes |> dplyr::filter(!grepl("FGCZCont", protein_Id))
totRes <- totRes |> dplyr::filter(!grepl("contam_", protein_Id))
totRes <- totRes |> dplyr::filter(!grepl(rev_pattern, protein_Id)) # here we take out all revs!


phosphoRes <- readxl::read_xlsx(path = phosphoXlsx, sheet = "diff_exp_analysis")
phosphoRes$protein_length |> is.na() |> mean()
phosphoRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id) # these are the protein IDs without length! # should be only rev_

# drop contaminants and rev sequences
phosphoRes <- phosphoRes |> dplyr::filter(!grepl("FGCZCont", protein_Id))
phosphoRes <- phosphoRes |> dplyr::filter(!grepl("cont_", protein_Id))
phosphoRes <- phosphoRes |> dplyr::filter(!grepl(rev_pattern, protein_Id))

# add here more phospho related things
# add things to phospho site
phosphoRes$BGScollapseKey <- sapply(strsplit(phosphoRes$site, "~"), `[`, 2)
phosphoRes$Multiplicity <- sapply(strsplit(phosphoRes$BGScollapseKey, "_"), `[`, 3)
phosphoRes$SiteNlocation <- sapply(strsplit(phosphoRes$BGScollapseKey, "_"), `[`, 2)
# needed downstream
phosphoRes$posInProtein <- as.numeric(gsub(x = phosphoRes$SiteNlocation, pattern = "[S|T|Y](\\d+)", replacement = "\\1"))

# some columns are still missing for downstream plotting
phosphoRes$startModSite <- phosphoRes$posInProtein - 1
phosphoRes$endModSite <- phosphoRes$posInProtein + 1
# get the first modified AA
phosphoRes$modAA <- gsub(x = phosphoRes$SiteNlocation, pattern = "([S|T|Y])\\d+", replacement = "\\1")

# num phos, numLoc and AllLocalized
phosphoRes$AllLocalized <- TRUE

# where is the fasta file (regular uniprot fasta)
fasta_file = "../2024-10-16-reviewed-contam-UP000000589.fasta"

# source helper functions
source("/Users/jonasgrossmann/GitHub/prophosqua/R/FP_phosphoHelperFunctions_v3_202310.R")

# join total and phospho analysis
join_column <- c("fasta.id" = "fasta.id", "protein_Id","contrast", "description", "protein_length")

# do diff-diff and additional propagation of errors as suggested by MSstatsPTM
combined_test_diff <- test_diff(phosphoRes, totRes, join_column = join_column)

# write out html
resDir
dir.create(resDir)
drumm <- prolfquapp::make_DEA_config_R6(
  PROJECTID = fgczProject,
  ORDERID = OIDfgcz,
  WORKUNITID = descri)


# render html
rmarkdown::render("_Overview_Ubi_Integration_WEW.Rmd", params = list(data = combined_test_diff, grp = drumm, phosres = phosphoRes), output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE))
file.copy(from = "_Overview_Ubi_Integration_WEW.html", to = file.path(resDir, "Result_phosphoAndIntegration.html"))


# write to excel
excelResultList <- list()
excelResultList$combinedStats <- combined_test_diff
#excelResultList$seq_window <- seq_window
writexl::write_xlsx(excelResultList, path = file.path(resDir,"Result_phosphoAndTotalIntegration.xlsx"))
#


# N-to-C plotting
# Function to determine the significance for plotting NtoC
fdrThreshold = 0.01 # only proteins where in at least one contrast the protein is significant (fdrThreshold) are plotted

# build up candidate matrix for plotting
candidateMat <- combined_test_diff[!is.na(combined_test_diff$FDR.site), ]
cand <- candidateMat[candidateMat$FDR.site < fdrThreshold ,]
nrow(cand)

# proteins with sites regulated in any of the contrasts.
mySigProteinHits <- unique(cand$protein_Id)
length(mySigProteinHits)

# look at one of the contrasts
comboMat <- candidateMat
comboMat <- comboMat |> dplyr::filter(protein_Id %in% mySigProteinHits)


# how many true are ok -> position parsed properly -> this is irrelevant for site centric approach here
table(!is.na(comboMat$posInProtein))

# these columns are used later for N-to-C plots
comboMat_min <- dplyr::select(
  comboMat,
  c("protein_Id",
    "contrast",
    "protein_length",
    "site",
    "diff.protein",
    "diff.site",
    "FDR.site",
    "posInProtein",
    "startModSite",
    "endModSite",
    "AllLocalized",
    "modAA",
    model_site = "modelName.site"
  ))

# nest all sites behind
comboMat_min <- comboMat_min |> dplyr::group_by(protein_Id, contrast, protein_length) |> tidyr::nest()
comboMat_min$plot <- vector(mode = "list", nrow(comboMat_min))

# how many phospho proteins at the beginning?
length(unique(phosphoRes$protein_Id))


# check how many pages are plotted
length(unique(comboMat_min$protein_Id))

# fill all slots with plots
for (i in 1:nrow(comboMat_min)) {
  print(i)
  comboMat_min$plot[[i]] <- prophosqua::N_to_C_plot(comboMat_min$data[[i]],
                                                    comboMat_min$protein_Id[[i]],
                                                    comboMat_min$protein_length[[i]],
                                                    comboMat_min$contrast[[i]])
}

# Do the plotting only for each contrast individually
for (j in 1:length(unique(comboMat_min$contrast))) {
  print(j)
  oneC_comboMat <- comboMat_min[comboMat_min$contrast == unique(comboMat_min$contrast)[j],]
  pdfFN <- paste0("SignificantProteins_",unique(comboMat_min$contrast)[j],"_NtoCplots.pdf")
  pdf(file.path(resDir, pdfFN))
  for (i in 1:nrow(oneC_comboMat)) {
    print(oneC_comboMat$plot[[i]])
    grid::grid.newpage()
    table <- oneC_comboMat$data[[i]]
    table <- table |> select(-all_of(c( "startModSite", "endModSite", "AllLocalized")))
    table_grob <- gridExtra::tableGrob(table, theme = gridExtra::ttheme_default(base_size=6))
    grid::grid.draw(table_grob)
  }
  dev.off()
}

#
#rFN <- paste0(",")
save.image(paste0(resDir, ".RData"))



