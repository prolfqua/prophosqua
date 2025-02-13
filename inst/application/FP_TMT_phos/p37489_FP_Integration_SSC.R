#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

#remotes::install_github('wolski/prolfquapp', dependencies = TRUE)
#remotes::install_github('fgcz/prolfqua', dependencies = TRUE)

# libs
library(tidyverse)
library(prolfqua)
library(prolfquapp)
library(prophosqua)
library(dplyr)
library(stringr)
library(openxlsx)
library(seqinr)

# integration of 2-plex phospho-TMT data
message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")


# parameters and thresholds
# variables
WUID <- "specifiedContrasts"
fgczProject <- "p37489"
OIDfgcz <- "o37489"
fracti <- "Integration"
descri <- WUID # "OneCondition"
(resDir <- paste0(fgczProject, "_",OIDfgcz,"_",fracti,"_", descri))

path <- "."


# find path manually
# read back in results
(totXlsx  <-  paste0("DEA_20250212_PIp37489_Oo37489_WUTotalProteome_none/Results_WU_TotalProteome/DE_DEA_20250212_PIp37489_Oo37489_WUTotalProteome_none_WUTotalProteome.xlsx"))
(phosXlsx <- paste0("p37489_PhosphoEnriched_PhosphoEnriched/DEA_20250212_PIp37489_Oo37489_WUPhosphoEnriched_none/Results_WU_PhosphoEnriched/DE_WUPhosphoEnriched.xlsx"))


totRes <- readxl::read_xlsx(path = totXlsx, sheet = "diff_exp_analysis")
totRes$protein_length |> is.na() |> mean() # check how many proteins have no length @ WeW how can this happen? these are all the fw proteins.. only rev proteins are properly parsed?
totRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id) # these are the protein IDs without length! # should be only rev_

# if smaller case rev_ are used and not fixed in yaml file before running DEA this here helps to still get the integration done
# totRes$protein_Id <- sapply(strsplit((totRes$protein_Id), split = "\\|"), function(x)x[2])

rev_pattern = "^rev_"
totRes <- totRes |> dplyr::filter(!grepl("FGCZCont", protein_Id))
totRes <- totRes |> dplyr::filter(!grepl("contam_", protein_Id))
totRes <- totRes |> dplyr::filter(!grepl(rev_pattern, protein_Id)) # here we take out all revs!


phosRes <- readxl::read_xlsx(path = phosXlsx, sheet = "diff_exp_analysis")
phosRes$protein_length |> is.na() |> mean()
phosRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id) # these are the protein IDs without length! # should be only rev_

# drop contaminants and rev sequences
phosRes <- phosRes |> dplyr::filter(!grepl("FGCZCont", protein_Id))
phosRes <- phosRes |> dplyr::filter(!grepl("cont_", protein_Id))
phosRes <- phosRes |> dplyr::filter(!grepl(rev_pattern, protein_Id))

# clean fasta.ids from leading "rev_" pattern -> comes  from but in script.sh-step
phosRes$fasta.id <-  gsub(x = phosRes$fasta.id, pattern = rev_pattern, replacement = "")

# work out site in protein AA on phosRes
colnames(phosRes)
# phospho specific parsing
phosRes$originalSite <- phosRes$site
phosRes$peptideSequence <- sapply(strsplit((phosRes$site), split = "~"), function(x)x[2])
phosRes$site <- sapply(strsplit((phosRes$site), split = "~"), function(x)x[1])
phosRes$AccFromSite <- sapply(strsplit((phosRes$site), split = "_"), function(x)x[1])
phosRes$startModSite <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[2]))
phosRes$endModSite <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[3]))
phosRes$NumPhos <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[4]))
phosRes$LocalizedNumPhos <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[5]))
phosRes$PhosSites <- sapply(strsplit((phosRes$site), split = "_"), function(x)x[6])
phosRes$SinglePhos_bool <- phosRes$NumPhos == 1
phosRes$AllLocalized <- phosRes$NumPhos == phosRes$LocalizedNumPhos
phosRes$SinglePhosLocalized_bool <- phosRes$NumPhos == 1 & phosRes$LocalizedNumPhos == 1

# parse residue
phosRes$AA <- gsub(x = phosRes$PhosSites, pattern = "\\d+", replacement = "")
# parse position in protein
phosRes$posInProtein <- as.numeric(gsub(x = phosRes$PhosSites, pattern = "[STY]", replacement = ""))



#
phosRes <- phosRes |> dplyr::mutate(
  AllLocalized = (NumPhos == LocalizedNumPhos)
)
phosRes$AllLocalized |> mean()

fasta_file = "../proteome/2024-11-28-decoys-contam-UP000006548.fasta"


# source helper functions
source("/Users/jonasgrossmann/GitHub/prophosqua/R/FP_phosphoHelperFunctions_v3_202310.R")
seq_window <- get_sequence_windows(phosRes, fasta_file, rev_pattern)

join_column <- c("fasta.id" = "protein_Id", "contrast", "description", "protein_length", "nr_tryptic_peptides")

combined_test_diff <- test_diff(phosRes, totRes, join_column = c("fasta.id" = "protein_Id", "contrast","description", "protein_length"))
combined_test_diff$AllLocalized |> mean(na.rm = TRUE)

phosRes <- phosRes |> dplyr::mutate(AA = substr(PhosSites, 1, 1))

# write out html
resDir
dir.create(resDir)
drumm <- prolfquapp::make_DEA_config_R6(
  PROJECTID = fgczProject,
  ORDERID = OIDfgcz,
  WORKUNITID = descri)
rmarkdown::render("_Overview_PhosphoAndIntegration_WEW.Rmd", params = list(data = combined_test_diff, grp = drumm, phosres = phosRes), output_format = bookdown::html_document2(toc = TRUE, toc_float = TRUE))
file.copy(from = "_Overview_PhosphoAndIntegration_WEW.html", to = file.path(resDir, "PhosphoAndIntegration.html"))

# write to excel
excelResultList <- list()
excelResultList$combinedStats <- combined_test_diff
excelResultList$seq_window <- seq_window
writexl::write_xlsx(excelResultList, path = file.path(resDir,"PhosphoAndIntegration.xlsx"))
#


# Function to determine the significance for plotting NtoC
fdrThreshold = 0.00000001


candidateMat <- combined_test_diff[!is.na(combined_test_diff$FDR.site), ]
cand <- candidateMat[candidateMat$FDR.site < fdrThreshold ,]
nrow(cand)

# proteins with sites regulated in any of the contrasts.
mySigProteinHits <- unique(cand$protein_Id)
length(mySigProteinHits)

# look at one of the contrasts
comboMat <- candidateMat
comboMat <- comboMat |> dplyr::filter(protein_Id %in% mySigProteinHits)
candidateMat$AllLocalized |> mean()

# add data frame with sequence windows
comboMat <- dplyr::inner_join(comboMat, seq_window)


# select columns necessary for plotting
comboMat_min <- dplyr::select(
  comboMat,
  c("protein_Id",
    "contrast",
    "protein_length",
    "site",
    "diff.protein",
    "diff.site",
    "FDR.site",
    "modAA" = "AA",
    "posInProtein",
    "startModSite",
    "endModSite",
    "AllLocalized",
    model_site = "modelName.site"
  ))


comboMat_min <- comboMat_min |> dplyr::group_by(protein_Id, contrast, protein_length) |> tidyr::nest()
comboMat_min$plot <- vector(mode = "list", nrow(comboMat_min))

# how many phospho proteins at the beginning?
length(unique(phosRes$protein_Id))


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


# all in one pdf
# pdf(file.path(resDir, "NtoCplots2.pdf"))
# for (i in 1:nrow(comboMat_min)) {
#   print(comboMat_min$plot[[i]])
#   grid::grid.newpage()
#   table <- comboMat_min$data[[i]]
#   table <- table |> select(-all_of(c( "startModSite", "endModSite", "AllLocalized")))
#   table_grob <- gridExtra::tableGrob(table, theme = gridExtra::ttheme_default(base_size=6))
#   grid::grid.draw(table_grob)
# }
# dev.off()




#
#rFN <- paste0(",")
save.image(paste0(resDir, ".RData"))



