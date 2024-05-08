message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")

# annotation and comparison based on:
# https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details


# libs
library(prolfqua)
library(prolfquapp)
library(readr)
library(dplyr)
library(stringr)
library(openxlsx)
library(seqinr)

# helperfunctions
source("FP_phosphoHelperFunctions_v3_202310.R")


# parameters and thresholds
# params ideally taken from yaml
fgczProject <- "pXXXX"
descri <- "TMTphospho_integration"
compari <- "_SimulatedOne"
WUID <- "WUxx"

# read back in results
#paste0(resDir,"Results_DEA_WU",WUID,"/DE_Groups_vs_Controls.xlsx" )
totRes <- read.xlsx(xlsxFile = "pXXXX_SimulationTMTphospho__TotalProteome/Results_DEA_WUWUxx/DE_Groups_vs_Controls_WUWUxx.xlsx", sheet = "diff_exp_analysis")
phosRes <- read.xlsx(xlsxFile = "pXXXX_SimulationTMTphospho__PhosphoEnriched/DE_Groups_vs_Controls/DE_Groups_vs_Controls.xlsx", sheet = "diff_exp_analysis")

(resDir <- paste0(fgczProject, "_",descri, compari))

# there is no fasta for simlulated data
# myFasta <- seqinr::read.fasta("../fgcz_9606_reviewed_cnl_20230330.fasta", seqtype = "AA", as.string = TRUE)


# work on phosRes
head(phosRes$site)

# Filter out FGCZContaminants -> problematic at sites
phosRes <- phosRes |> filter(!grepl("FGCZCont", protein_Id))


#
phosRes$originalSite <- phosRes$site
#phosRes$peptideSequence <- sapply(strsplit((phosRes$site), split = "~"), function(x)x[2])
phosRes$site <- sapply(strsplit((phosRes$site), split = "~"), function(x)x[1])
phosRes$AccFromSite <- sapply(strsplit((phosRes$site), split = "_"), function(x)x[1])
phosRes$startModSite <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[2]))
phosRes$endModSite <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[3]))
#phosRes$NumPhos <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[4]))
#phosRes$LocalizedNumPhos <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[5]))
phosRes$NumPhos <- 1
phosRes$LocalizedNumPhos <- 1
phosRes$PhosSites <- gsub(x = phosRes$site, pattern = "_",replacement = "")
phosRes$SinglePhos_bool <- phosRes$NumPhos == 1
phosRes$AllLocalized <- phosRes$NumPhos == phosRes$LocalizedNumPhos
phosRes$SinglePhosLocalized_bool <- phosRes$NumPhos == 1 & phosRes$LocalizedNumPhos == 1
phosRes$peptideSequence <- "PEPTIDEK"

# work on phospho to
phosRes$AA <- gsub(x = phosRes$PhosSites, pattern = "\\d+", replacement = "")
# parse position in protein
phosRes$posInProtein <- as.numeric(gsub(x = phosRes$PhosSites, pattern = "[STY]", replacement = ""))

# add info if found in more than one unique peptide
phosRes$SiteFoundInManyPeptides <- 1

table(phosRes$AA[phosRes$SinglePhosLocalized_bool])
head(phosRes$posInProtein[phosRes$SinglePhosLocalized_bool])

# find mapping for sequence window
phosRes |> select(protein_Id, peptideSequence, posInProtein) |> distinct() |> dim_desc()
phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct() |> dim_desc()

uniqueProtPepSeq <- phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct()

# we do not have a fasta
# for (i in 1:nrow(uniqueProtPepSeq)) {
#   uniqueProtPepSeq$SequenceWindows[i] <- getSequenceWindowForLogo(poi = uniqueProtPepSeq$protein_Id[i], soi = uniqueProtPepSeq$posInProtein[i], fastaObject = myFasta)
# }

uniqueProtPepSeq$SequenceWindows <- "MYPEPTIDEKWINDWWFAKE"
tail(uniqueProtPepSeq)

# join sequence windows back!
phosRes <- left_join(x = phosRes, y = uniqueProtPepSeq)

# join sheets
# get rid  of decoys
phosRes  <- phosRes |> filter(!grepl("REV_", protein_Id))
totRes  <- totRes |> filter(!grepl("REV_", protein_Id))

# parse  middle part from totRes$proteinID -> sp|A0A0D9S1R0|APOE_CHLSB
#totRes$protAcc <- sapply(strsplit((totRes$protein_Id), split = "\\|"), function(x)x[2])
totRes$protAcc <- sapply(strsplit((totRes$protein_Id), split = "\\|"), function(x)x[1])
phosRes$pNr <- sapply(strsplit((phosRes$protein_Id), split = "_"), function(x)x[2])
phosRes$pNam <- sapply(strsplit((phosRes$protein_Id), split = "_"), function(x)x[1])
phosRes$cleanPhosProtein <- paste(phosRes$pNam, phosRes$pNr, sep = "_")
# combine
combo <- left_join(x = phosRes, y = totRes, join_by("cleanPhosProtein" == "protAcc", "contrast" == "contrast"))


# MS stats like adjustment for protein change
comboWithAdj <- doMSstatsLikeSiteNormalizationUsingProteinStatsOnComboObject(combo)
colnames(comboWithAdj)
head(comboWithAdj)

# save.image("Integration_o33038_allIn.RData")
# load("Integration_vs5hpiNoGSK.RData")

# for RMD report
GRP2 <- prolfquapp::make_DEA_config(PROJECTID = fgczProject, ORDERID = fgczProject, WORKUNITID = "WUxxx")

# Idea for reporting and downstream processing

# Number crunching
# number of sites split STY
# all of them
#barplot(table(comboWithAdj$NumPhos[comboWithAdj$AllLocalized]), main = "Number of identified Phosphorylations per Peptide")

# only single phos-localized -> in RMD
# table(comboWithAdj$AA[comboWithAdj$SinglePhosLocalized_bool])
#barplot(table(comboWithAdj$AA[comboWithAdj$SinglePhosLocalized_bool]))


# Abundances and STY -> in RMD
# p <- ggplot2::ggplot(na.omit(comboWithAdj[comboWithAdj$SinglePhosLocalized_bool,]), ggplot2::aes(x=AA, y=avgAbd.x)) +
#   ggplot2::geom_boxplot()
# p


# separate for each constrast
(NumCoi <- length(unique(comboWithAdj$contrast)))


# N-to-C plotting per contrast
# fdrThreshold <- 0.1 # for plotting N-to-C
#
# for (i in 1:length(unique(comboWithAdj$contrast))) {
#   comboMat <- comboWithAdj[comboWithAdj$contrast == unique(comboWithAdj$contrast)[i],]
#   # Towards N-to-C plots
#   sum(comboWithAdj$FDR.x < fdrThreshold)
#   candidateMat <- comboWithAdj[comboWithAdj$FDR.x < fdrThreshold,]
#   generateNtoCProteinPDFsWithPhosphoPeptides_FragPipeTMT(globalNphosphoCombinedNResultMatrix = comboMat, candidateMatrix = candidateMat, expName = unique(comboWithAdj$contrast)[i])
# }


# render integration html
(resultPath <- resDir)
(htmlFN <- paste0("Integration",compari))
comboWithAdj$protein_Id <- comboWithAdj$protein_Id.x
prolfquapp::render_DEA(GRP2 = comboWithAdj, outpath = resultPath, htmlname = htmlFN, word = FALSE, markdown = "_Overview_PhosphoAndIntegration.Rmd")

# write to excel
#prolfquapp::write_DEA(GRP2 = comboWithAdj, outpath = "IntegrationResults", xlsxname = "IntegrationPhosphoCentric", write = TRUE)
excelResultList <- list()
excelResultList$combinedStats <- comboWithAdj
writexl::write_xlsx(excelResultList, path = paste0(resultPath, "/",htmlFN,".xlsx"))
