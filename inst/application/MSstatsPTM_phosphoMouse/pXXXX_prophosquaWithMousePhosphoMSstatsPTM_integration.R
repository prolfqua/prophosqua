message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")

# annotation and descrison based on:
# https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details


# for developping:
# setwd("~/FGCZ/Interdisziplinar/pXXXX_Prophosqua/mouse_phospho")

# libs
library(prolfqua)
library(prolfquapp)
library(prophosqua)
library(dplyr)
library(stringr)
library(openxlsx)
library(seqinr)


# parameters and thresholds
# params ideally taken from yaml
fgczProject <- "pXXX"
fracti <- "Integration"
descri <- "_phosphoMouse"
(resDir <- paste0(fgczProject, "_",fracti, descri))



# read back in results
(totXlsx  <-  paste0("FIRST_TOTAL_fN_OI_oxxx_WU_/Results_DEA_WU/DE_Groups_vs_Controls_WU.xlsx"))
(phosXlsx  <-  paste0("FIRST_PHOSPHO_pIDxx_MSstatsPTM_phosphoMouse_WUID_PhosphoEnriched_PI_pIDxx_OI_oxxx_WU__none/Results_DEA_WU/DE_Groups_vs_Controls_WU.xlsx"))
totRes <- read.xlsx(xlsxFile = totXlsx, sheet = "diff_exp_analysis")
phosRes <- read.xlsx(xlsxFile = phosXlsx, sheet = "diff_exp_analysis")

# some little adaptations is necessairy
colnames(phosRes)[4] <- "site"


# now here we read in decoy and only return fw entries
myFasta <- readDecoyFastaNreturnFwOnly("fgcz_10090_1spg_d_20230405.fasta")

# Filter out FGCZContaminants -> problematic at sites
phosRes <- phosRes |> filter(!grepl("FGCZCont", protein_Id))


# # chatgpt:
# # Example strings
# strings <- c("A2AUK8_S437_S441", "A2AR50_S311", "A2AUK8_S437_S441")
# # Use sub() to extract the part after the first underscore
# parsed_strings <- sub("^[^_]*_", "", strings)
# # Print the result
# print(parsed_strings)
#
# # or
# # Example strings
# strings <- c("A2AUK8_S437_S441", "A2AR50_S311", "A2AUK8_S437_S441")
# # Use strsplit() to split the strings at the first underscore and then sapply() to extract the second part
# parsed_strings <- sapply(strsplit(strings, "_", fixed = TRUE), function(x) paste(x[-1], collapse = "_"))
# # Print the result
# print(parsed_strings)
#


# phospho specific parsing some adaptions necessairy for phosphoMouse
phosRes$originalSite <- phosRes$site
phosRes$peptideSequence <- "PEPTIDEK"
phosRes$site <- sapply(strsplit(phosRes$originalSite, "_", fixed = TRUE), function(x) paste(x[-1], collapse = "_"))
phosRes$AccFromSite <- sapply(strsplit((phosRes$originalSite), split = "_"), function(x)x[1])
phosRes$FirstSite <- sapply(strsplit((phosRes$site), split = "_"), function(x)x[1])
phosRes$startModSite <- 100
phosRes$endModSite <- 110
phosRes$NumPhos <- as.numeric(str_count(string = phosRes$site, pattern = "_")+1)
phosRes$LocalizedNumPhos <- phosRes$NumPhos
phosRes$PhosSites <- phosRes$FirstSite
phosRes$SinglePhos_bool <- phosRes$NumPhos == 1
phosRes$AllLocalized <- phosRes$NumPhos == phosRes$LocalizedNumPhos
phosRes$SinglePhosLocalized_bool <- phosRes$NumPhos == 1 & phosRes$LocalizedNumPhos == 1

# parse residue
phosRes$AA <- gsub(x = phosRes$PhosSites, pattern = "\\d+", replacement = "")
# parse position in protein
phosRes$posInProtein <- as.numeric(gsub(x = phosRes$PhosSites, pattern = "[STY]", replacement = ""))

# add info if found in more than one unique peptide for TMT-integrator this is not necessairy
phosRes$SiteFoundInManyPeptides <- str_count(string = phosRes$peptide, pattern = ";") + 1

# localized or not
table(phosRes$SinglePhosLocalized_bool)
table(phosRes$AA[phosRes$SinglePhosLocalized_bool])

colnames(phosRes)
# # find mapping for sequence window
# phosRes |> select(protein_Id, peptideSequence, posInProtein) |> distinct() |> dim_desc()
# phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct() |> dim_desc()
#
# uniqueProtPepSeq <- phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct()
#
# # parse sequence window from fasta givn position
# for (i in 1:nrow(uniqueProtPepSeq)) {
#   #print(i)
#   uniqueProtPepSeq$SequenceWindows[i] <- getSequenceWindowForLogo(poi = uniqueProtPepSeq$protein_Id[i], soi = uniqueProtPepSeq$posInProtein[i], fastaObject = myFasta)
# }
# tail(uniqueProtPepSeq)

# join sequence windows back!
# phosRes <- left_join(x = phosRes, y = uniqueProtPepSeq)
phosRes$SequenceWindows <- "ARBITRARARYWINDOWTHATDOESNOTMAKESENSE"


# join sheets
# get rid  of decoys
phosRes  <- phosRes |> filter(!grepl("REV_", protein_Id))
totRes  <- totRes |> filter(!grepl("REV_", protein_Id))

head(totRes)
head(phosRes)
# parse  middle part from totRes$proteinID -> sp|A0A0D9S1R0|APOE_CHLSB
combo <- left_join(x = phosRes, y = totRes, join_by("IDcolumn" == "IDcolumn", "contrast" == "contrast"))


# MS stats like adjustment for protein change
comboWithAdj <- prophosqua::doMSstatsLikeSiteNormalizationUsingProteinStatsOnComboObject(combo)
colnames(comboWithAdj)
head(comboWithAdj)

# some
comboWithAdj$protein_Id <- comboWithAdj$protein_Id.x
comboWithAdj$protein_Id.x <- NULL
# for RMD report
GRP2 <- prolfquapp::make_DEA_config(PROJECTID = fgczProject, ORDERID = fgczProject, WORKUNITID = "WUxxx")

# Idea for reporting and downstream processing

# render integration html
(resultPath <- resDir)
(htmlFN <- paste0("Integration",descri))
prolfquapp::render_DEA(GRP2 = comboWithAdj, outpath = resultPath, htmlname = htmlFN, word = FALSE, markdown = "_Overview_PhosphoAndIntegration.Rmd")

# separate for each constrast
(NumCoi <- length(unique(comboWithAdj$contrast)))

# for how many do we have global protein
colnames(comboWithAdj)
totProtFromSites <- comboWithAdj |> select(protein_Id, protein_Id.y) |> distinct() |> nrow()
totProtWGlobalProt <-  comboWithAdj |> select(protein_Id, protein_Id.y) |> distinct() |> filter(!is.na(protein_Id.y)) |> nrow()
round(totProtWGlobalProt/totProtFromSites, 2)

# N-to-C plotting per contrast
fdrThreshold <- 0.1 # for plotting N-to-C

for (i in 1:length(unique(comboWithAdj$contrast))) {
  comboMat <- comboWithAdj[comboWithAdj$contrast == unique(comboWithAdj$contrast)[i],]
  # Towards N-to-C plots
  sum(comboWithAdj$FDR.x < fdrThreshold)
  candidateMat <- comboWithAdj[comboWithAdj$FDR.x < fdrThreshold,]
  generateNtoCProteinPDFsWithPhosphoPeptides_FragPipeTMT(globalNphosphoCombinedNResultMatrix = comboMat, candidateMatrix = candidateMat, expName = unique(comboWithAdj$contrast)[i])
}

# move N-to-C plots to Integration result Dir
(myNtoCplotsFiles <- list.files(pattern = "SignificantProtein_.*\\.pdf"))
file.copy(myNtoCplotsFiles, to = resultPath)
file.remove(myNtoCplotsFiles)


# write to excel
#prolfquapp::write_DEA(GRP2 = comboWithAdj, outpath = "IntegrationResults", xlsxname = "IntegrationPhosphoCentric", write = TRUE)
excelResultList <- list()
excelResultList$combinedStats <- comboWithAdj
writexl::write_xlsx(excelResultList, path = paste0(resultPath, "/",htmlFN,".xlsx"))

