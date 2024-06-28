message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")

# annotation and descrison based on:
# https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details


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
fgczProject <- "pIDxx"
fracti <- "Integration"
descri <- "_vsCond1"
(resDir <- paste0(fgczProject, "_",fracti, descri))



# read back in results
paste0(GRP2$get_result_dir())
(totXlsx  <-  paste0(GRP2$get_result_dir(),"/DE__WUWUID.xlsx"))
(phosXlsx <- paste0(GRP2_phos$get_result_dir() , "/DE_Groups_vs_Controls_WUWUID.xlsx"))
totRes <- read.xlsx(xlsxFile = totXlsx, sheet = "diff_exp_analysis")
phosRes <- read.xlsx(xlsxFile = phosXlsx, sheet = "diff_exp_analysis")

# now here we read in decoy and only return fw entries
myFasta <- readDecoyFastaNreturnFwOnly(files$fasta)

# Filter out FGCZContaminants -> problematic at sites
phosRes <- phosRes |> filter(!grepl("FGCZCont", protein_Id))

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

# add info if found in more than one unique peptide for TMT-integrator this is not necessairy
phosRes$SiteFoundInManyPeptides <- str_count(string = phosRes$peptide, pattern = ";") + 1

# localized or not
table(phosRes$SinglePhosLocalized_bool)
table(phosRes$AA[phosRes$SinglePhosLocalized_bool])

# find mapping for sequence window
phosRes |> select(protein_Id, peptideSequence, posInProtein) |> distinct() |> dim_desc()
phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct() |> dim_desc()

uniqueProtPepSeq <- phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct()

# parse sequence window from fasta givn position
for (i in 1:nrow(uniqueProtPepSeq)) {
  #print(i)
  uniqueProtPepSeq$SequenceWindows[i] <- getSequenceWindowForLogo(poi = uniqueProtPepSeq$protein_Id[i], soi = uniqueProtPepSeq$posInProtein[i], fastaObject = myFasta)
}
tail(uniqueProtPepSeq)

# join sequence windows back!
phosRes <- left_join(x = phosRes, y = uniqueProtPepSeq)

# join sheets
# get rid  of decoys
phosRes  <- phosRes |> filter(!grepl("REV_", protein_Id))
totRes  <- totRes |> filter(!grepl("REV_", protein_Id))

# parse  middle part from totRes$proteinID -> sp|A0A0D9S1R0|APOE_CHLSB
totRes$protAcc <- sapply(strsplit((totRes$protein_Id), split = "\\|"), function(x)x[2])
combo <- left_join(x = phosRes, y = totRes, join_by("protein_Id" == "protAcc", "contrast" == "contrast"))


# MS stats like adjustment for protein change
comboWithAdj <- doMSstatsLikeSiteNormalizationUsingProteinStatsOnComboObject(combo)
colnames(comboWithAdj)
head(comboWithAdj)

# for RMD report
GRP2 <- prolfquapp::make_DEA_config(PROJECTID = fgczProject, ORDERID = fgczProject, WORKUNITID = "WUxxx")

# Idea for reporting and downstream processing
# number of proteins
uniqueProtPepSeqSite <- comboWithAdj |> select(protein_Id, peptideSequence, site, originalSite) |> distinct()
length(unique(uniqueProtPepSeqSite$protein_Id))
myPhosphoProts <- rle(sort(uniqueProtPepSeqSite$protein_Id))

table(myPhosphoProts$lengths)
hist(myPhosphoProts$lengths, breaks=1000)

myPhosphoProts$values[which(myPhosphoProts$lengths == max(myPhosphoProts$lengths))]
uniqueProtPepSeqSite$originalSite[which(uniqueProtPepSeqSite$protein_Id == myPhosphoProts$values[which(myPhosphoProts$lengths == max(myPhosphoProts$lengths))])]

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
#
save.image(paste0(htmlFN, ".RData"))
