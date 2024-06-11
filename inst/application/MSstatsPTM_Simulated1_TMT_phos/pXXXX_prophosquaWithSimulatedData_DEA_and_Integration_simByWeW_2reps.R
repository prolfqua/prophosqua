#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# here we want to try prophosqua with simlated data from the MSstatsPTM paper
# https://www.sciencedirect.com/science/article/pii/S1535947622002857#tbl1
# https://github.com/devonjkohler/MSstatsPTM_simulations/tree/main/data

library(tidyverse)
library(prolfqua)
library(prolfquapp)
library(readr)
library(openxlsx)


#load(file = "simulation1_data.rda") # downloaded from github page -> this one is flawed .. 10 features also in the PTM data
load(file="simulation1_data_newByWeW.rda") # this one is fixed")

simulation1_data[[1]]$PTM$Run |> table() # 1 is the one with 4 Runs -> 2 reps
simulation1_data[[1]]$PTM$Condition |> table() # 1 has 2 conditions
#
#
# # check which idx to use from simulated_data
# simulation1_data[[20]]$PTM$Run |> table() # 12 files
# simulation1_data[[20]]$PTM$Condition |> table() # 4 conditions
#
# simulation1_data[[3]]$PTM$Run |> table()  # triplicates
# simulation1_data[[3]]$PTM$Condition |> table()
#
# simulation1_data[[5]]$PTM$Run |> table() # quintuplicates
# simulation1_data[[5]]$PTM$Condition |> table()

idxOfInterest <- 1

# we can generate them ourselves with the provided code at least!

# params ideally taken from yaml
fgczProject <- "pXXXX"
OIDfgcz <- "oYYYY"
descri <- "simTwoGrps2Reps_idOne"
fracti <- "TotalProteome"
WUID <- "WUxx"

# work on GRP for having better folder name
(fN <- paste0(descri,fracti))
GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = fN,PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)
psm <- data.frame(simulation1_data[[idxOfInterest]]$PROTEIN)

# get an overview
(annotable <- psm |> select(BioReplicate, Condition, Run) |> distinct())
psm |> group_by(Condition, BioReplicate) |> summarise(n = n())


# prepare annot table
annotable <- annotable |> rename(group = Condition)
annotable$CONTROL <- "T"
annotable$CONTROL[annotable$group == "G_1"] <- "C"

annotable$raw <- annotable$BioReplicate
(annotable <- annotable |> rename(Name = BioReplicate))
write_tsv(annotable, file = "annotation_forSimulatedData_2reps_rr.tsv")

# final annotation table
annot <- prolfquapp::read_annotation(annotable)

# dropping bs and add missing stuff
colnames(psm)
psm$PrecursorCharge <- NULL
psm$FragmentIon <- NULL
psm$ProductCharge <- NULL

psm <- psm |> rename(raw = BioReplicate)
psm[["Condition"]] <- NULL
psm[["Run"]] <- NULL

# reshaping
psmX <- psm |> group_by(ProteinName,PeptideSequence,raw) |> summarize(n = n(), Intensity = sum(Intensity)) |> ungroup()
psmX$n |> table()


# adding things
psmX$PeptideProphet.Probability <- 1
psmX$qValue <- 0.001
psmX$oldID <- psmX$ProteinName
#psmX$ProteinName <- gsub(x = psmX$ProteinName, pattern = "\\|.*", replacement = "") # we need to keep all separate

# Setup configuration
atable <- annot$atable
atable$ident_Score = "PeptideProphet.Probability"
atable$ident_qValue = "qValue"
atable$fileName <- "raw"
atable$hierarchy[["protein_Id"]] <- c("ProteinName")
atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
atable$hierarchyDepth <- 1
atable$set_response("Intensity")


# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
psm2 <- dplyr::inner_join(annot$annot, psmX, multiple = "all")
adata <- prolfqua::setup_analysis(psm2, config)


# lfq data object
lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities(threshold = 1)
lfqdata$hierarchy_counts()

# protein annotation
pa <- data.frame(protein_Id = unique(lfqdata$data$protein_Id))
pa$description <- "description needed"
pa$IDcolumn <- pa$protein_Id

protAnnot <- ProteinAnnotation$new(lfqdata, pa, cleaned_ids =  "IDcolumn")
protAnnot$row_annot

lfqdata$config$table$hkeysDepth()

GRP2$processing_options$aggregate
# aggregate from psm to protein level here
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$processing_options$aggregate)

# grp object
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, protAnnot, Contrasts = annot$contrasts)

copy_DEA_DIANN()
# write reports
dir.create(grp$zipdir)
prolfquapp::write_DEA_all(grp2 = grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")






# now we go for the phospho enriched data
#
#
#
#       PTM
#
#
#
#

# load(file = "simulation1_data.rda") # downloaded from github page
fracti <- "PhosphoEnriched"
#fgczProject <- "pXXXX"
#OIDfgcz <- "oYYYY"
#descri <- "Simulation_redone_ONETMTphospho_"
#WUID <- "WUxx"

## Sim -> from: https://github.com/devonjkohler/MSstatsPTM_simulations/blob/main/code/simulate_model_data.R
# 250 proteins where site is changed but not protein! no pipe in protein_id
# 250 proteins where site is changed and protein is changed in identical way -> no_change1 these should not be differentially expressed
# 500 proteins where nothing is changed -> no_change2
# sim <- PTMsimulateExperiment(
#   nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=250, nSite=1, nFeature=2, nFeature_prot = 10,
#   logAbundance=list(
#     PTM=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=.25),
#     PROTEIN=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], sPeak=0.25))
# )
# sim_no_change1 <- PTMsimulateExperiment(
#   nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=250, nSite=1, nFeature=2, nFeature_prot = 10,
#   logAbundance=list(
#     PTM=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=0.25),
#     PROTEIN=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=0.25))
# )
# sim_no_change2 <- PTMsimulateExperiment(
#   nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=500, nSite=1, nFeature=2, nFeature_prot = 10,
#   logAbundance=list(
#     PTM=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], sPeak=0.25),
#     PROTEIN=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], sPeak=0.25))
#


multiSite_long <- data.frame(simulation1_data[[idxOfInterest]]$PTM)

# work on GRP for having better folder name
GRP2_phos <- prolfquapp::make_DEA_config_R6(ZIPDIR = "fN",PROJECTID = fgczProject,
                                       ORDERID = fracti)

# get an overview
(annotable_phos <- multiSite_long |> select(BioReplicate, Condition, Run) |> distinct())
# 4 files, 2 conditions
multiSite_long |> group_by(Condition, BioReplicate) |> summarise(n = n())


# prepare annot table
annotable_phos <- annotable_phos |> rename(group = Condition)
annotable_phos$CONTROL <- "T"
annotable_phos$CONTROL[annotable_phos$group == "G_1"] <- "C"

annotable_phos$raw <- annotable_phos$BioReplicate
(annotable_phos <- annotable_phos |> rename(Name = BioReplicate))
write_tsv(annotable_phos, file = "annotation_PTM_forSimulatedData_2reps.tsv")

# final annotation table
annot_phos <- prolfquapp::read_annotation(annotable_phos)

# dropping bs and add missing stuff
colnames(multiSite_long)
multiSite_long$PrecursorCharge <- NULL
multiSite_long$FragmentIon <- NULL
multiSite_long$ProductCharge <- NULL

multiSite_long <- multiSite_long |> rename(raw = BioReplicate)
multiSite_long[["Condition"]] <- NULL
multiSite_long[["Run"]] <- NULL


head(unique(multiSite_long$ProteinName))
tail(unique(multiSite_long$ProteinName))

# adding things
multiSite_long$PeptideProphet.Probability <- 1
multiSite_long$qValue <- 0.001
multiSite_long$oldID <- multiSite_long$ProteinName
#multiSite_long$protNsite <- gsub(x = multiSite_long$ProteinName, pattern = "\\|.*", replacement = "")
multiSite_long$protNsite <- multiSite_long$ProteinName

# Setup configuration
atable_phos <- annot_phos$atable
atable_phos$ident_Score = "PeptideProphet.Probability"
atable_phos$ident_qValue = "qValue"
atable_phos$fileName <- "raw"
atable_phos$hierarchy[["protein_Id"]] <- c("ProteinName")
atable_phos$hierarchy[["site"]] <- c("ProteinName","PeptideSequence","protNsite")

atable_phos$hierarchyDepth <- 1
atable_phos$set_response("Intensity")


# Preprocess data - aggregate proteins.
config_phos <- prolfqua::AnalysisConfiguration$new(atable_phos)
psm2 <- dplyr::inner_join(annot_phos$annot, multiSite_long, multiple = "all")
adata_phos <- prolfqua::setup_analysis(psm2, config_phos)


# lfq data object
lfqdata_phos <- prolfqua::LFQData$new(adata_phos, config_phos)
lfqdata_phos$hierarchy_counts()
lfqdata_phos$data
lfqdata_phos$remove_small_intensities(threshold = 1)
lfqdata_phos$hierarchy_counts()

# protein annotation
pa_phos <- data.frame(protein_Id = unique(lfqdata_phos$data$protein_Id))
#pa <- data.frame(protein_Id = unique(gsub(x = psmX$ProteinName, pattern = "\\|.*", replacement = "")))
#pa <- tidyr::separate(pa, protein_Id , c(NA, "IDcolumn"), sep = "\\|",remove = FALSE)
#pa$IDcolumn <- pa$protein_Id
pa_phos$description <- "description needed"
pa_phos$IDcolumn <- pa_phos$protein_Id

protAnnot_phos <- ProteinAnnotation$new(lfqdata_phos, pa_phos, cleaned_ids =  "IDcolumn")
protAnnot_phos$row_annot

# roll up to site level
lfqdata_phos$config$table$hierarchyDepth <- 1
lfqdata_phos$config$table$hkeysDepth()

GRP2_phos$processing_options$aggregate <- "medpolish" # anyway we do not aggregate here
# aggregate from psm to protein level here
# now we also aggregate since it is the same feature in the same protein
lfqdata_phos <- prolfquapp::aggregate_data(lfqdata_phos, agg_method = GRP2_phos$processing_options$aggregate) # LHS error

# grp object
grp_phos <- prolfquapp::generate_DEA_reports2(lfqdata_phos, GRP2_phos, protAnnot_phos, Contrasts = annot_phos$contrasts)

logger::log_info("DONE WITH DEA REPORTS")

# result dir
(GRP2_phos$zipdir <- paste0(descri,fracti,"_", Sys.Date()))
dir.create(GRP2_phos$zipdir)

# need helper functions to properly write reports not on protein but peptide level
copy_DEA_DIANN()
library(prophosqua)
copy_phosphoDEA_FragPipe_TMT()

# writing reports
GRP2_phos$RES$lfqData$data$site <- GRP2_phos$RES$lfqData$data$protein_Id
colnames(GRP2_phos$RES$contrastsData)
GRP2_phos$RES$contrastsData$site <- GRP2_phos$RES$contrastsData$protein_Id
GRP2_phos$RES$contrastsData_signif$site <- GRP2_phos$RES$contrastsData_signif$protein_Id

GRP2 <- GRP2_phos
prolfquapp::write_DEA_all(grp2 = grp_phos, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")
prolfquapp::write_DEA_all(grp2 = grp_phos, boxplot = FALSE, markdown = "_DiffExpQC_Phospho_V2.Rmd")

# integration starting here
#rm(list=ls())
# parameters and thresholds
# params ideally taken from yaml
fgczProject <- "pXXXX"
descri <- "integration"
compari <- "_2repsReRedone"
WUID <- "WUxx"

# read back in results
totRes <- read.xlsx(xlsxFile = "simTwoGrps2Reps_idOneTotalProteome_PI_pXXXX_OI_oYYYY_WU__none/Results_DEA_WU/DE_Groups_vs_Controls_WU.xlsx", sheet = "diff_exp_analysis")
phosRes <- read.xlsx(xlsxFile = "simTwoGrps2Reps_idOnePhosphoEnriched_2024-06-11//Results_DEA_WU/DE_Groups_vs_Controls_WU.xlsx", sheet = "diff_exp_analysis")

# site missing since we rolled up to protein
phosRes$site <- phosRes$protein_Id

(resDir <- paste0(fgczProject, "_",descri, compari))

# some parsing .. here often not useful but necessairy to not run into errors
phosRes$originalSite <- phosRes$site
phosRes$site <- sapply(strsplit((phosRes$site), split = "~"), function(x)x[1])
phosRes$AccFromSite <- sapply(strsplit((phosRes$site), split = "_"), function(x)x[1])
phosRes$startModSite <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[2]))
phosRes$endModSite <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[3]))
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

# # add info if found in more than one unique peptide
# phosRes$SiteFoundInManyPeptides <- 1
#
# table(phosRes$AA[phosRes$SinglePhosLocalized_bool])
#
# uniqueProtPepSeq <- phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct()
#
# # we do not have a fasta
# # for (i in 1:nrow(uniqueProtPepSeq)) {
# #   uniqueProtPepSeq$SequenceWindows[i] <- getSequenceWindowForLogo(poi = uniqueProtPepSeq$protein_Id[i], soi = uniqueProtPepSeq$posInProtein[i], fastaObject = myFasta)
# # }
#
# uniqueProtPepSeq$SequenceWindows <- "MYPEPTIDEKWINDWWFAKE"
# tail(uniqueProtPepSeq)
#
# # join sequence windows back!
# phosRes <- left_join(x = phosRes, y = uniqueProtPepSeq)


# this is important for the integration of the simulated data
# parse  phosRes from "Protein_128|NoChange1_S_1" back to totRes$protein_Id[2] -> "Protein_1|NoChange1"
phosRes$fProt <- sapply(strsplit((phosRes$protein_Id), split = "\\|"), function(x)x[1])
# cut away after _.* -> "Protein_1"
# q: how can I do gsub with back reference in R? -> \\1
phosRes$fProt <- gsub(x = phosRes$fProt, pattern = "(Protein_\\d+).*", replacement = "\\1")
phosRes$Accpart <- sapply(strsplit((phosRes$protein_Id), split = "\\|"), function(x)x[2])
# replace NA with ""
phosRes$Accpart[is.na(phosRes$Accpart)] <- ""
phosRes$Accpart <- gsub(x = phosRes$Accpart, pattern = "(NoChange\\d+).*", replacement = "\\1")

phosRes$cleanPhosProtein <- paste(phosRes$fProt, phosRes$Accpart, sep = "|")

#problematic Protein_1| for no Accpart
phosRes$cleanPhosProtein3 <- gsub(x = phosRes$cleanPhosProtein, pattern = "\\|$", replacement = "")

# combine
combo <- left_join(x = phosRes, y = totRes, join_by("cleanPhosProtein3" == "protein_Id", "contrast" == "contrast"))


# MS stats like adjustment for protein change
comboWithAdj <- doMSstatsLikeSiteNormalizationUsingProteinStatsOnComboObject(combo)

# for RMD report
GRP2 <- prolfquapp::make_DEA_config(PROJECTID = fgczProject, ORDERID = fgczProject, WORKUNITID = "WUxxx")


# # render integration html
(resultPath <- resDir)
(htmlFN <- paste0("Integration",compari))
comboWithAdj$protein_Id <- comboWithAdj$IDcolumn.x
comboWithAdj$protein_Id.y <- comboWithAdj$IDcolumn.x

prolfquapp::render_DEA(GRP2 = comboWithAdj, outpath = resultPath, htmlname = htmlFN, word = FALSE, markdown = "_Overview_PhosphoAndIntegration.Rmd")

# write to excel
#rm(GRP2)
#prolfquapp::write_DEA(GRP2 = comboWithAdj, outpath = "IntegrationResults", xlsxname = "IntegrationPhosphoCentric", write = TRUE)

excelResultList <- list()
excelResultList$combinedStats <- comboWithAdj
writexl::write_xlsx(excelResultList, path = paste0(resultPath, "/",htmlFN,".xlsx"))


#
#
#  Evaluation of results to MSstatsPTM results
#
#

# functions
# some functions
# write function to get spezificity and sensitivity
# get_spezificity_sensitivity <- function(objectWithProteinIDs){
#   # get the true positives
#   TP <- sum(str_count(string = objectWithProteinIDs, pattern = "NoChange")==0)
#   # get the false positives
#   FP <- sum(str_count(string = objectWithProteinIDs, pattern = "NoChange")>0)
#   # get the false negatives
#   FN <- TP - 250
#   # get the true negatives
#   TN <- 750 - FP
#   # length
#   len <- length(objectWithProteinIDs)
#   # get the spezificity
#   spezificity <- TN / (TN + FP)
#   # get the sensitivity
#   sensitivity <- TP / (TP + FN)
#   # get precision
#   precision <- TP / (TP + FP)
#   # get recall
#   recall <- TP / (TP + FN)
#   return(c(spezificity, sensitivity, precision, recall, len))
# }

# write function to get spezificity and sensitivity
get_empirical_FDR <- function(objectWithProteinIDs){
  # get the true positives
  TP <- sum(str_count(string = objectWithProteinIDs, pattern = "NoChange")==0)
  # get the false positives
  FP <- sum(str_count(string = objectWithProteinIDs, pattern = "NoChange")>0)
  # length
  len <- length(objectWithProteinIDs)
  # get eFDR
  eFDR <- FP / (TP + FP)
  return(c(eFDR, TP, FP, len))
}




# read in msstatsPTM results
load("adj_limma_models_sim1.rda")
res_MSstatsPTM <- adj_limma_sim1[[idxOfInterest]]
head(res_MSstatsPTM)
res_MSstatsPTM$adj.P.Val <- p.adjust(res_MSstatsPTM$pvalue, method = "BH")


sigThreshold <- 0.05
# get the significant results
sig_prophosqua <- comboWithAdj[comboWithAdj$MSstatsPTMadj_FDR < sigThreshold,]
sig_MSstatsPTM <- res_MSstatsPTM[res_MSstatsPTM$adj.P.Val < sigThreshold,]

# get the empirical FDR
round(get_empirical_FDR(objectWithProteinIDs = sig_MSstatsPTM$PTM),2)
round(get_empirical_FDR(objectWithProteinIDs = sig_prophosqua$IDcolumn.x),2)

# get the overlap
overlap <- intersect(sig_prophosqua$IDcolumn.x, sig_MSstatsPTM$PTM)
length(overlap)




