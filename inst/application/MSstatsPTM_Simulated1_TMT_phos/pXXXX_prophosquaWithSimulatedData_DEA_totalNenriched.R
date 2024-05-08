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


load(file = "simulation1_data.rda") # downloaded from github page
# we can generate them ourselves with the provided code at least!

# params ideally taken from yaml
fgczProject <- "pXXXX"
OIDfgcz <- "oYYYY"
descri <- "SimulationTMTphospho_"
fracti <- "TotalProteome"
WUID <- "WUxx"

# look into structure
# 1000 proteins with 10 peptides each are simluated in 4 files
# 2 conditions G_1 and G_2 with 2 reps each ->  no missing data
simulation1_data[[1]]$PROTEIN



#simulate annotation file
str(simulation1_data[[1]]$PROTEIN)
simulation1_data[[1]]$PROTEIN |> select(Condition, Run, feature, BioReplicate) |> distinct()
simulation1_data[[1]]$PROTEIN |> select(Condition, Run, feature, BioReplicate, PeptideSequence) |> distinct()
simulation1_data[[1]]$PROTEIN |> select(Condition, Run, BioReplicate, PeptideSequence) |> distinct()

# go for psm
psm <- data.frame(simulation1_data[[1]]$PROTEIN)

# add missing things
psm$"PeptideProphet.Probability" <- 1
psm$qValue <- 0.001
psm$rawfile <- psm$BioReplicate
psm$GroupingVar <- psm$Condition
psm$BioReplicate <- NULL

# C or T
psm$CONTROL <- NA
psm$CONTROL[psm$Condition == "G_1"] <- "C"
psm$CONTROL[psm$Condition == "G_2"] <- "T"

# now we prefer it to have it directly in the workingDir
(resDir <- paste(fgczProject, descri, fracti, sep="_"))
# GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = resDir, Normalization = "robscale", PROJECTID = fgczProject, ORDERID = OIDfgcz, WORKUNITID = WUID)
GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = resDir, Normalization = "none", PROJECTID = fgczProject, ORDERID = OIDfgcz, WORKUNITID = WUID) # should be available in new version

# Setup configuration
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$ident_Score = "PeptideProphet.Probability"
atable$ident_qValue = "qValue"
atable$fileName = "rawfile"
atable$hierarchy[["protein_Id"]] <- c("ProteinName")
atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
#atable$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide","Assigned.Modifications")
#atable$hierarchy[["Spectrum"]] <- c("Spectrum")
atable$set_response("Intensity")

#
#debug(dataset_set_factors_deprecated)
tmp <- prolfquapp::dataset_set_factors_deprecated(atable, psm)
atable <- tmp$atable
atable$factors
psm <- tmp$msdata
psm$desc <- "myDescription"
psm$nrPeps <- 1

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
colnames(psm)
adata <- prolfqua::setup_analysis(psm, config)
colnames(adata)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities(threshold = 1)


#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata$config$table$hkeysDepth()
lfqdata$config$table$hierarchyDepth <- 1
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)

#logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
#lfqdata$factors()

logger::log_info("END OF DATA TRANSFORMATION.")

# Build protein annot w/ new function
protAnnot <- build_protein_annot(
  lfqdata,
  psm,
  c("protein_Id" = "ProteinName"),
  cleaned_protein_id = "ProteinName",
  protein_description = "desc",
  nr_children = "nrPeps",
  more_columns = NULL)



grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, protAnnot) # this is taking quite a while

prolfquapp::copy_DEA_DIANN()
dir.create(GRP2$zipdir)

for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp2 = grp[[i]], name = names(grp)[i], ZIPDIR = GRP2$zipdir, boxplot = FALSE)
}


# PTM
# load(file = "simulation1_data.rda") # downloaded from github page
fracti <- "PhosphoEnriched"
#fgczProject <- "pXXXX"
#OIDfgcz <- "oYYYY"
#descri <- "SimulationTMTphospho_"
#WUID <- "WUxx"

multiSite_long <- data.frame(simulation1_data[[1]]$PTM)

# zip specify here already
(resDir <- paste(fgczProject, descri, fracti, sep="_"))

GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = resDir, Normalization = "none", PROJECTID = fgczProject, ORDERID = fgczProject, WORKUNITID = "WUxxx_Phospho")
# GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2) # we do not have proper annot file here

# parse annotation from input directly
# add missing things
multiSite_long$"PeptideProphet.Probability" <- 1
multiSite_long$qValue <- 0.001
multiSite_long$rawfile <- multiSite_long$BioReplicate
multiSite_long$GroupingVar <- multiSite_long$Condition
multiSite_long$BioReplicate <- NULL

# C or T
multiSite_long$CONTROL <- NA
multiSite_long$CONTROL[multiSite_long$Condition == "G_1"] <- "C"
multiSite_long$CONTROL[multiSite_long$Condition == "G_2"] <- "T"

# check
multiSite_long |> select(rawfile, GroupingVar) |> distinct() |> group_by(GroupingVar) |> summarise(n())

# Setup configuration
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$ident_Score = "PeptideProphet.Probability"
atable$ident_qValue = "qValue"
atable$fileName = "rawfile"
atable$hierarchy[["protein_Id"]] <- c("ProteinName")
atable$hierarchy[["site"]] <- c("site", "PeptideSequence")
atable$set_response("Intensity")

atable$hierarchyDepth <- 2

tmp <- prolfquapp::dataset_set_factors_deprecated(atable, multiSite_long)
atable <- tmp$atable
atable$factors
multiSite_long <- tmp$msdata
multiSite_long$desc <- "myDescription"
multiSite_long$nrPeps <- 1

# Configuration
config <- prolfqua::AnalysisConfiguration$new(atable)

adata <- prolfqua::setup_analysis(multiSite_long, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities(threshold = 1)
lfqdata$hierarchy_counts()

lfqdata$config$table$hierarchyDepth <- 2
# some checks
lfqdata$to_wide()
lfqdata$response()
lfqdata$hierarchy_counts()


#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata$config$table$hierarchy_keys()

# Build protein annot w/ new function
protAnnot <- build_protein_annot(
  lfqdata,
  multiSite_long,
  c("protein_Id" = "ProteinName"),
  cleaned_protein_id = "ProteinName",
  protein_description = "desc",
  nr_children = "nrPeps",
  more_columns = NULL)

# important for phospho
GRP2$pop$nr_peptdes <- 1
GRP2$pop$aggregate <- "none"

logger::log_info("GENERATING DEA REPORTS")
logger::log_info("starting modelling")

grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, protAnnot)

logger::log_info("DONE WITH DEA REPORTS")

source("FP_phosphoHelperFunctions_v3_202310.R")

# generate result folders
dir.create(GRP2$zipdir)

# write DEAs
for (i in seq_along(grp)) {
  #prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = FALSE)
  write_phosphoDEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = FALSE)
}



