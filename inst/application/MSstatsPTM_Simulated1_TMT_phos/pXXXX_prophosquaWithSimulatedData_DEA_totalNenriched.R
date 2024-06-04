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

# work on GRP for having better folder name
GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = "fN",PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)

# look into structure
# 1000 proteins with 10 peptides each are simluated in 4 files
# 2 conditions G_1 and G_2 with 2 reps each ->  no missing data
tail(simulation1_data[[1]]$PROTEIN, n=100)



#simulate annotation file
head(simulation1_data[[1]]$PROTEIN)
simulation1_data[[1]]$PROTEIN |> select(Condition, Run, feature, BioReplicate) |> distinct()
simulation1_data[[1]]$PROTEIN |> select(Condition, Run, feature, BioReplicate, PeptideSequence) |> distinct()
simulation1_data[[1]]$PROTEIN |> select(Condition, Run, BioReplicate, PeptideSequence) |> distinct()

unique(simulation1_data[[1]]$PROTEIN$ProteinName) |> length()
unique(simulation1_data[[1]]$PROTEIN$ProteinName)

psm <- data.frame(simulation1_data[[1]]$PROTEIN)
colnames(psm)

# get an overview
(annotable <- psm |> select(BioReplicate, Condition, Run) |> distinct())
# 4 files, 2 conditions
psm |> group_by(Condition, BioReplicate) |> summarise(n = n())


# prepare annot table
annotable <- annotable |> rename(group = Condition)
annotable$CONTROL <- "T"
annotable$CONTROL[annotable$group == "G_1"] <- "C"

annotable$raw <- annotable$BioReplicate
annotable <- annotable |> rename(Name = BioReplicate)
write_tsv(annotable, file = "annotation_forSimulatedData.tsv")

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
psmX$ProteinName <- gsub(x = psmX$ProteinName, pattern = "\\|.*", replacement = "")

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
#pa <- data.frame(protein_Id = unique(gsub(x = psmX$ProteinName, pattern = "\\|.*", replacement = "")))
#pa <- tidyr::separate(pa, protein_Id , c(NA, "IDcolumn"), sep = "\\|",remove = FALSE)
#pa$IDcolumn <- pa$protein_Id
pa$description <- "description needed"
pa$IDcolumn <- pa$protein_Id

protAnnot <- ProteinAnnotation$new(lfqdata, pa, cleaned_ids =  "IDcolumn")
protAnnot$row_annot

lfqdata$config$table$hkeysDepth()

GRP2$processing_options$aggregate
# aggregate from psm to protein level here
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$processing_options$aggregate)

# what is in the game
lfqdata$factors()
lfqdata$to_wide()

# grp object
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, protAnnot, Contrasts = annot$contrasts)

# write reports
prolfquapp::write_DEA_all(grp2 = grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")






# now we go for the phospho enriched data

# PTM
# load(file = "simulation1_data.rda") # downloaded from github page
fracti <- "PhosphoEnriched"
#fgczProject <- "pXXXX"
#OIDfgcz <- "oYYYY"
#descri <- "SimulationTMTphospho_"
#WUID <- "WUxx"

multiSite_long <- data.frame(simulation1_data[[1]]$PTM)

# work on GRP for having better folder name
GRP2_phos <- prolfquapp::make_DEA_config_R6(ZIPDIR = "fN",PROJECTID = fgczProject,
                                       ORDERID = fracti)

# look into structure
# 1000 proteins with 10 peptides each are simluated in 4 files
# 2 conditions G_1 and G_2 with 2 reps each ->  no missing data


multiSite_long |> select(Condition, Run, feature, BioReplicate) |> distinct()
multiSite_long |> select(Condition, Run, feature, BioReplicate, PeptideSequence) |> distinct()

# get an overview
(annotable_phos <- multiSite_long |> select(BioReplicate, Condition, Run) |> distinct())
# 4 files, 2 conditions
multiSite_long |> group_by(Condition, BioReplicate) |> summarise(n = n())


# prepare annot table
annotable_phos <- annotable_phos |> rename(group = Condition)
annotable_phos$CONTROL <- "T"
annotable_phos$CONTROL[annotable_phos$group == "G_1"] <- "C"

annotable_phos$raw <- annotable_phos$BioReplicate
annotable_phos <- annotable_phos |> rename(Name = BioReplicate)
write_tsv(annotable_phos, file = "annotation_PTM_forSimulatedData.tsv")

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

# reshaping
#(multiSite_longx <- multiSite_long |> group_by(ProteinName,PeptideSequence,raw) |> summarize(n = n(), Intensity = sum(Intensity)) |> ungroup())
#multiSite_longx$n |> table()

dim(multiSite_longx)
dim(multiSite_long)

# adding things
multiSite_long$PeptideProphet.Probability <- 1
multiSite_long$qValue <- 0.001
multiSite_long$oldID <- multiSite_long$ProteinName
multiSite_long$protNsite <- gsub(x = multiSite_long$ProteinName, pattern = "\\|.*", replacement = "")
multiSite_long$ProteinName <- gsub(x = multiSite_long$protNsite, pattern = "_S.*", replacement = "")

# Setup configuration
atable_phos <- annot_phos$atable
atable_phos$ident_Score = "PeptideProphet.Probability"
atable_phos$ident_qValue = "qValue"
atable_phos$fileName <- "raw"
atable_phos$hierarchy[["protein_Id"]] <- c("ProteinName")
atable_phos$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
atable_phos$hierarchy[["site"]] <- c("protNsite")


atable_phos$hierarchyDepth <- 2
atable_phos$set_response("Intensity")


# Preprocess data - aggregate proteins.
config_phos <- prolfqua::AnalysisConfiguration$new(atable_phos)
psm2 <- dplyr::inner_join(annot_phos$annot, multiSite_long, multiple = "all")
adata_phos <- prolfqua::setup_analysis(psm2, config_phos)


# lfq data object
lfqdata_phos <- prolfqua::LFQData$new(adata_phos, config_phos)
lfqdata_phos$hierarchy_counts()
lfqdata_phos$remove_small_intensities(threshold = 1)
lfqdata_phos$hierarchy_counts()

# protein annotation
pa_phos <- data.frame(protein_Id = unique(lfqdata_phos$data$protein_Id))
#pa <- data.frame(protein_Id = unique(gsub(x = psmX$ProteinName, pattern = "\\|.*", replacement = "")))
#pa <- tidyr::separate(pa, protein_Id , c(NA, "IDcolumn"), sep = "\\|",remove = FALSE)
#pa$IDcolumn <- pa$protein_Id
pa_phos$description <- "description needed"
pa_phos$IDcolumn <- pa$protein_Id

protAnnot <- ProteinAnnotation$new(lfqdata_phos, pa, cleaned_ids =  "IDcolumn")
protAnnot$row_annot

# roll up to site level
lfqdata_phos$config$table$hierarchyDepth <- 2
lfqdata_phos$config$table$hkeysDepth()



GRP2_phos$processing_options$aggregate <- "none"

# aggregate from psm to protein level here
lfqdata_phos <- prolfquapp::aggregate_data(lfqdata_phos, agg_method = GRP2$processing_options$aggregate) # LHS error

# what is in the game
lfqdata_phos$factors()
lfqdata_phos$to_wide()

# grp object
grp_phos <- prolfquapp::generate_DEA_reports2(lfqdata_phos, GRP2, protAnnot, Contrasts = annot_phos$contrasts)

# here we need more parsing w/ site!
# pa_phos <- data.frame(protein_Id = unique(lfqdata_phos$data$protein_Id))
# pa_phos$description <- "description needed"
# pa_phos$IDcolumn <- pa_phos$protein_Id
# protAnnot_phos <- prolfquapp::ProteinAnnotation$new(lfqdata_phos, pa_phos, cleaned_ids = "IDcolumn")
# protAnnot_phos$row_annot

# aggregating
lfqdata_phos <- prolfquapp::aggregate_data(lfqdata_phos, agg_method = GRP2_phos$processing_options$aggregate) # LHS error

# modelling
grp_phos <- prolfquapp::generate_DEA_reports2(lfqdata_phos, GRP2_phos, protAnnot, Contrasts = annot$contrasts)

# all fine? some checks
grp_phos$RES$lfqData$to_wide()
grp_phos$RES$contrastsData_signif
myResPlotter <- grp_phos$RES$contrMerged$get_Plotter()
myResPlotter$volcano()

logger::log_info("DONE WITH DEA REPORTS")

# result dir
GRP2_phos$zipdir
dir.create(GRP2_phos$zipdir)

# need helper functions to properly write reports not on protein but peptide level
# source("FP_phosphoHelperFunctions_v3_202310.R")

# writing reports
prolfquapp::write_DEA_all(grp2 = grp_phos, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")
prolfquapp::write_DEA_all(grp2 = grp_phos, boxplot = FALSE, markdown = "_DiffExpQC_Phospho_V2.Rmd")




# important for phospho
GRP2_phos$proup$nr_peptdes <- 1
GRP2$pop$aggregate <- "none"

logger::log_info("GENERATING DEA REPORTS")
logger::log_info("starting modelling")

source("FP_phosphoHelperFunctions_v3_202310.R")

# generate result folders
dir.create(GRP2$zipdir)

# write DEAs
for (i in seq_along(grp)) {
  #prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = FALSE)
  write_phosphoDEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = FALSE)
}

















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



