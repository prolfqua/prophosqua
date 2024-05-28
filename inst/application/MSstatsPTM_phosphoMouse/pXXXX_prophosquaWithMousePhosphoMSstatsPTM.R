################################################################################
#
#
#      psm -> Total part
#
#
################################################################################
library(tidyverse)
library(prolfqua)
library(prolfquapp)
library(prophosqua)

# params ideally taken from yaml
fgczProject <- "pIDxx"
OIDfgcz <- "oxxx"
descri <- "MSstatsPTM_phosphoMouse"
fracti <- "TotalProteome"
WUID <- "WUID"

# v3
prolfquapp::copy_DEA_DIANN()

#
path = "."

# work on GRP for having better folder name
GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = "fN",PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)


# fromTSV <- read_tsv("global_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
fromTSV <- read_tsv("global_ProteoSAFe-ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results/global_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
colnames(fromTSV) <- c("id", "ProteinName", "PeptideSequence", "z", "pepSeqNcharge", "plex", "TechRep", "Run", "channel", "Condition", "BiolRep", "Intensity")

# idea: filter here for only 2 or 4 conditions to slim it down
psm <- as.data.frame(fromTSV)

# get an overview
(annotable <- psm |> select(BiolRep, Condition, Run, TechRep, channel, plex) |> distinct())

annotable$BiolRep |> table() |> table()

# prepare annot table
annotable <- annotable |> rename(group = Condition)
annotable$CONTROL <- "T"
annotable$CONTROL[annotable$group == "WT_Uninfect"] <- "C"

annotable$plex <- NULL
annotable$channel <- NULL
annotable$TechRep <- NULL
annotable$Batch <- paste0("B",annotable$Run)
annotable$Run <- NULL
annotable$raw <- annotable$BiolRep
annotable <- annotable |> rename(Name = BiolRep)
write_tsv(annotable, file = "annotation.tsv")

# final annotation table
annot <- prolfquapp::read_annotation(annotable)



# dropping bs and add missing stuff
psm <- psm |> rename(raw = BiolRep)
psm[["Condition"]] <- NULL
psm[["Run"]] <- NULL
psm[["TechRep"]] <- NULL
psm[["plex"]] <- NULL
psm[["channel"]] <- NULL
psm$id <- NULL
psm$z <- NULL
psmX <- psm |> group_by(ProteinName,PeptideSequence,pepSeqNcharge,raw) |> summarize(n = n(), Intensity = sum(Intensity)) |> ungroup()
psmX$n |> table() ## WTF?

# adding things
psmX$PeptideProphet.Probability <- 1
psmX$qValue <- 0.001


# Setup configuration
atable <- annot$atable
atable$ident_Score = "PeptideProphet.Probability"
atable$ident_qValue = "qValue"
atable$fileName <- "raw"
atable$hierarchy[["protein_Id"]] <- c("ProteinName")
#atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
atable$hierarchy[["site"]] <- c("PeptideSequence")
atable$hierarchy[["precursor"]] <- c("pepSeqNcharge")
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
pa <- tidyr::separate(pa, protein_Id , c(NA, "IDcolumn"), sep = "\\|",remove = FALSE)
pa$description <- "description needed"

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



#
# PTM (ST and Y enriched)
#


# params ideally taken from yaml
#fgczProject <- "pIDxx"
#OIDfgcz <- "oxxx"
#descri <- "MSstatsPTM_phosphoMouse"
fracti <- "PhosphoEnriched"
#WUID <- "WUID"


# also copy the phospho specific Rmd files from prophosqua
prophosqua::copy_phosphoDEA_FragPipe_TMT()
#
path = "."


# work on GRP for having better folder name
(fN <- paste0(fgczProject,"_", descri, "_", WUID,"_",fracti))
GRP2_phos <- prolfquapp::make_DEA_config_R6(ZIPDIR = fN,PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)
# adjust roll-up
GRP2_phos$processing_options$aggregate <- "topN"


# read in ST and Y results
fromTSV_ST <- read_tsv("phosphoST-ProteoSAFe-ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results/pST_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
fromTSV_Y <- read_tsv("phosphoY_ProteoSAFe-ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results (1)/pY_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")

stypsm <- data.frame(rbind(fromTSV_ST, fromTSV_Y))
colnames(stypsm) <- c("rowid", "ProteinName", "PeptideSequence", "z", "pepSeqNcharge", "plex", "TechRep", "Run", "channel", "Condition", "raw", "Intensity", "id")

# dropping bs and add missing stuff
stypsm[["Condition"]] <- NULL
stypsm[["Run"]] <- NULL
stypsm[["TechRep"]] <- NULL
stypsm[["plex"]] <- NULL
stypsm[["channel"]] <- NULL
stypsm$id <- NULL
stypsm$rowid <- NULL
stypsm$z <- NULL
head(stypsm)
stypsmX <- stypsm |> group_by(ProteinName,PeptideSequence,pepSeqNcharge,raw) |> summarize(n = n(), Intensity = sum(Intensity)) |> ungroup()
stypsmX$n |> table() ## some are measured twice therefore intensity aggregated w/ sum

# add missing and modify
stypsmX$PeptideProphet.Probability <- 1
stypsmX$qValue <- 0.001

# working on site and protein
stypsmX$GeneName <- sapply(strsplit((stypsmX$ProteinName), split = "\\|"), function(x)x[1])
stypsmX$ProtNsite <- sapply(strsplit((stypsmX$ProteinName), split = "\\|"), function(x)x[2])
#unique(stypsmX$ProtNsite) # can have multiple sites: "Q9QZQ1_S193_S196_Y203"
stypsmX$Acc <- sapply(strsplit((stypsmX$ProtNsite), split = "_"), function(x)x[1])

# Setup configuration
atable_phos <- annot$atable # take same annot table from total
atable_phos$ident_Score = "PeptideProphet.Probability"
atable_phos$ident_qValue = "qValue"
atable_phos$fileName <- "raw"
atable_phos$hierarchy[["protein_Id"]] <- c("Acc") # change here if different roll-up id should be taken.. e.g. GeneName
atable_phos$hierarchy[["site"]] <- c("ProtNsite")
atable_phos$hierarchy[["precursor"]] <- c("pepSeqNcharge")
atable_phos$hierarchyDepth <- 2 # no roll-up to protein but to protNsite

atable_phos$set_response("Intensity")


# Preprocess data - aggregate proteins.
config_phos <- prolfqua::AnalysisConfiguration$new(atable_phos)
phospsm2 <- dplyr::inner_join(annot$annot, stypsmX, multiple = "all")
adata_phos <- prolfqua::setup_analysis(phospsm2, config_phos)


# get lfq object
lfqdata_phos <- prolfqua::LFQData$new(adata_phos, config_phos)
lfqdata_phos$hierarchy_counts()
lfqdata_phos$remove_small_intensities(threshold = 1)
lfqdata_phos$hierarchy_counts()

# here we need more parsing w/ site!
pa_phos <- data.frame(protein_Id = unique(lfqdata_phos$data$protein_Id))
pa_phos$description <- "description needed"
pa_phos$IDcolumn <- pa_phos$protein_Id
protAnnot_phos <- prolfquapp::ProteinAnnotation$new(lfqdata_phos, pa_phos, cleaned_ids = "IDcolumn")
protAnnot_phos$row_annot

# aggregating
lfqdata_phos <- prolfquapp::aggregate_data(lfqdata_phos, agg_method = GRP2_phos$processing_options$aggregate)

# modelling
grp_phos <- prolfquapp::generate_DEA_reports2(lfqdata_phos, GRP2_phos, protAnnot_phos, Contrasts = annot$contrasts)

# all fine? some checks
#grp_phos$RES$lfqData$to_wide()
#grp_phos$RES$contrastsData_signif
#myResPlotter <- grp_phos$RES$contrMerged$get_Plotter()
#myResPlotter$volcano()

logger::log_info("DONE WITH DEA REPORTS")

# result dir
GRP2_phos$zipdir
dir.create(GRP2_phos$zipdir)

# need helper functions to properly write reports not on protein but peptide level
# source("FP_phosphoHelperFunctions_v3_202310.R")

# writing reports
prolfquapp::write_DEA_all(grp2 = grp_phos, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")
prolfquapp::write_DEA_all(grp2 = grp_phos, boxplot = FALSE, markdown = "_DiffExpQC_Phospho_V2.Rmd")

