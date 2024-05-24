#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# global data from here:
# https://www.dropbox.com/scl/fi/xm5779z007qmhjycd2kxj/global_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv?rlkey=1djm8yajrxizj87qbnzv9e9b8&dl=0

# here we want to try prophosqua with published data MSstatsPTM paper
# https://www.sciencedirect.com/science/article/pii/S1535947622002857#tbl1
# https://github.com/devonjkohler/MSstatsPTM_simulations/tree/main/data

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
prophosqua::copy_phosphoDEA_FragPipe_TMT()
#
path = "."

# work on GRP for having better folder name
GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = "fN",PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)


fromTSV <- read_tsv("global_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
colnames(fromTSV) <- c("id", "ProteinName", "PeptideSequence", "z", "pepSeqNcharge", "plex", "TechRep", "Run", "channel", "Condition", "BiolRep", "Intensity")

# idea: filter here for only 2 or 4 conditions to slim it down

psm <- as.data.frame(fromTSV)

# get an overview
(annotable <- psm |> select(BiolRep, Condition, Run, TechRep, channel, plex) |> distinct())

annotable$BiolRep |> table() |> table()


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
annot <- prolfquapp::read_annotation(annotable)



# add missing stuff
psm <- psm |> rename(raw = BiolRep)
psm[["Condition"]] <- NULL
psm[["Run"]] <- NULL
psm[["TechRep"]] <- NULL
psm[["plex"]] <- NULL
psm[["channel"]] <- NULL
psm$id <- NULL
psm$z <- NULL
head(psm)
psmX <- psm |> group_by(ProteinName,PeptideSequence,pepSeqNcharge,raw) |> summarize(n = n(), Intensity = sum(Intensity)) |> ungroup()
psmX$n |> table() ## WTF?


psmX$PeptideProphet.Probability <- 1
psmX$qValue <- 0.001


# Setup configuration
atable <- annot$atable
atable$ident_Score = "PeptideProphet.Probability"
atable$ident_qValue = "qValue"
atable$fileName <- "raw"
atable$hierarchy[["protein_Id"]] <- c("ProteinName")
atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
atable$hierarchy[["precursor"]] <- c("pepSeqNcharge")
atable$hierarchyDepth <- 1

atable$set_response("Intensity")


# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
psm2 <- dplyr::inner_join(annot$annot, psmX, multiple = "all")
colnames(psm2)
adata <- prolfqua::setup_analysis(psm2, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities(threshold = 1)


pa <- data.frame(protein_Id = unique(lfqdata$data$protein_Id))
pa <- tidyr::separate(pa, protein_Id , c(NA, "IDcolumn"), sep = "\\|",remove = FALSE)
pa$description <- "description needed"

protAnnot <- prolfquapp::ProteinAnnotation$new(lfqdata, pa, cleaned_ids = "IDcolumn")
protAnnot$row_annot

lfqdata$config$table$hkeysDepth()
GRP2$processing_options$aggregate
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$processing_options$aggregate)

#logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
lfqdata$factors()
lfqdata$to_wide()
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, protAnnot, Contrasts = annot$contrasts)

debug(write_DEA_all)
prolfquapp::write_DEA_all(grp2 = grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")


#
# PTM (ST and Y enriched)
#


# params ideally taken from yaml
fgczProject <- "pIDxx"
OIDfgcz <- "oxxx"
descri <- "MSstatsPTM_phosphoMouse"
fracti <- "PhosphoEnriched"
WUID <- "WUID"


# also copy the phospho specific Rmd files from prophosqua
# prophosqua::copy_phosphoDEA_FragPipe_TMT() # not yet working, package not built?
#
path = "."


# work on GRP for having better folder name
(fN <- paste0(fgczProject,"_", descri, "_", WUID,"_",fracti))
GRP2_phos <- prolfquapp::make_DEA_config_R6(ZIPDIR = fN,PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)
GRP2_phos$processing_options$aggregate



#fromTSV_ST <- read_tsv("phosphoST-ProteoSAFe-ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results/pST_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
#str(fromTSV_ST)

fromTSV_Y <- read_tsv("phosphoY_ProteoSAFe-ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results (1)/pY_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
str(fromTSV_Y)

#stypsm <- data.frame(rbind(fromTSV_ST, fromTSV_Y))
#head(stypsm)

stypsm <- fromTSV_Y

# make compatible first by summarizing plexes
stypsm <- stypsm |> select(ProteinName, PeptideSequence, Charge, PSM, BioReplicate, Condition, Intensity) |> group_by(ProteinName, PeptideSequence, PSM, BioReplicate, Condition,)  |> summarise(SumIntensity=sum(Intensity, na.rm = TRUE))
head(stypsm)

# add missing and modify
# stypsm$idx <- 1:nrow(stypsm)
stypsm$PeptideProphet.Probability <- 1
stypsm$qValue <- 0.001
stypsm$GroupingVariable <- stypsm$Condition
stypsm$CONTROL <- "T"
stypsm$CONTROL[stypsm$Condition == "WT_Uninfect"] <- "C"
#stypsm$channel <- paste(stypsm$Channel, stypsm$Mixture, sep = "_")
stypsm$sampleName <- stypsm$BioReplicate
stypsm$ProtNsite <- sapply(strsplit((stypsm$ProteinName), split = "\\|"), function(x)x[2])
stypsm$Acc <- sapply(strsplit((stypsm$ProtNsite), split = "_"), function(x)x[1]) # parse clean accession
stypsm$desc <- "desc"
stypsm$ProtNpepSeq <- paste(stypsm$ProtNsite, stypsm$PeptideSequence, sep="~")
stypsm$desc <- "myDescription"
stypsm$nrPeps <- 1
stypsm$BioReplicate <- NULL
stypsm$Channel <- stypsm$sampleName

head(stypsm)
table(stypsm$sampleName)
table(stypsm$Channel)

# Setup configuration
atable_phos <- prolfqua::AnalysisTableAnnotation$new()
atable_phos$ident_Score = "PeptideProphet.Probability"
atable_phos$ident_qValue = "qValue"
atable_phos$fileName = "sampleName"
atable_phos$hierarchy[["protein_Id"]] <- c("Acc")
#atable_phos$hierarchy[["peptide_Id"]] <- c("ProtNpepSeq")
atable_phos$hierarchy[["site"]] <- c("ProtNsite")
atable_phos$hierarchy[["Spectrum"]] <- c("PSM")

atable_phos$hierarchyDepth <- 2
atable_phos$set_response("SumIntensity")
atable_phos$factors

#
#debug(dataset_set_factors_deprecated)
tmp_phos <- prolfquapp::dataset_set_factors_deprecated(atable_phos, stypsm)
atable_phos <- tmp_phos$atable
atable_phos$factors
atable_phos$hierarchy
stypsm <- tmp_phos$msdata
str(stypsm)

# Preprocess data - aggregate proteins.
config_phos <- prolfqua::AnalysisConfiguration$new(atable_phos)
colnames(stypsm)
adata_phos <- prolfqua::setup_analysis(stypsm, config_phos)
colnames(adata_phos)

lfqdata_phos <- prolfqua::LFQData$new(adata_phos, config_phos)
lfqdata_phos$hierarchy_counts()
lfqdata_phos$remove_small_intensities(threshold = 1)


#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata_phos$config$table$hkeysDepth()
lfqdata_phos$config$table$hierarchyDepth <- 2
lfqdata_phos$to_wide()
lfqdata_phos$hierarchy_counts()

lfqdata_phos$response()
GRP2_phos$pop$aggregate
lfqdata_phos <- prolfquapp::aggregate_data(lfqdata_phos, agg_method = GRP2_phos$processing_options$aggregate)
#LFQDataAggregator$debug("medpolish")
#agg <- lfqdata_phos$get_Aggregator()
#agg$medpolish()

#

lfqdata_phos$to_wide() # strange.. many NAs
#logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
#lfqdata$factors()

logger::log_info("END OF DATA TRANSFORMATION.")

# Build protein annot w/ new function
colnames(stypsm)

protAnnot <- build_protein_annot(
  lfqdata_phos,
  stypsm,
  c("protein_Id" = "ProteinName"),
  cleaned_protein_id = "Acc",
  protein_description = "desc",
  nr_children = "nrPeps",
  more_columns = NULL)

lfqdata_phos$to_wide()

# we get errors here
grp <- prolfquapp::generate_DEA_reports(lfqdata_phos, GRP2, protAnnot) # this is taking quite a while


