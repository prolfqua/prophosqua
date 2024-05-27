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


fromTSV <- read_tsv("global_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
#fromTSV <- read_tsv("global_ProteoSAFe-ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results/global_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
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
lfqdata$hierarchy_counts()


pa <- data.frame(protein_Id = unique(lfqdata$data$protein_Id))
pa <- tidyr::separate(pa, protein_Id , c(NA, "IDcolumn"), sep = "\\|",remove = FALSE)
pa$description <- "description needed"

# problem
#protAnnot <- prolfquapp::ProteinAnnotation$new(lfqdata, pa, cleaned_ids = "IDcolumn")
protAnnot <- ProteinAnnotation$new(lfqdata, pa, cleaned_ids =  "IDcolumn")

protAnnot$row_annot

lfqdata$config$table$hkeysDepth()
GRP2$processing_options$aggregate
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$processing_options$aggregate)

#logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
lfqdata$factors()
lfqdata$to_wide()
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, protAnnot, Contrasts = annot$contrasts)

#debug(write_DEA_all)
prolfquapp::write_DEA_all(grp2 = grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")


