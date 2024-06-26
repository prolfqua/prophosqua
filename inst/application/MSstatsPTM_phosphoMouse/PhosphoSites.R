#
# PTM (ST and Y enriched)
#
library(tidyverse)

# params ideally taken from yaml
fgczProject <- "pIDxx"
OIDfgcz <- "oxxx"
descri <- "MSstatsPTM_phosphoMouse"
fracti <- "PhosphoEnriched"
WUID <- "WUID"


# also copy the phospho specific Rmd files from prophosqua
prophosqua::copy_phosphoDEA_FragPipe_TMT()
#
path = "."


# work on GRP for having better folder name
(fN <- paste0(fgczProject,"_", descri, "_", WUID,"_",fracti))
GRP2_phos <- prolfquapp::make_DEA_config_R6(ZIPDIR = fN,PROJECTID = fgczProject,
                                            ORDERID = OIDfgcz)
GRP2_phos$processing_options$aggregate <- "topN"


# read in ST and Y results
fromTSV_ST <- read_tsv("pST_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")
fromTSV_Y <- read_tsv("pY_ADD-MASSIVE-REANALYSIS-e9ceea92-display_quant_results-main.tsv")

stypsm <- data.frame(rbind(fromTSV_ST, fromTSV_Y))
#head(stypsm)
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
# stypsm$idx <- 1:nrow(stypsm)
stypsmX$PeptideProphet.Probability <- 1
stypsmX$qValue <- 0.001

# working on site and protein
stypsmX$GeneName <- sapply(strsplit((stypsmX$ProteinName), split = "\\|"), function(x)x[1])
stypsmX$ProtNsite <- sapply(strsplit((stypsmX$ProteinName), split = "\\|"), function(x)x[2])
unique(stypsmX$ProtNsite) # can have multiple sites: "Q9QZQ1_S193_S196_Y203"
stypsmX$Acc <- sapply(strsplit((stypsmX$ProtNsite), split = "_"), function(x)x[1])

annotable <- read_tsv(file = "annotation.tsv")
annot <- prolfquapp::read_annotation(annotable)


# Setup configuration
atable_phos <- annot$atable
atable_phos$ident_Score = "PeptideProphet.Probability"
atable_phos$ident_qValue = "qValue"
atable_phos$fileName <- "raw"
atable_phos$hierarchy[["protein_Id"]] <- c("Acc")
atable_phos$hierarchy[["peptide_Id"]] <- c("ProtNsite")
atable_phos$hierarchy[["precursor"]] <- c("pepSeqNcharge")
atable_phos$hierarchyDepth <- 2 # no roll-up to protein
atable_phos$set_response("Intensity")


# Preprocess data - aggregate proteins.
config_phos <- prolfqua::AnalysisConfiguration$new(atable_phos)
phospsm2 <- dplyr::inner_join(annot$annot, stypsmX, multiple = "all")
adata_phos <- prolfqua::setup_analysis(phospsm2, config_phos)


# get lfq object
lfqdata_phos <- prolfqua::LFQData$new(adata_phos, config_phos)
lfqdata_phos$remove_small_intensities(threshold = 1)

# here we need more parsing w/ site!
pa_phos <- data.frame(protein_Id = unique(lfqdata_phos$data$protein_Id))
#pa_phos <- tidyr::separate(pa_phos, protein_Id , c(NA, "IDcolumn"), sep = "_",remove = FALSE) # done before
pa_phos$description <- "description needed"
pa_phos$IDcolumn <- pa_phos$protein_Id


protAnnot_phos <- prolfquapp::ProteinAnnotation$new(lfqdata_phos, pa_phos, cleaned_ids = "IDcolumn")
protAnnot_phos$row_annot

lfqdata <- prolfquapp::aggregate_data(lfqdata_phos, agg_method = GRP2_phos$processing_options$aggregate)
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2_phos, protAnnot_phos, Contrasts = annot$contrasts)

#debug(write_DEA_all)
prolfquapp::write_DEA_all(grp2 = grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")













##


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


