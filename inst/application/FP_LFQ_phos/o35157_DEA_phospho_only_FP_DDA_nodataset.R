dmessage("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquapp Version :", packageVersion("prolfquapp"), "\n")

library(prolfqua)
library(dplyr)
library(tidyr)
library(tibble)

# some fgcz specific settings
# params ideally taken from yaml
fgczProject <- "p35157"
OIDfgcz <- "o35157"
descri <- "phosphoOnly"
fracti <- "enriched"
WUID <- "WU304087"


# what part
whatpart <- "enriched"

# get LocalizationSiteProbThreshold
LocProbThresh <- 0.75

# do prolfqua on FragPipe phospho peptides!
###
source("FP_phosphoHelperFunctions_v3_202310.R")
proteinf <- "../o35157_FPtables_finalExp/combined_site_STY_79.9663.tsv"

#protein <- prolfqua::tidy_FragPipe_combined_protein(proteinf, spcnames = c(), intnames = c("Intensity"), maxlfqnames = c())
protein <- tidy_FragPipe_phospho_STYfile(styphosphopeptideTable = proteinf, locProbThreshold = LocProbThresh)



# build up annotatation from raw file names
unique(protein$raw.file) # 20 samples in total
table(gsub(x = unique(protein$raw.file) , pattern = "\\d$", replacement = ""))
table(gsub(x = gsub(x = unique(protein$raw.file) , pattern = "\\d$", replacement = ""), pattern = "^\\d+_\\d+_S\\d+_", replacement = ""))


# initalize myAnnot
myAnnot <- data.frame(raw.file = unique(protein$raw.file), Fraction = whatpart)
myAnnot$'Grouping.Var' <- gsub(x = gsub(x = unique(protein$raw.file) , pattern = "\\d$", replacement = ""), pattern = "^\\d+_\\d+_S\\d+_", replacement = "")
myAnnot$CONTROL <- "T"

# what is control condition
# from: https://fgcz-bfabric.uzh.ch/bfabric/order/show.html?id=35157 # -> see table
# starved = mock
table(myAnnot$Grouping.Var)
myAnnot$CONTROL[myAnnot$Grouping.Var == "HepG2Starvation"] <- "C"
table(myAnnot$Grouping.Var, myAnnot$CONTROL)


# split annot here in enriched and total part
annotfull <- myAnnot
(annot_enriched <- annotfull |> filter(Fraction == whatpart))

nr <- sum(annot_enriched$raw.file %in% unique(protein$raw.file))
logger::log_info("nr : ", nr, " files annotated")

# add missing info to protein
protein$"PeptideProphet.Probability" <- 1
protein$qValue <- 0

# join annotation and site centric output
protein <- dplyr::inner_join(annot_enriched, protein, multiple = "all")
head(protein)

# Grp 2
# work on GRP for having better folder name
(fN <- paste0(fgczProject,"_", descri, "_", WUID,"_",fracti))
GRP2_phos <- prolfquapp::make_DEA_config_R6(ZIPDIR = fN,PROJECTID = fgczProject,
                                            ORDERID = OIDfgcz)
# adjust roll-up
GRP2_phos$processing_options$aggregate <- "topN"

# final annotation table
annot <- prolfquapp::read_annotation(annot_enriched)


# Setup configuration
atable_phos <- annot$atable # take same annot table from total
atable_phos$ident_Score = "PeptideProphet.Probability"
atable_phos$ident_qValue = "qValue"
atable_phos$fileName <- "raw.file"
atable_phos$hierarchy[["protein_Id"]] <- c("proteinid")
atable_phos$hierarchy[["site"]] <- c("index")
#atable_phos$hierarchy[["peptide"]] <- c("peptide")
atable_phos$hierarchyDepth <- 2 # no roll-up to protein but to protNsite
atable_phos$set_response("intensity")


# Preprocess data - aggregate proteins.
config_phos <- prolfqua::AnalysisConfiguration$new(atable_phos)
#phospsm2 <- dplyr::inner_join(annot$annot, stypsmX, multiple = "all")
adata_phos <- prolfqua::setup_analysis(protein, config_phos)


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

# aggregating not needed since already site centric report (LFQ FP: combined_site_STY)
# lfqdata_phos <- prolfquapp::aggregate_data(lfqdata_phos, agg_method = GRP2_phos$processing_options$aggregate)

# modelling
grp_phos <- prolfquapp::generate_DEA_reports2(lfqdata_phos, GRP2_phos, protAnnot_phos, Contrasts = annot$contrasts)

# all fine? some checks
#grp_phos$RES$lfqData$to_wide()
#grp_phos$RES$contrastsData_signif
#myResPlotter <- grp_phos$RES$contrMerged$get_Plotter()
#myResPlotter$volcano()

logger::log_info("DONE WITH DEA REPORTS")

# result dir
(GRP2_phos$zipdir)
GRP2_phos$zipdir <- "p35157_phosphoOnly_WU304087"
dir.create(GRP2_phos$zipdir)

# need helper functions to properly write reports not on protein but peptide level
# source("FP_phosphoHelperFunctions_v3_202310.R")


copy_DEA_DIANN()
library(prophosqua)
copy_phosphoDEA_FragPipe_TMT()


GRP2 <- GRP2_phos

# writing reports
prolfquapp::write_DEA_all(grp2 = grp_phos, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")
prolfquapp::write_DEA_all(grp2 = grp_phos, boxplot = FALSE, markdown = "_DiffExpQC_Phospho_V2.Rmd")







