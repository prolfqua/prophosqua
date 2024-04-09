# 2023-10-19: adapting for phosphoTMT total proteome
# 2024-03-26: running with latest prolfqua (1.1.5) and prlfquapp (0.1.9) before going prophosqua

# annotation based  on: https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details


# This script starts from psm.tsv
# does basic filtering using purity.
#
# remotes::install_github('wolski/prolfquapp', dependencies = TRUE, force = TRUE)
#remotes::install_github('fgcz/prolfqua', dependencies = TRUE)


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

# params ideally taken from yaml
fgczProject <- "o34441"
OIDfgcz <- "o34441"
descri <- "FPguiTMTphospho_"
fracti <- "TotalProteome"
WUID <- "WUxx"

source("FP_phosphoHelperFunctions_v3_202310.R") # for tidy_psm_v2


# work on annotation
# structure of a dataset fix this in a template?
# Relative Path,Name,Grouping Var,channel, CONTROL
# for full proteome channel encodes the sample name like relative path in the dataset
ds_file <- "../o34441_FP_tsvFiles/experiment_annotation.tsv"

(annot <- read_tsv(ds_file))
# drop some this is necessairy here because run via GUI
annot <- annot |> select(sample, channel, condition)

# since not structured in BF some work needed
annot$genotype <- c(rep("WT", 6), rep("cl2", 6), rep("cl20", 6))
annot$treatment <- c(rep(c("untr", "S"), 3),rep(c("untr", "S"), 3),rep(c("untr", "S"), 3))
annot$`Grouping Var` <- paste(annot$genotype, annot$treatment, sep = "_")
annot$SampleName <- annot$sample
#CONTROL column with WT_untr as C (all others as T)
annot$CONTROL <- "T"
annot$CONTROL[annot$`Grouping Var` == "WT_untr"] <- "C"

(dsFN <- paste(fgczProject, descri,"Dataset_TotalNEnriched_better.tsv", sep="_"))
write_tsv(x = annot, file = dsFN)


#
#  Total Proteome starting from psm.tsv
#
#

# sanitize peptide csv.
psm_file <- dir(path = "../o34441_FP_tsvFiles/philosopher_prot/", pattern = "psm.tsv", recursive = TRUE, full.names = TRUE)
xx <- readr::read_tsv(psm_file)

# other columns that might be of interest at one point? at the moment we are NOT using these at all simply trusting on FP output!
# NoMC, NoET, ProteinStart, ProteinEnd, Purity, IsUnique, EntryName, Gene, Mapped Genes, Mapped Proteins 
colnames(xx)[c(22:25,29,30,33:37)]
head(xx[,c(22:25,29,30,33:37)])

# overview on mod tables -> useful to check labeling efficiency? nope -> keep?
xa <- xx$`Assigned Modifications`
tmp <- gsub(" ", "", unlist(str_split(xa, ",")))
tmp <- gsub("^[0-9]+","", tmp)
table(tmp)

# look at some scores first?
# scores <- xx |> select(all_of(c("Expectation","Hyperscore","Nextscore","PeptideProphet Probability"))) # keep?
# image(cor(scores))

# properly parse psm file too long
psm <- tidy_FragPipe_psm_v2(psm_file, lastColumnBeforeQuants = "Mapped Proteins")
psm$qValue <- 1 - psm$PeptideProphet.Probability


# now we prefer it to have it directly in the workingDir
(resDir <- paste(fgczProject, descri, fracti, sep="_"))
GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = resDir, Normalization = "robscale", PROJECTID = fgczProject, ORDERID = OIDfgcz, WORKUNITID = WUID)
dir.create(GRP2$zipdir)


# find  all contrasts
GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2)

channelCol  <- grep("^channel", names(annot), ignore.case = TRUE, value = TRUE)
if (channelCol != "channel") {
  annot[["channel"]] <- annot[[channelCol]]
  annot[[channelCol]] <- NULL
}

colnames(psm)
colnames(annot)
unique(psm$channel)
# join psm and annot by channel or by SampleName
# combo <- left_join(x = phosRes, y = totRes, join_by("proteinid" == "protein_Id", "contrast" == "contrast"))
#psm <- dplyr::inner_join(x = annot, y = psm, by = "channel")
psm$SampleName <- psm$channel # channels were annotated in the GUI 
psm <- dplyr::inner_join(x = annot, y = psm, by = "SampleName")
nr <- sum(annot$SampleName %in% unique(psm$SampleName))
logger::log_info("nr : ", nr, " files annotated")

#work on some columns
colnames(psm)
colnames(psm)[23]
psm <- psm[,-23]
psm$channel <- psm$channel.x
colnames(psm)
psm <- psm[,-2]
psm$Name <-  psm$sample
unique(psm$Name)

psm |> select(Name, `Grouping Var`) |> distinct() |> group_by(`Grouping Var`) |> summarise(n())

# filter for all that are in the game
psm |> filter(!is.na(CONTROL)) |> select(Name, `Grouping Var`) |> distinct() |> group_by(`Grouping Var`) |> summarise(n())
psm <- psm |> filter(!is.na(CONTROL))
# there is NO control column, we get all against all?


# Setup configuration
atable <- prolfqua::AnalysisTableAnnotation$new()

atable$ident_Score = "PeptideProphet.Probability"
atable$ident_qValue = "qValue"
atable$fileName = "channel"
atable$hierarchy[["protein_Id"]] <- c("Protein")
atable$hierarchy[["peptide_Id"]] <- c("Peptide")
atable$hierarchy[["mod_peptide_Id"]] <- c("Modified.Peptide","Assigned.Modifications")
atable$hierarchy[["Spectrum"]] <- c("Spectrum")

#
tmp <- prolfquapp::dataset_set_factors_deprecated(atable, psm)
atable <- tmp$atable
atable$factors
psm <- tmp$msdata

head(psm)

# CREATE protein annotation.
prot_annot <- prolfquapp::dataset_protein_annot(
  psm,
  c("protein_Id" = "Protein"),
  protein_annot = "Protein.Description",
  more_columns = "nrPeptides")


atable$set_response("abundance")
# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
colnames(psm)
adata <- prolfqua::setup_analysis(psm, config)
colnames(adata)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities(threshold = 1)
lfqdata$hierarchy_counts()
xx <- lfqdata$summarize_hierarchy()
xx$Spectrum_n |> max()

lfqdata$config$table$hierarchyDepth <- 3
lfqdata$config$table$hierarchy_keys_depth()

lfqdata$config$table$ident_Score
lfqdata$config$table$ident_qValue

ag <- lfqdata$get_Aggregator()

ag$sum_topN(N = 10000)

lfqdata <- ag$lfq_agg
lfqdata$data
lfqdata$response()
lfqdata$hierarchy_counts()


#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata$config$table$hkeysDepth()
lfqdata$config$table$hierarchyDepth <- 1
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)

#logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
#lfqdata$factors()

logger::log_info("END OF DATA TRANSFORMATION.")
                     #debug(prolfquapp::generate_DEA_reports)
protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)

grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, protAnnot) # this is taking quite a while

# change output Dir
# "../o32778_Anouk_FP-TMT-phospho/WU294264_FragPipeTMT_m1_prot/dataset_from_combined_annotation_VS_5hpiNOGSK"
# -> originally - right next to the dataset


for (i in seq_along(grp)) {
  prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = FALSE)
}

(imageFN <- paste(fgczProject, descri,"total.RData", sep="_"))
save.image(imageFN)

#prolfquapp::render_DEA(grp2, outpath = ".", htmlname = "check.html")

################################################################################
#
#
#      Phospho Enriched Part
#
#
################################################################################
rm(list=ls())
# params ideally taken from yaml
fgczProject <- "o34441"
descri <- "FPguiTMTphospho_"
fracti <- "PhosphoEnriched"
source("FP_phosphoHelperFunctions_v3_202310.R")

multiSite_file <- dir(path = "../o34441_FP_tsvFiles/tmt-report_phos//", pattern = "abundance_multi-site_None.tsv", recursive = TRUE, full.names = TRUE)
xx <- readr::read_tsv(multiSite_file)
colnames(xx)

# assuming that ReferenceIntensity is ALWAYS the first column before the individual quant channels are reported
(quant_idx_start <- grep(pattern = "ReferenceIntensity", x = colnames(xx))  + 1)

multiSiteAnnot <- xx[,1:(quant_idx_start-1)]
multiSiteQuant <- xx[,c(1,quant_idx_start:ncol(xx))]
# go long
multiSite_long <- tidyr::pivot_longer(data = multiSiteQuant, cols = 2:ncol(multiSiteQuant), values_to = "abundance", names_to = "SampleName")
multiSite_long <- dplyr::inner_join(x = multiSiteAnnot, y = multiSite_long)
head(multiSite_long)

annot <- read_tsv("o34441_FPguiTMTphospho__Dataset_TotalNEnriched.tsv")


# in case of multi-site all the filtering has been done by FragPipe
# since TMT annotation is the same maybe the join has to be adapted!
#ds_file <- "../o32778_Anouk_FP-TMT-phospho/WU294902_phos_noNorm_correctAnnotation/dataset_from_combined_annotation_phospho_VS_16hpiNOGSK.tsv"
#(annot <- read_tsv(ds_file))

# zip specify here already
(resDir <- paste(fgczProject, descri, fracti, sep="_"))

GRP2 <- prolfquapp::make_DEA_config(ZIPDIR = resDir, Normalization = "robscale", PROJECTID = fgczProject, ORDERID = fgczProject, WORKUNITID = "WU295758_Phospho")


# what happens here?
GRP2 <- prolfquapp::dataset_extract_contrasts(annot, GRP2)

# join with anno again this should work now with Name
multiSite_long <- dplyr::inner_join(x = annot, y = multiSite_long)


# again filter things that are in the comparisons
colnames(multiSite_long)
#multiSite_long$`Grouping Var` <- multiSite_long$Group
multiSite_long |> select(SampleName, `Grouping Var`) |> distinct() |> group_by(`Grouping Var`) |> summarise(n())

# filter for all that are in the game
multiSite_long |> filter(!is.na(CONTROL)) |> select(SampleName, `Grouping Var`) |> distinct() |> group_by(`Grouping Var`) |> summarise(n())

multiSite_long <- multiSite_long |> filter(!is.na(CONTROL))


# here we use Name to match annot and multiSite_long
nr <- sum(annot$SampleName %in% unique(multiSite_long$SampleName))
logger::log_info("nr : ", nr, " files annotated")

colnames(multiSite_long)


# add missing required parameters (qvalue)
multiSite_long$qvalue <- 1 - multiSite_long$MaxPepProb

# Setup configuration
atable <- prolfqua::AnalysisTableAnnotation$new()

atable$ident_Score = "MaxPepProb"
atable$ident_qValue = "qValue"
atable$fileName = "Name"
atable$hierarchy[["protein_Id"]] <- c("ProteinID")
atable$hierarchy[["site"]] <- c("Index", "Peptide")

atable$hierarchyDepth <- 2

#
tmp <- prolfquapp::dataset_set_factors_deprecated(atable, multiSite_long)
atable <- tmp$atable
atable$factors
multiSite_long <- tmp$msdata

head(multiSite_long)

# CREATE protein annotation.
prot_annot <- prolfquapp::dataset_protein_annot(
  multiSite_long,
  c("protein_Id" = "ProteinID"),
  protein_annot = "Gene",
  more_columns = c())


atable$set_response("abundance")

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
colnames(multiSite_long)
adata <- prolfqua::setup_analysis(multiSite_long, config)
colnames(adata)
nrow(adata)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities(threshold = 1)
lfqdata$hierarchy_counts()

lfqdata$config$table$hierarchyDepth <- 2
lfqdata$config$table$hierarchy_keys_depth() # @Witold -> how to model this? set here hierarchy depth to peptide?
lfqdata$config$table$hierarchyKeys()
lfqdata$config$table$ident_Score
lfqdata$config$table$ident_qValue

lfqdata$data
lfqdata$response()
lfqdata$hierarchy_counts()


#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata$config$table$hierarchy_keys()
#lfqdata$config$table$hierarchyDepth <- 2
#lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method = GRP2$pop$aggregate)
#logger::log_info("data aggregated: {GRP2$pop$aggregate}.")
#lfqdata$factors()

#debug(prolfquapp::generate_DEA_reports)

#prolfqua::ProteinAnnotation$debug("initialize")
protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)
# important for phospho
GRP2$pop$nr_peptdes <- 1
GRP2$pop$aggregate <- "none"

logger::log_info("GENERATING DEA REPORTS")
logger::log_info("starting modelling")

grp <- prolfquapp::generate_DEA_reports(lfqdata, GRP2, protAnnot)

logger::log_info("DONE WITH DEA REPORTS")

# result dir
# now we prefer it to have it directly in the workingDir -> specified before for enriched
dir.create(GRP2$zipdir)

# need helper functions to properly write reports not on protein but peptide level
source("FP_phosphoHelperFunctions_v3_202310.R")


# write DEAs
for (i in seq_along(grp)) {
  #prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = FALSE)
  write_phosphoDEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = TRUE)
}

(imageFN <- paste(fgczProject, descri, "enriched",".RData", sep="_"))
save.image(imageFN)
