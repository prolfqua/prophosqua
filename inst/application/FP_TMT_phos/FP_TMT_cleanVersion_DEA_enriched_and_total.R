# 2023-10-19: adapting for phosphoTMT total proteome and enriched
# 2024-03-26: running with latest prolfqua (1.1.5) and prlfquapp (0.1.9) before going prophosqua
# 2024-04-10: next version of DEA w/ prolfquapp 0.1.9
# 2024-04-25: towards clean version in prophosqua

# data and annotation based  on: https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details


# This script starts from psm.tsv for TotalProteome
# does basic filtering using purity.
#
# remotes::install_github('wolski/prolfquapp', dependencies = TRUE, force = FALSE)
# remotes::install_github('fgcz/prolfqua', dependencies = TRUE)


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
fgczProject <- "pIDxx"
OIDfgcz <- "oxxx"
descri <- "FPguiTMTphospho_vsOneCondition"
fracti <- "TotalProteome"
WUID <- "WUID"


# data from
# # https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=301538&tab=details


# v3
prolfquapp::copy_DEA_DIANN()
# also copy the phospho specific Rmd files from prophosqua
# prophosqua::copy_phosphoDEA_FragPipe_TMT() # not yet working, package not built?
#
path = "."

#ymlfile <- file.path(path,"config_FP_dea.yaml")
#GRP2 <- prolfquapp::read_BF_yamlR6(ymlfile, application = "DIANN")

# work on GRP for having better folder name
(fN <- paste0(fgczProject,WUID,"_",fracti))
GRP2 <- prolfquapp::make_DEA_config_R6(ZIPDIR = fN,PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)

dsf = file.path(path,"../RscriptNReports_nextLevel/o34441_FPguiTMTphospho__Dataset_TotalNEnriched_better.tsv")
#dsf <- readr::read_csv(dsf) # in BF it is csv -> make dsf from Bf
# make minimal dsf
dsf <- readr::read_tsv(dsf)
dsf$condition <- NULL
dsf$genotype <- NULL
dsf$treatment <- NULL
dsf$SampleName <- NULL
dsf$channel <- dsf$sample

# introduce a check? if all is there and correct?


# all contrasts based on CONTROL column
annotation <- read_annotation(dsf, prefix = "Group_")

path = "../o34441_FP_tsvFiles/philosopher_prot/"
dir(path)
# get psm and fasta file
files <- get_FP_PSM_files(path = path)

# do config and preprocessing from psm for total
xd <- prolfquapp::preprocess_FP_PSM(quant_data = files$data,
                                    fasta_file = files$fasta,
                                    annotation = annotation,
                                    column_before_quants = "Mapped Proteins")


lfqdata <- xd$lfqdata
lfqdata$hierarchy_counts()
lfqdata$config$table$hierarchyDepth <- 3
lfqdata$config$table$hierarchy_keys_depth()

lfqdata$config$table$ident_Score
lfqdata$config$table$ident_qValue

logger::log_info("AGGREGATING PSM to peptidoforms w/ Top1000!")
ag <- lfqdata$get_Aggregator()
ag$sum_topN(N = 1000)

# write aggregated out to lfqdata to continue!
lfqdata <- ag$lfq_agg
lfqdata$data
lfqdata$response()
lfqdata$hierarchy_counts()

lfqdata$config$table$hierarchyDepth <- 1

logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method =
                                        GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$processing_options$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

lfqdata$hierarchy_counts()


logger::log_info("run analysis")
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2,
                                         xd$protein_annotation, annotation$contrasts)

logger::log_info("write results and html reports")
prolfquapp::write_DEA_all(grp[[1]], names(grp)[1], GRP2$zipdir , boxplot
                          = FALSE, markdown = "_Grp2Analysis_V2.Rmd")

logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp[[1]])
saveRDS(SE, file = file.path(GRP2$zipdir,
                             paste0("Results_DEA_WU", grp[[1]]$project_spec$workunitID) ,
                             paste0("SummarizedExperiment",".rds") ))

# write out experimental design (annotation$annot)
(dsFN <- paste0("ExperimentAnnotation_DEA_WU", GRP2$project_spec$workunitID,".tsv"))
write_tsv(x = annotation$annot, file = file.path(GRP2$zipdir,dsFN))



#
# # work on annotation
# # structure of a dataset fix this in a template?
# # Relative Path,Name,Grouping Var,channel, CONTROL
# # for full proteome channel encodes the sample name like relative path in the dataset
# ds_file <- "../o34441_FP_tsvFiles/experiment_annotation.tsv"
#
# (annot <- read_tsv(ds_file))
#
# # drop some this is necessairy here because run via GUI and we need to generate DS from file names
# annot <- annot |> select(sample, channel, condition)
#
# # since not structured in BF some work needed
# annot$genotype <- c(rep("WT", 6), rep("cl2", 6), rep("cl20", 6))
# annot$treatment <- c(rep(c("untr", "S"), 3),rep(c("untr", "S"), 3),rep(c("untr", "S"), 3))
# annot$`Grouping Var` <- paste(annot$genotype, annot$treatment, sep = "_")
# annot$SampleName <- annot$sample
# #CONTROL column with WT_untr as C (all others as T)
# annot$CONTROL <- "T"
# #annot$CONTROL[annot$`Grouping Var` == "WT_untr"] <- "C"
# unique(annot$`Grouping Var`)
# annot$CONTROL[annot$`Grouping Var` == "WT_S"] <- "C"
# table(annot$CONTROL)
#
# (dsFN <- paste(fgczProject, descri,"Dataset_TotalNEnriched_better.tsv", sep="_"))
# write_tsv(x = annot, file = dsFN)
#


################################################################################
#
#
#      Phospho Enriched Part
#
#
################################################################################


fracti <- "PhosphoEnriched"
source("FP_phosphoHelperFunctions_v3_202310.R")

multiSite_file <- dir(path = "../o34441_FP_tsvFiles/tmt-report_phos//", pattern = "abundance_multi-site_None.tsv", recursive = TRUE, full.names = TRUE)
xx <- readr::read_tsv(multiSite_file)

# assuming that ReferenceIntensity is ALWAYS the first column before the individual quant channels are reported
(quant_idx_start <- grep(pattern = "ReferenceIntensity", x = colnames(xx))  + 1)

multiSiteAnnot <- xx[,1:(quant_idx_start-1)]
multiSiteQuant <- xx[,c(1,quant_idx_start:ncol(xx))]

# go long
multiSite_long <- tidyr::pivot_longer(data = multiSiteQuant, cols = 2:ncol(multiSiteQuant), values_to = "abundance", names_to = "sample")
multiSite_long <- dplyr::inner_join(x = multiSiteAnnot, y = multiSite_long)
head(multiSite_long)

# check
table(annotation$annot$CONTROL, annotation$annot$Grouping.Var)
table(annotation$annot$sample)

# in case of multi-site all the filtering has been done by FragPipe
# since TMT annotation is the same maybe the join has to be adapted!
#ds_file <- "../o32778_Anouk_FP-TMT-phospho/WU294902_phos_noNorm_correctAnnotation/dataset_from_combined_annotation_phospho_VS_16hpiNOGSK.tsv"
#(annot <- read_tsv(ds_file))
# zip specify here already
#(resDir <- paste(fgczProject, descri, fracti, sep="_"))

# already done!
myPhosZip <- paste0(fgczProject,WUID,"_",fracti)
GRP2_phos <- prolfquapp::make_DEA_config_R6(ZIPDIR = myPhosZip,PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)

# join with anno again this should work now with Name
multiSite_long <- dplyr::inner_join(x = annotation$annot, y = multiSite_long)


# again filter things that are in the comparisons
colnames(multiSite_long)
#multiSite_long$`Grouping Var` <- multiSite_long$Group
multiSite_long |> select(sample, `Grouping.Var`) |> distinct() |> group_by(`Grouping.Var`) |> summarise(n())

# filter for all that are in the game
multiSite_long |> filter(!is.na(CONTROL)) |> select(sample, `Grouping.Var`) |> distinct() |> group_by(`Grouping.Var`) |> summarise(n())
# potentially filter out here files not  in the game (specified in CONTROL == NA)
multiSite_long <- multiSite_long |> filter(!is.na(CONTROL))


# here we use Name to match annot and multiSite_long
nr <- sum(annotation$annot$sample %in% unique(multiSite_long$sample))
logger::log_info("nr : ", nr, " files annotated")

# add missing required parameters (qvalue)
multiSite_long$qvalue <- 1 - multiSite_long$MaxPepProb


# Do I really need this? Or it is already in the grp2
#tmp <- prolfquapp::dataset_set_factors_deprecated(atable, multiSite_long)
#atable <- tmp$atable
#atable$factors
#multiSite_long <- tmp$msdata


# CREATE protein annotation. -> 2024-04-11: deprecated! use build_protein_annot
# prot_annot <- prolfquapp::dataset_protein_annot(
#   multiSite_long,
#   c("protein_Id" = "ProteinID"),
#   protein_annot = "Gene",
#   more_columns = c())

# fasta and protein annotation part!
fasta_annot <- get_annot_from_fasta(files$fasta)
# reshape fasta_annot for matching protein
# fasta_annot$proteinAcc <- sapply(strsplit(fasta_annot$fasta.id, split = "\\|"), function(x)x[2])
colnames(multiSite_long)
multiSite_long$nrPeptides <- 1
multiSite_long <- dplyr::left_join(multiSite_long, fasta_annot, by = c(ProteinID = "proteinname"), multiple = "all")

# Setup configuration
#atable <- annotation$atable # this one is from total and has things in we dont want
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$factors <- annotation$atable$factors
atable$ident_Score = "MaxPepProb"
atable$ident_qValue = "qValue"
atable$fileName = "channel"
atable$hierarchy[["protein_Id"]] <- c("ProteinID")
atable$hierarchy[["site"]] <- c("Index", "Peptide")
atable$set_response("abundance")
atable$hierarchyDepth <- 2
atable$get_response()

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
#multiSite_long$Modified.Peptide <- NA
#multiSite_long$Assigned.Modifications <- NA

adata <- prolfqua::setup_analysis(multiSite_long, config)
colnames(adata)
nrow(adata)
colnames(multiSite_long)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$hierarchy_counts()
lfqdata$remove_small_intensities(threshold = 1)
lfqdata$hierarchy_counts()

lfqdata$config$table$hierarchyDepth <- 2
lfqdata$config$table$hierarchy_keys_depth() # @Witold -> how to model this? set here hierarchy depth to peptide?
lfqdata$config$table$hierarchyKeys() #checks if all is here
lfqdata$config$table$ident_Score #checks if all is here
lfqdata$config$table$ident_qValue #checks if all is here

lfqdata$data #checks if all is here
lfqdata$response() #checks if all is here


#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata$config$table$hierarchy_keys()

#prolfqua::ProteinAnnotation$debug("initialize")
colnames(multiSite_long)
head(multiSite_long)
#prot_annot <- prolfquapp::dataset_protein_annot(multiSite_long, c(protein_Id = "ProteinID"), protein_annot = "fasta.header", more_columns = "nrPeptides")
# protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)  # again here ;( Error in `sym()`: ! Can't convert a character vector to a symbol.
prot_annot <- prolfquapp::build_protein_annot(lfqdata, multiSite_long,
                                              idcol = c(protein_Id = "ProteinID"), cleaned_protein_id = "ProteinID",
                                              protein_description = "fasta.header", nr_children = "nrPeptides",
                                              more_columns = "fasta.id")


# important for phospho
GRP2_phos$pop$nr_peptdes <- 1
GRP2_phos$pop$aggregate <- "none"
GRP2_phos$pop$contrasts <- annotation$contrasts


logger::log_info("GENERATING DEA REPORTS")
logger::log_info("starting modelling")

# fix for old grps
GRP2_phos$pop$aggregate
GRP2_phos$processing_options$transform <- "robscale"
GRP2_phos$pop$transform <- GRP2_phos$processing_options$transform
GRP2_phos$pop$Diffthreshold <- GRP2_phos$processing_options$diff_threshold
GRP2_phos$pop$FDRthreshold <- GRP2_phos$processing_options$FDR_threshold



# something is  not processing properly
#grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)
grp_phos <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2_phos, prot_annot, GRP2_phos$pop$contrasts)

# all fine?
grp_phos$Groups_vs_Controls$RES$lfqData$to_wide()
grp_phos$Groups_vs_Controls$RES$contrastsData_signif
myResPlotter <- grp_phos$Groups_vs_Controls$RES$contrMerged$get_Plotter()
myResPlotter$volcano()

logger::log_info("DONE WITH DEA REPORTS")
# result dir
# now we prefer it to have it directly in the workingDir -> specified before for enriched
dir.create(GRP2_phos$zipdir)

# need helper functions to properly write reports not on protein but peptide level
source("FP_phosphoHelperFunctions_v3_202310.R")

#GRP2_phos$pop$DiffThreshold
#grp$Groups_vs_Controls$pop$DiffThreshold
# write DEAs
# mydebug
#params = list(grp = grp_phos[[1]])
# missing before render
#grp2$pop$LocProbThresh <- 0.75

#
for (i in seq_along(grp_phos)) {
  #prolfquapp::write_DEA_all(grp[[i]], names(grp)[i], GRP2$zipdir, boxplot = FALSE)
  write_phosphoDEA_all(grp_phos[[i]], names(grp_phos)[i], GRP2_phos$zipdir, boxplot = FALSE)
}

# Save RData from enriched and total (only lfqdata is overwritten?) # keep lfqdata, grp, adata separate for phos and total!
(imageFN <- paste(fgczProject, descri, "total_and_enriched",".RData", sep="_"))
save.image(imageFN)



