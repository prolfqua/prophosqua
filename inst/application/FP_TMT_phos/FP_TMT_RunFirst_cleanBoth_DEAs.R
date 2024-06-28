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

library(tidyverse)
library(prolfqua)
library(prolfquapp)
library(prophosqua)

# params ideally taken from yaml
fgczProject <- "pIDxx"
OIDfgcz <- "oxxx"
descri <- "FPguiTMTphospho_vsOneCondition"

################################################################################
#
#
#      psm -> Total part
#
#
################################################################################

fracti <- "TotalProteome"
WUID <- "WUID"


# data from
# # https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=301538&tab=details


# v3
prolfquapp::copy_DEA_DIANN()

# also copy the phospho specific Rmd files from prophosqua
prophosqua::copy_phosphoDEA_FragPipe_TMT()
#
path = "."

#ymlfile <- file.path(path,"config_FP_dea.yaml")
#GRP2 <- prolfquapp::read_BF_yamlR6(ymlfile, application = "DIANN")

# work on GRP for having better folder name
(fN <- paste0(fgczProject,WUID,"_",fracti))
GRP2 <- prolfquapp::make_DEA_config_R6(PATH = path,WORKUNITID = WUID, PROJECTID = fgczProject,
                                       ORDERID = OIDfgcz)


# for Witold: 28-06-2024
# this part is in the future taken from BF and the dataset
# for the time being we do manual parsing from sample names for condition and CONTROL column

dsf = file.path(path,"o34441_FP_tsvFiles/experiment_annotation.tsv") # this is directly from LK via FP-TMT GUI
#dsf <- readr::read_csv(dsf) # in BF it is csv -> make dsf from Bf

# make minimal dsf
dsf <- readr::read_tsv(dsf)
dsf
dsf$condition <- NULL
dsf$genotype <- NULL
dsf$treatment <- NULL
dsf$SampleName <- NULL
dsf$sample_name <- NULL
dsf$channel <- dsf$sample

# we further need CONTROL column and Grouping Var
# Grouping Var
table(gsub(x = dsf$channel, pattern = "_\\d+$", replacement = ""))
dsf$`Grouping Var` <- gsub(x = dsf$channel, pattern = "_\\d+$", replacement = "")
dsf$`Grouping Var` <- gsub(x = dsf$`Grouping Var`, pattern = "_\\d+_", replacement = "_")
dsf$`Grouping Var` <- gsub(x = dsf$`Grouping Var`, pattern = "WT2", replacement = "WT1")
dsf$`Grouping Var` <- gsub(x = dsf$`Grouping Var`, pattern = "_s$", replacement = "_st")
dsf$`Grouping Var` <- gsub(x = dsf$`Grouping Var`, pattern = "_S$", replacement = "_st")
table(dsf$`Grouping Var`)

# Control
dsf$CONTROL <- "T"
dsf$CONTROL[dsf$`Grouping Var` == "YZ135_WT1_untr"] <- "C"
table(dsf$CONTROL, dsf$`Grouping Var`)

# introduce a check? if all is there and correct?
# all contrasts based on CONTROL column
annotation <- read_annotation(dsf, prefix = "Group_")
table(dsf$channel)

path = "o34441_FP_tsvFiles/philosopher_prot/"
dir(path)
# get psm and fasta file
files <- get_FP_PSM_files(path = path)

# do config and preprocessing from psm for total
xd <- prolfquapp::preprocess_FP_PSM(quant_data = files$data,
                                    fasta_file = files$fasta,
                                    annotation = annotation,
                                    column_before_quants = "Mapped Proteins")


#
lfqdata <- xd$lfqdata
lfqdata$hierarchy_counts()
lfqdata$config$table$hierarchyDepth <- 2
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
grp_total <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2,
                                               xd$protein_annotation, annotation$contrasts)

logger::log_info("write results and html reports")

outpath <- prolfquapp::write_DEA_all(grp_total,  GRP2$zipdir , boxplot
                                     = FALSE, markdown = "_Grp2Analysis_V2.Rmd")

logger::log_info("write results and summarized experiment")
SE <- prolfquapp::make_SummarizedExperiment(grp_total)
saveRDS(SE, file = file.path(outpath ,
                             paste0("SummarizedExperiment",".rds") ))

# write out experimental design (annotation$annot)
(dsFN <- paste0("ExperimentAnnotation_DEA_WU", GRP2$project_spec$workunitID,".tsv"))
write_tsv(x = annotation$annot, file = file.path(outpath,dsFN))



################################################################################
#
#
#      Phospho Enriched Part
#
#
################################################################################


fracti <- "PhosphoEnriched"

multiSite_file <- dir(path = "o34441_FP_tsvFiles/tmt-report_phos//", pattern = "abundance_multi-site_None.tsv", recursive = TRUE, full.names = TRUE)
xx <- readr::read_tsv(multiSite_file)

# assuming that ReferenceIntensity is ALWAYS the first column before the individual quant channels are reported
(quant_idx_start <- grep(pattern = "ReferenceIntensity", x = colnames(xx))  + 1)

multiSiteAnnot <- xx[,1:(quant_idx_start-1)]
multiSiteQuant <- xx[,c(1,quant_idx_start:ncol(xx))]

# go long
multiSite_long <- tidyr::pivot_longer(data = multiSiteQuant, cols = 2:ncol(multiSiteQuant), values_to = "abundance", names_to = "sample")
multiSite_long <- dplyr::inner_join(x = multiSiteAnnot, y = multiSite_long)

# check
table(annotation$annot$CONTROL, annotation$annot$Grouping.Var)
table(annotation$annot$sample)

# in case of multi-site all the filtering has been done by FragPipe
# since TMT annotation is the same maybe the join has to be adapted!

# already done!
(myPhosPath <- paste0(fgczProject,WUID,"_",fracti))
GRP2_phos <- prolfquapp::make_DEA_config_R6(PATH = myPhosPath,WORKUNITID = WUID, PROJECTID = fgczProject,
                                            ORDERID = OIDfgcz)

# join with anno again this should work now with Name # if not all samples are used in the dataset they would be removed here (to be tested)
multiSite_long <- dplyr::inner_join(x = annotation$annot, y = multiSite_long)


# potentially filter out here files not  in the game (specified in CONTROL == NA)
multiSite_long <- multiSite_long |> filter(!is.na(CONTROL))
# show effect
multiSite_long |> select(sample, `Grouping.Var`) |> distinct() |> group_by(`Grouping.Var`) |> summarise(n())


# here we use Name to match annot and multiSite_long
nr <- sum(annotation$annot$sample %in% unique(multiSite_long$sample))
logger::log_info("nr : ", nr, " files annotated")

# add missing required parameters (qvalue)
multiSite_long$qvalue <- 1 - multiSite_long$MaxPepProb



# fasta and protein annotation part!
fasta_annot <- get_annot_from_fasta(files$fasta)
# reshape fasta_annot for matching protein
colnames(multiSite_long)
multiSite_long$nrPeptides <- 1
multiSite_long <- dplyr::left_join(multiSite_long, fasta_annot, by = c(ProteinID = "proteinname"), multiple = "all")
colnames(multiSite_long)

# Setup configuration manually for peptide analysis (phospho)
atable_phos <- prolfqua::AnalysisTableAnnotation$new()
atable_phos$factors <- annotation$atable$factors
atable_phos$ident_Score = "MaxPepProb"
atable_phos$ident_qValue = "qValue"
atable_phos$fileName = "channel"
atable_phos$hierarchy[["protein_Id"]] <- c("ProteinID")
atable_phos$hierarchy[["site"]] <- c("Index", "Peptide")
atable_phos$set_response("abundance")
atable_phos$hierarchyDepth <- 2
atable_phos$get_response()

# Preprocess data - aggregate proteins.
config_phos <- prolfqua::AnalysisConfiguration$new(atable_phos)
#multiSite_long$Modified.Peptide <- NA
#multiSite_long$Assigned.Modifications <- NA

adata_phos <- prolfqua::setup_analysis(multiSite_long, config_phos)
colnames(adata_phos)
nrow(adata_phos)
colnames(multiSite_long)

lfqdata_phos <- prolfqua::LFQData$new(adata_phos, config_phos)
lfqdata_phos$hierarchy_counts()
lfqdata_phos$remove_small_intensities(threshold = 1)
lfqdata_phos$hierarchy_counts()

lfqdata_phos$config$table$hierarchyDepth <- 2
lfqdata_phos$config$table$hierarchy_keys_depth() #
lfqdata_phos$config$table$hierarchyKeys() #checks if all is here
lfqdata_phos$config$table$ident_Score #checks if all is here
lfqdata_phos$config$table$ident_qValue #checks if all is here

lfqdata_phos$data #checks if all is here
lfqdata_phos$response() #checks if all is here


#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata_phos$config$table$hierarchy_keys()

#prolfqua::ProteinAnnotation$debug("initialize")
colnames(multiSite_long)
head(multiSite_long)

prot_annot_phos <- prolfquapp::build_protein_annot(lfqdata = lfqdata_phos, multiSite_long,
                                                   idcol = c(protein_Id = "ProteinID"), cleaned_protein_id = "ProteinID",
                                                   protein_description = "fasta.header", exp_nr_children = "nrPeptides",
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



# all fine 2024-04-30:: till here
#grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)
grp_phos <- prolfquapp::generate_DEA_reports2(lfqdata_phos, GRP2_phos, prot_annot_phos, GRP2_phos$pop$contrasts)

# all fine?
grp_phos$Groups_vs_Controls$RES$lfqData$to_wide()
grp_phos$Groups_vs_Controls$RES$contrastsData_signif
myResPlotter <- grp_phos$Groups_vs_Controls$RES$contrMerged$get_Plotter()
myResPlotter$volcano()

logger::log_info("DONE WITH DEA REPORTS")
# result dir
#GRP2_phos$zipdir_name <- paste0("DEA_",myPhosPath)
dir.create(GRP2_phos$path) 
dir.create(GRP2_phos$get_zipdir())


outpath <- prolfquapp::write_DEA_all(grp_phos, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho.Rmd")

# Save RData from enriched and total (only lfqdata is overwritten?) # keep lfqdata, grp, adata separate for phos and total!
(imageFN <- paste(fgczProject, descri, "total_and_enriched",".RData", sep="_"))
save.image(imageFN)
