#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

#remotes::install_github('wolski/prolfquapp', dependencies = TRUE)
#remotes::install_github('fgcz/prolfqua', dependencies = TRUE)

library(prolfqua)
library(prolfquapp)
library(prophosqua)
library(readr)
library(dplyr)
library(openxlsx)

# source FP phospho helper Functions to use with preprocess multiplex
# source("prophosqua/R/FP_phospho_helper_functions_v3_202310.R")

# Integration of 2-plex phospho-TMT data
# how does the multi-site export look like -> for two plexes?

# variables
fgczProject <- "pxxxx"
WUID  <- "vs_oneCondition"
OIDfgcz <- "oyyy"


path <- "."
# v3
prolfquapp::copy_DEA_Files()
 # also copy the phospho specific Rmd files from prophosqua
prophosqua::copy_phosphoDEA_FragPipe_TMT()


# Look at annotation file provided and adapt it to the needs of prolfquapp
dsf <- "SampleAnnotation_2plex.txt"

# make minimal dsf
dsf <- readr::read_tsv(dsf)
dsf
dsf$Condition <- NULL
dsf$Timing <- NULL
dsf$channel <- dsf$SampleID
dsf$fgczLabeling <- dsf$channel_plex
dsf$G_ch <- NULL
dsf$plex <- NULL
dsf$channel_plex <- NULL
dsf$Control <- NULL
dsf$sample <- dsf$SampleID

# we further need CONTROL column and Grouping Var
# Grouping Var
dsf$`Grouping Var` <- dsf$GroupingVar
dsf$GroupingVar <- NULL
table(dsf$`Grouping Var`)

# Control
dsf$CONTROL <- "T"
dsf$CONTROL[dsf$`Grouping Var` == WUID] <- "C"
table(dsf$CONTROL, dsf$`Grouping Var`)

# get annotation objectt
annotation <- read_annotation(dsf, prefix = "Group_")

# write out annotation table
(fN <- paste0(fgczProject,"_",WUID,"_annotationTable.txt"))
write_tsv(x = annotation$annot, file = fN)

# also write to xlsx for generating yaml file afterwards
(fNxls <- paste0(fgczProject,"_vs_",WUID,"_annotationTable.xlsx"))
openxlsx::write.xlsx(annotation$annot[,-grep(x = colnames(annotation$annot),pattern = "idx")], file = fNxls)


################################################################################
#
#
#      Total Proteome
#
#
################################################################################
# starting here with total (psm.tsv) -> if multiple plexes are present there are multiple tsvs

fracti <- "total"
# work on GRP for having better folder name
(fN <- paste0(fgczProject,"_",WUID,"_",fracti))
GRP2 <- prolfquapp::make_DEA_config_R6(PATH = path,WORKUNITID = WUID, PROJECTID = paste0(fracti,"_",fgczProject),
                                       ORDERID = "")
# check result folder
GRP2$zipdir_name

# get psm for each plex and fasta file
psmF1 <- "../Fragpipe_o35920_proteome/proteome_psm_plex1.tsv"
psmF2 <- "../Fragpipe_o35920_proteome/proteome_psm_plex2.tsv"

# read each psm file individually (mit witold besprochen!)
# then feed it back to adapted function
psm1 <- tidy_FragPipe_psm(psm_file = psmF1)
psm2 <- tidy_FragPipe_psm(psm_file = psmF2)

#join psm objects (how to handle nr_Peptides_exp?) when same proteins are identified?
# Combine tibbles into a list
nrPep_tibble_list <- list(psm1$nrPeptides_exp, psm2$nrPeptides_exp)

# only take max peptide if found in both plexes
uNrProtein_with_maxPep <- bind_rows(nrPep_tibble_list) %>%
  group_by(Protein) %>%
  summarize(nrPeptides = max(nrPeptides)) %>%
  ungroup()

psm_all <- list(data = rbind(psm1$data, psm2$data), nrPeptides_exp = uNrProtein_with_maxPep)

nrPeptides_exp <- psm_all$nrPeptides
psm_all$data$qValue <- 1 - psm_all$data$Probability

# fasta file
fastaf <- "../Fragpipe_phospho/myUP000005640.fasta"

# how does it look
unique(psm_all$data$channel)
annotation$annot$channel
annotation$atable$fileName

# do preprocessing and go long
xd <- preprocess_FP_multiplexPSM(psm = psm_all, fasta_file = fastaf, annotation = annotation, column_before_quants = "ReferenceIntensity", pattern_decoys = "rev_")

# hand over to prolfqua and do config by hand
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

# check results before rendering
# totRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id)
dim(grp_total$RES$contrastsData)
grp_total$RES$contrastsData |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id) # this should only be decoy proteins

logger::log_info("write results and html reports")
outpath <- prolfquapp::write_DEA_all(grp_total,  GRP2$zipdir_name , boxplot
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

# analyzing multi-site phospho with shell scripts
# https://github.com/prolfqua/prolfquapp/tree/master

# get shell scripts here -> use Terminal
# navigate first to proper directory
# R --vanilla -e "prolfquapp::copy_shell_script()"
# do not forget to add executable permissions
# chmod a+x CMD_*

# get help for CMD_MAKE_YAML.R --help
# generate yaml file
# mkdir PhosphoDEA_viaSH # not necessairy
# bash CMD_MAKE_YAML.sh --trans robscale --outdir PhosphoDEA_outDir --workunit WUxxx -p xxx -o yyy -s FP_multisite --yaml minimalPhosphoAnalysis.yaml

# dataset
# we do this dataset before but this would also be a good starting point
# bash CMD_MAKE_DATASET.sh --help


# generate output directory
mydate <- format(Sys.Date(), "%Y%m%d")
(outdir <- paste0("DEA_", mydate,"_",fracti, "_",fgczProject, "_WU_",WUID))


# generate yaml script line
#bash CMD_MAKE_YAML.sh --norm robscale --outdir PhosphoDEA_outDir --workunit WUID -p xxx -O yyy -s FP_multisite --yaml minimalPhosphoAnalysis22.yaml
(ymlF <- paste0("minimalYaml_vs",WUID,".yaml"))
(ymlCMD <- paste0("bash CMD_MAKE_YAML.sh --norm robscale --outdir ",outdir," --workunit ", WUID," -p ",fgczProject," -O ",OIDfgcz," -s FP_multisite --yaml ",ymlF))

# run yaml script
system(ymlCMD)

# Potential issue here:
#    pattern_decoys: ^REV_ ->  pattern_decoys: ^rev_ # pattern to identify decoy proteins cannot be changed?
#    do this manually in yaml file -> otherwise in protein annotation fasta.id contains rev_


# run DEA
# run shell script
#bash CMD_DEA.sh -i ../Fragpipe_phospho/ -d "pxxx_WUxxx_annotationTableXLSX.xlsx" -s FP_multisite -o PhosphoDEA_outDir -y PhosphoDEA_outDir/minimalPhosphoAnalysis.yaml
# build up shell script in R
firstPartStable <- "bash CMD_DEA.sh -i ../Fragpipe_phospho/ -d "
middlePart <- " -s FP_multisite -o "
(deaCMD <- paste0(firstPartStable, fNxls, middlePart, outdir, " -y ", outdir,"/",ymlF))

# run DEA
system(deaCMD)








