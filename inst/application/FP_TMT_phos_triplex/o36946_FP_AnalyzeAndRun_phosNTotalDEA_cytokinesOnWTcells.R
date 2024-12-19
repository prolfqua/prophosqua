#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

#remotes::install_github('fgcz/prolfqua', dependencies = TRUE)
#remotes::install_github('wolski/prolfquapp', dependencies = TRUE)

library(prolfqua)
library(prolfquapp)
library(prophosqua)
library(prolfquappPTMreaders)
library(readr)
library(dplyr)
library(openxlsx)

# source FP phospho helper Functions to use with preprocess multiplex
source("~/GitHub/prophosqua/R/FP_phosphoHelperFunctions_v3_202310.R")

# variables
fgczProject <- "p36946"
WUID  <- "_CytokinesOnWTCells"
OIDfgcz <- "o36946"


path <- "."
# v3
# prolfquapp::copy_DEA_Files()
# also copy the phospho specific Rmd files from prophosqua
# prophosqua::copy_phosphoDEA_FragPipe_TMT() # also the _wew is needed

# this only needs to be run once!
# generate dsf from annotation in headers
# dat <- read_tsv("../o36946_FP_TMTi_enriched/tmt-report/abundance_multi-site_None.tsv")
# colnames(dat)
# myFiles <- colnames(dat)[10:63]
# plex <- gsub(x = myFiles, pattern = "(rep[A|B|C])_.*", replacement = "\\1")
# SampleID <- myFiles
# allConditions <- gsub(x = myFiles, pattern = "(rep[A|B|C])_", replacement = "")
# genotype <- gsub(x=gsub(x = allConditions, pattern = "(.*)_.*", replacement = "\\1"), pattern = "(.*)_.*", replacement = "\\1")
# condition <- gsub(x = allConditions, pattern = ".*_(.*)", replacement = "\\1")
# GroupingVar <- allConditions

# generate dataframe for sample annotation
# dsf <- data.frame(SampleID = SampleID, plex = plex, genotype = genotype, GroupingVar = GroupingVar)
# table(dsf$GroupingVar)
# dsfN <- "SampleAnnotation_multiplexed.txt"
# write_tsv(x = dsf, file = dsfN)

# Look at annotation file provided and adapt it to the needs of prolfquapp

# dsfN <- "SampleAnnotation_multiplexed.txt"
# # # make minimal dsf
# dsf <- readr::read_tsv(dsfN)
# dsf
# dsf$channel <- dsf$SampleID
# dsf$fgczLabeling <- dsf$plex
# dsf$Subject <- dsf$plex
# dsf$channel_plex <- NULL
# dsf$Control <- NULL
# dsf$sample <- dsf$SampleID
#
# we further need CONTROL column and Grouping Var
# # Grouping Var
# dsf$`Grouping Var` <- dsf$GroupingVar
# dsf$GroupingVar <- NULL
# table(dsf$`Grouping Var`)
#
# # Control
# dsf$CONTROL <- "T"
# table(dsf$`Grouping Var`)
# dsf$CONTROL[dsf$`Grouping Var` == "WT_ctrl"] <- "C"
# table(dsf$CONTROL, dsf$`Grouping Var`)

# inject here the xls file that is modified with self specified contrasts
myFNxlsx <- "p36946_selfspecifiedContrasts_annotationTable_cytokinesEffect.xlsx"
mydsf <- openxlsx::read.xlsx(xlsxFile = myFNxlsx, sheet = 1)
str(mydsf)

# get annotation objectt
annotation <- read_annotation(mydsf, prefix = "G_")
annotation$annot
annotation$contrasts

# write out annotation table
(fN <- paste0(fgczProject,"_",WUID,"_annotationTable.txt"))
write_tsv(x = annotation$annot, file = fN)

# also write to xlsx for generating yaml file afterwards
(fNxls <- paste0(fgczProject,WUID,"_annotationTable_rr.xlsx"))
#openxlsx::write.xlsx(annotation$annot[,-grep(x = colnames(annotation$annot),pattern = "idx")], file = fNxls)
openxlsx::write.xlsx(annotation$annot, file = fNxls)


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
psmF1 <- "../o36946_FP_TMTi_PSMtables_input/repA_input/psm.tsv"
psmF2 <- "../o36946_FP_TMTi_PSMtables_input/repB_input/psm.tsv"
psmF3 <- "../o36946_FP_TMTi_PSMtables_input/repC_input/psm.tsv"

# read each psm file individually (mit witold besprochen!)
# then feed it back to adapted function
psm1 <- tidy_FragPipe_psm(psm_file = psmF1)
psm2 <- tidy_FragPipe_psm(psm_file = psmF2)
psm3 <- tidy_FragPipe_psm(psm_file = psmF3)

#join psm objects (how to handle nr_Peptides_exp?) when same proteins are identified?
# Combine tibbles into a list
nrPep_tibble_list <- list(psm1$nrPeptides_exp, psm2$nrPeptides_exp, psm3$nrPeptides_exp)

# only take max peptide if found in both plexes
uNrProtein_with_maxPep <- bind_rows(nrPep_tibble_list) %>%
  group_by(Protein) %>%
  summarize(nrPeptides = max(nrPeptides)) %>%
  ungroup()

psm_all <- list(data = rbind(psm1$data, psm2$data, psm3$data), nrPeptides_exp = uNrProtein_with_maxPep)

nrPeptides_exp <- psm_all$nrPeptides
psm_all$data$qValue <- 1 - psm_all$data$Probability

# fasta file
fastaf <- "../o36946_FP_TMTi_enriched/2024-12-02-decoys-contam-UP000000589.fasta"

# how does it look
unique(psm_all$data$channel)
annotation$annot$channel
annotation$atable$fileName

# do preprocessing and go long # for this function the r helper file should be sourced
xd <- preprocess_FP_multiplexPSM(psm = psm_all, fasta_file = fastaf, annotation = annotation, column_before_quants = "Mapped Proteins", pattern_decoys = "rev_")

# hand over to prolfqua and do config by hand
lfqdata <- xd$lfqdata
lfqdata$hierarchy_counts()
lfqdata$config$table$hierarchyDepth <- 2
lfqdata$config$table$hierarchy_keys_depth()

lfqdata$config$table$ident_Score
lfqdata$config$table$ident_qValue

logger::log_info("AGGREGATING PSM to peptidoforms w/ Top1000!")
ag <- lfqdata$get_Aggregator()
ag$sum_topN(N = 10000)

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

grp_total$RES$contrastsData$contrast |> unique()

# check results before rendering
# totRes |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id)
dim(grp_total$RES$contrastsData)
grp_total$RES$contrastsData |> dplyr::filter(is.na(protein_length)) |> dplyr::pull(protein_Id) # this should only be decoy proteins

logger::log_info("write results and html reports")
outpath <- prolfquapp::write_DEA_all(grp_total,  GRP2$zipdir_name , boxplot
                                     = FALSE, markdown = "_Grp2Analysis_V2.Rmd")

# logger::log_info("write results and summarized experiment")
# SE <- prolfquapp::make_SummarizedExperiment(grp_total)
# saveRDS(SE, file = file.path(outpath ,
#                              paste0("SummarizedExperiment",".rds") ))

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
# system('R --vanilla -e "prolfquapp::copy_shell_script()"')
# do not forget to add executable permissions
# chmod a+x prolfqua_*
# system("chmod a+x prolfqua_*.sh")

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
mkDirCMD <- paste("mkdir", outdir)
system(mkDirCMD)
softwareHere <- "FP_multisite"
(ymlCMD <- paste0("bash prolfqua_yaml.sh --norm robscale --outdir ",outdir," --workunit ", WUID," -p ",fgczProject," -O ",OIDfgcz," -s ", softwareHere, " --yaml ",ymlF))

# run yaml script
system(ymlCMD)

# Potential issue here:
#    pattern_decoys: ^REV_ ->  pattern_decoys: ^rev_ # pattern to identify decoy proteins cannot be changed?
#    do this manually in yaml file -> otherwise in protein annotation fasta.id contains rev_


# run DEA
# run shell script
#bash CMD_DEA.sh -i ../Fragpipe_phospho/ -d "pxxx_WUxxx_annotationTableXLSX.xlsx" -s FP_multisite -o PhosphoDEA_outDir -y PhosphoDEA_outDir/minimalPhosphoAnalysis.yaml
# build up shell script in R
myPhosphoFolder <- "../o36946_FP_TMTi_enriched/"
firstPartStable <- paste0("bash prolfqua_dea.sh -i ",myPhosphoFolder ," -d ")
middlePart <- paste0(" -s ",softwareHere," -o ")
(deaCMD <- paste0(firstPartStable, myFNxlsx, middlePart, outdir, " -y ", outdir,"/",ymlF))

# run DEA
# system(deaCMD)
# At this point we have to switch to the CMD_DEA and customize slightly in order to have the correct RMD report

#
#
#
#
#
#

# here copy paste from CMD_DEA_custom.R
if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse", dependencies = TRUE)
}

option_list <- list(
  optparse::make_option(c("-i", "--indir"), type = "character", default = ".",
                        help = "folder containing fasta and diann-output files",
                        metavar = "path"),
  optparse::make_option(c("-d", "--dataset"), type = "character", default = "dataset.csv",
                        help = "file with annotation",
                        metavar = "character"),
  optparse::make_option(c("-y", "--yaml"), type = "character", default = "config.yaml",
                        help = "yaml configuration file",
                        metavar = "character"),
  optparse::make_option(c("-w","--workunit"), type = "character", default = NULL,
                        help = "yaml configuration file",
                        metavar = "character"),
  optparse::make_option(c("-s", "--software"), type = "character", default = NULL,
                        help = paste0("possible options: ", paste(names(prolfquapp::prolfqua_preprocess_functions), collapse = ", ")),
                        metavar = "character"),
  optparse::make_option(c("-o", "--outdir"), type = "character", default = NULL,
                        help = "output directory",
                        metavar = "character"),
  optparse::make_option(c("--libPath"), type = "character", default = NULL,
                        help = " (optional) R library path",
                        metavar = "string")
)

parser <- optparse::OptionParser(usage = "%prog config.yaml --software DIANN --indir .", option_list = option_list)

if (length(commandArgs(TRUE)) == 0) {
  optparse::print_help(parser)
  #quit(status = 1)
}

arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
ymlfile <- arguments$args

logger::log_appender(logger::appender_console)
logger::log_info("LIBRARY PATHS (.libPaths()):",paste(.libPaths(), collapse = "\n"))

prolfquapp::set_lib_path(opt$libPath);

library(prolfquapp)
logger::log_info("using : ", system.file(package = "prolfqua"))
logger::log_info("using : ", system.file(package = "prolfquapp"))

if (FALSE) {
  ymlfile <- "p35593_uniprot_paired/WholeProtUniprot.yaml"
  opt$dataset <- "p35593_uniprot_paired/dataset.xlsx"
  opt$indir <- "o35593_prot_ionquant/"
}
if (FALSE) {
  ymlfile <- "minimalPhosphoAnalysis.yaml"
  opt$dataset <- "p35540_WU313409_annotationTableXLSX.xlsx"
  opt$indir <- "Fragpipe_o35920_phospho/"
  opt$software <- "FP_multisite"
}
if (FALSE) {
  ymlfile <- "minimal.yaml"
  opt$indir <- "."
  opt$software <- "MSSTATS"
  opt$dataset <- "dataset.csv"
}
if (FALSE) {
  ymlfile <- "config.yaml"
  opt$indir <- "."
  opt$software <- "DIANN"
  opt$dataset <- "dataset.csv"
}
if (FALSE) {
  ymlfile <- "FragPipe_f20/msstats20.yaml"
  opt$indir <- "FragPipe_f20"
  opt$software <- "MSSTATS_FP_DIA_PEPTIDE"
  opt$dataset <- "FragPipe_f20/dataset_msstats20_parallel.xlsx"
  opt$workunit <- "testing_peptide"
}

# jg (2024-12-17)
# here set processing options
#bash prolfqua_dea.sh -i ../o36946_FP_TMTi_enriched/ -d p36946_vs_vs_WTctrlv2_annotationTable.xlsx -s FP_multisite -o DEA_20241213_PhosphoEnriched_p36946_WU_vs_WTctrlv2 -y DEA_20241213_PhosphoEnriched_p36946_WU_vs_WTctrlv2/minimalYaml_vsvs_WTctrlv2.yaml"
if (TRUE) {
  ymlfile <- "DEA_20241213_PhosphoEnriched_p36946_WU_vs_WTctrlv2/minimalYaml_vsvs_WTctrlv2.yaml"
  opt$indir <- "../o36946_FP_TMTi_enriched/"
  opt$software <- "BGS_DEFAULT_PROTEIN" # wouldnt matter at this point
  #opt$dataset <- "p36946_vs_vs_WTctrlv2_annotationTable.xlsx"
  opt$dataset <- myFNxlsx
  opt$workunit <- "DEA_20241218_PhosphoEnriched_CytokineEffectsOnWT"
  opt$outdir <- opt$workunit
}
if (FALSE) {
  ymlfile <- "config_2_or_more.yaml"
  opt$indir <- "FragPipe_f20_diann"
  opt$software <- "DIANN"
  opt$dataset <- "dataset_all_interaction_no_Subject.xlsx"
  opt$workunit <- "Diet_Subgroup_2PEP"
  #./prolfqua_dea.sh -s DIANN -i FragPipe_f20_diann \
  #-d dataset_all_interaction_no_Subject.xlsx \
  #-y config_2_or_more.yaml -w f20_diann_with_interaction
}

ymlfile <- if ( length(ymlfile) == 0 ) { opt$yaml } else { ymlfile }

logger::log_info("YAML file read: ", ymlfile)
stopifnot(file.exists(ymlfile))

GRP2 <- prolfquapp::get_config(ymlfile)


res <- prolfquapp::sync_opt_config(opt, GRP2)
opt <- res$opt
GRP2 <- res$config

dir.create(opt$outdir)
dir.create(GRP2$get_zipdir())

current_time <- Sys.time()

formatted_time <- format(current_time, "%Y%m%d%H%M")
logfile <- paste0("prolfqua_", formatted_time, ".log")
appender_combined <- logger::appender_tee(file.path(GRP2$get_zipdir(), logfile))
logger::log_appender(appender_combined)
logger::log_info(prolfquapp::capture_output(quote(lobstr::tree(opt))))
logger::log_info("Writing to output directory : ", GRP2$get_zipdir(), " and file :", logfile)



logger::log_info("prolfquapp paramters : ")
logger::log_info( prolfquapp::capture_output( quote(lobstr::tree(R6_extract_values(GRP2)))))


annotation <- file.path(opt$dataset) |>
  prolfquapp::read_table_data() |> prolfquapp::read_annotation(prefix = GRP2$group)

logger::log_info("Contrasts: \n", paste(annotation$contrasts, collapse = "\n"))

logger::log_info("Factors : ",paste(annotation$atable$factor_keys_depth(), collapse = "\n"))
prolfquapp::copy_DEA_Files()
logger::log_info("Software: ", opt$software)


#debug(prolfquappPTMreaders::preprocess_FP_multi_site)
#undebug(prolfquappPTMreaders::preprocess_FP_multi_site)

result <- tryCatch({
  # Attempt to run the function

  procsoft <- preprocess_software(
    opt$indir,
    annotation,
    software = opt$software,
    preprocess_functions_str = GRP2$ext_reader,
    pattern_contaminants = GRP2$processing_options$pattern_contaminants,
    pattern_decoys = GRP2$processing_options$pattern_decoys
  )
  # Return the result if successful
  list(value = procsoft, error = NULL, stack_trace = NULL)
}, error = function(e) {
  # On error, capture the stack trace as text
  stack_trace <- capture.output(traceback())
  # Return the error message and stack trace
  list(
    value = NULL,
    error = conditionMessage(e),
    stack_trace = paste(stack_trace, collapse = "\n")
  )
})

if (!is.null(result$error)) {
  logger::log_error(result$error, "\n")
  logger::log_error("Stack trace:\n")
  logger::log_error(result$stack_trace, "\n")
  stop("error occured")
} else {
  xd <- result$value$xd
  files <- result$value$files
}


logger::log_info("Processing done:", opt$software)
logger::log_info(paste(c("Protein Annotation :\n",capture.output( print(xd$protein_annotation$get_summary()))),collapse = "\n"))
logger::log_info("AGGREGATING PEPTIDE DATA: {GRP2$processing_options$aggregate}.")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)

#check -> works
xd$lfqdata$summarize_hierarchy()
xd$lfqdata$hierarchy_counts()

logger::log_info("END OF PROTEIN AGGREGATION")
logger::log_info("RUN ANALYSIS")

# important set nr_children filter to one
GRP2$processing_options$nr_peptides # it s at 1

grp <- prolfquapp::generate_DEA_reports2(
  lfqdata,
  GRP2,
  xd$protein_annotation,
  annotation$contrasts)
logger::log_info("Writing results to: " ,  GRP2$get_zipdir())


sr <- lfqdata$get_Summariser()
sr$plot_missingness_per_group_cumsum()
sr$plot_hierarchy_counts_sample()

# for sure here we should have a different markdown? _Grp2Analysis_Phospho_V2
# let's try
outdir <- prolfquapp::write_DEA_all(
  grp, name = "", boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")

lfqdataIB <- xd$lfqdata$get_subset(xd$protein_annotation$clean(
  contaminants = GRP2$processing_options$remove_cont,
  decoys = GRP2$processing_options$remove_decoys))

# do not write when peptide level analysis --> also for phospho analysis this generates errors and does not make sense
# if (length(xd$protein_annotation$pID) == 1) {
#   ibaq <- compute_IBAQ_values(lfqdataIB, xd$protein_annotation)
#   writexl::write_xlsx(
#     ibaq$to_wide()$data,
#     path = file.path(grp$get_result_dir(), paste0("IBAQ_",opt$workunit,".xlsx")))
# }

logger::log_info("Writing summarized experiment.")
SE <- prolfquapp::make_SummarizedExperiment(grp)
saveRDS(SE, file = file.path( grp$get_result_dir(), "SummarizedExperiment.rds"))

# here we specify and write out input files for the analysis
logger::log_info("Creating directory with input files :", GRP2$get_input_dir())
dir.create(GRP2$get_input_dir())

prolfquapp::copy_DEA_Files(workdir = GRP2$get_input_dir())
prolfquapp::copy_shell_script(workdir = GRP2$get_input_dir())

file.copy(c(files$data, files$fasta, ymlfile, opt$dataset), GRP2$get_input_dir())

logger::log_info("Write yaml with parameters: ", file.path(GRP2$get_input_dir(), "minimal.yaml"))

GRP2$RES <- NULL
GRP2$pop <- NULL
yaml::write_yaml(prolfquapp::R6_extract_values(GRP2), file = file.path(GRP2$get_input_dir(), "minimal.yaml"))











