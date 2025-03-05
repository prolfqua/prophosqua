# 2023-10-19: adapting for phosphoTMT total proteome and enriched
# 2024-03-26: running with latest prolfqua (1.1.5) and prlfquapp (0.1.9) before going prophosqua
# 2024-04-10: next version of DEA w/ prolfquapp 0.1.9
# 2024-04-25: towards clean version in prophosqua

# data and annotation based  on: https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details


# This script
# 1) starts from psm.tsv for TotalProteome
# 2) from multisite_abundance_none.tsv (from TMT-integrator) for the phospho-enriched part


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
source("/Users/jonasgrossmann/GitHub/prophosqua/R/FP_phosphoHelperFunctions_v3_202310.R")

# Integration of 2-plex phospho-TMT data
# how does the multi-site export look like -> for two plexes?

# variables
fgczProject <- "p37688"
WUID  <- "specifiedContrasts"
OIDfgcz <- "o37688"


path <- "."
# v3
#prolfquapp::copy_DEA_Files()
# also copy the phospho specific Rmd files from prophosqua
#prophosqua::copy_phosphoDEA_FragPipe_TMT()


# Look at annotation file provided and adapt it to the needs of prolfquapp
dsf <- "../proteome_FP22_prot/experiment_annotation.tsv"

# make minimal dsf
dsf <- readr::read_tsv(dsf)
dsf
dsf$channel <- dsf$sample
dsf$sample
dsf$plex <- NULL
dsf$condition <- NULL
dsf$replicate <- NULL
dsf$sample_name <- NULL

# we further need CONTROL column and Grouping Var
# Grouping Var
dsf$`Grouping Var` <- gsub(x = dsf$sample, pattern = "(.*)_\\d+", replacement = "\\1")
table(dsf$`Grouping Var`)

# specify Contrasts
# from paper
#Nine total comparisons were made, namely:
# 1) KO Early-WT Early
# 2) KO Late-WT Late,
# 3) KO Uninfected-WT Uninfected,
# 4) KO Early-KO Uninfected
# 5) KO Late-KO Uninfected
# 6) WT Early-WT Un-infected
# 7) WT Late-WT Uninfected
# 8) Infected-Uninfected
# 9) KO-WT.
# Since the dataset was a biological investigation, the true positive
# modi cations were unknown.

unique(dsf$`Grouping Var`)
# selfspecifiedContrasts
ContrastName <- c("KOearly_vs_WTearly","KOlate_vs_WTLate", "KOuninf_vs_WTunif", "KOearly_vs_KOuninf", "KOlate_vs_KOuninf", "WTearly_vs_WTuninf", "WTlate_vs_WTuninf", "Infected_vs_Uninfected", "KO_vs_WT")
nCont <- length(ContrastName)
ContrastName <- c(ContrastName, rep(NA, nrow(dsf) - nCont))
Contrast <- c("G_KO_Early - G_WT_Early", "G_KO_Late - G_WT_Late", "G_KO_Uninfect - G_WT_Uninfect", "G_KO_Early - G_KO_Uninfect", "G_KO_Late - G_KO_Uninfect", "G_WT_Early - G_WT_Uninfect", "G_WT_Late - G_WT_Uninfect", "(G_KO_Early + G_KO_Late + G_WT_Early + G_WT_Late)/4 - (G_KO_Uninfect + G_WT_Uninfect)/2", "(G_KO_Early + G_KO_Late + G_KO_Uninfect)/3 - (G_WT_Early + G_WT_Late + G_WT_Uninfect)/3")
Contrast <- c(Contrast, rep(NA, nrow(dsf) - nCont))
dsf$ContrastName <- ContrastName
dsf$Contrast <- Contrast

# get annotation objectt
annotation <- read_annotation(dsf, prefix = "G_")

# write out annotation table
(fN <- paste0(fgczProject,"_",WUID,"_annotationTable.txt"))
write_tsv(x = annotation$annot, file = fN)

# also write to xlsx for generating yaml file afterwards
(fNxls <- paste0(fgczProject,"_",WUID,"_annotationTable.xlsx"))
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
psmF1 <- "../proteome_FP22_prot/o37688_P1_global/psm.tsv"
psmF2 <- "../proteome_FP22_prot/o37688_P2_global/psm.tsv"

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
fastaf <- "../phospho_tmt-report/p37688_db3_MusNShigella_20250219.fasta"

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

# summarize PSMS to peptideforms using TopN=1000 -> Sum all psms from same peptide
logger::log_info("AGGREGATING PSM to peptidoforms w/ Top1000!")
ag <- lfqdata$get_Aggregator()
ag$sum_topN(N = 1000)

# write aggregated out to lfqdata to continue!
lfqdata <- ag$lfq_agg
lfqdata$data
lfqdata$response()
lfqdata$hierarchy_counts()

# Roll up to protein usingusing medpolish
lfqdata$config$table$hierarchyDepth <- 1
logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method =
                                        GRP2$processing_options$aggregate)
logger::log_info("data aggregated: {GRP2$processing_options$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

# how many proteins are identified
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
WUID <- fracti
# read in data file
multiSite_file <- dir(path = "../phospho_tmt-report/", pattern = "abundance_multi-site_None.tsv", recursive = TRUE, full.names = TRUE)
xx <- readr::read_tsv(multiSite_file)

# assuming that ReferenceIntensity is ALWAYS the first column before the individual quant channels are reported
(quant_idx_start <- grep(pattern = "ReferenceIntensity", x = colnames(xx))  + 1)

multiSiteAnnot <- xx[,1:(quant_idx_start-1)]
multiSiteQuant <- xx[,c(1,quant_idx_start:ncol(xx))]

# go long
multiSite_long <- tidyr::pivot_longer(data = multiSiteQuant, cols = 2:ncol(multiSiteQuant), values_to = "abundance", names_to = "sample")
multiSite_long <- dplyr::inner_join(x = multiSiteAnnot, y = multiSite_long)

# check
#table(annotation$annot$CONTROL, annotation$annot$Grouping.Var)
table(annotation$annot$sample)
unique(multiSite_long$sample) %in% annotation$annot$sample
multiSite_long$SampleName <- multiSite_long$sample

# in case of multi-site all the filtering has been done by FragPipe
path = "../phospho_tmt-report/"
dir(path)
# get psm and fasta file
files <- get_FP_PSM_files(path = path)
files$data # can be empty here.. is looking for psm.tsv we are using multisite here
files$fasta # we only need this later
# since TMT annotation is the same maybe the join has to be adapted!

# already done!
(myPhosPath <- paste0(fgczProject,"_",fracti,"_",WUID))
GRP2_phos <- prolfquapp::make_DEA_config_R6(PATH = myPhosPath,WORKUNITID = fracti, PROJECTID = fgczProject,
                                            ORDERID = OIDfgcz)

# join with anno again this should work now with Name # if not all samples are used in the dataset they would be removed here (to be tested)
multiSite_long <- dplyr::inner_join(x = annotation$annot, y = multiSite_long)

# all fine with our conditions
unique(multiSite_long$Grouping.Var)

# potentially filter out here files not  in the game (specified in CONTROL == NA)
# multiSite_long <- multiSite_long |> filter(!is.na(CONTROL))
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

# add here phospho relevant columns w/ parsing site?

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
                                                   more_columns = c("fasta.id", "protein_length"))

prot_annot_phos$row_annot
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
grp_phos$RES$lfqData$to_wide()
grp_phos$RES$contrastsData_signif
myResPlotter <- grp_phos$RES$contrMerged$get_Plotter()
myResPlotter$volcano()

logger::log_info("DONE WITH DEA REPORTS")

# result dir
#GRP2_phos$zipdir_name <- paste0("DEA_",myPhosPath)
GRP2_phos$path
GRP2_phos$get_zipdir()
dir.create(GRP2_phos$path)
dir.create(GRP2_phos$get_zipdir())

outpath <- prolfquapp::write_DEA_all(grp_phos, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")

# Save RData from enriched and total (only lfqdata is overwritten?) # keep lfqdata, grp, adata separate for phos and total!
(imageFN <- paste(fgczProject, "DEA_total_and_enriched.RData", sep="_"))
save.image(imageFN)
