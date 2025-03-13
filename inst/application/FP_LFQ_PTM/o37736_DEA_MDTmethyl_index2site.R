#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#

# here we use FragPipe "site centric" export from the LFQ workflow
# Rmd are customized for methylation analysis and can be found in prophosqua (https://github.com/prolfqua/prophosqua/)
message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquapp Version :", packageVersion("prolfquapp"), "\n")

library(prolfqua)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(openxlsx)
library(prophosqua)
library(prolfquapp)

# some fgcz specific settings
# params ideally taken from yaml
fgczProject <- "p37736"
OIDfgcz <- "o37736"
descri <- "MultiMethyl"
fracti <- "enriched_onlyToxoP"
WUID <- "WU304087"
datetoday <- format(Sys.Date(), "%Y%m%d")

# get LocalizationSiteProbThreshold
LocProbThresh <- 0.75
pattern_decoys <- "rev_"
pattern_contaminants <-"Cont_"

# do prolfqua on FragPipe methylated peptides! Search is done w/ 3 varmods
# -> site specific files are reported individually
###

# since Toxopl g. + host (Human) samples but only ToxoG is intersting, we do filtering before reading in


# source("/Users/jonasgrossmann/GitHub/prophosqua/R/FP_phosphoHelperFunctions_v3_202310.R")
mono_methyl <- "../o37736_TTF_FPtables_MaxQuantLFQ/combined_site_KR_14.0156.tsv"
di_methyl <- "../o37736_TTF_FPtables_MaxQuantLFQ/combined_site_KR_28.0313.tsv"
tri_methyl <- "../o37736_TTF_FPtables_MaxQuantLFQ/combined_site_K_42.0470.tsv"


peptide_14 <- tidy_FragPipe_phospho_STYfile(styphosphopeptideTable = mono_methyl, locProbThreshold = LocProbThresh)
peptide_28 <- tidy_FragPipe_phospho_STYfile(styphosphopeptideTable = di_methyl, locProbThreshold = LocProbThresh)
peptide_42 <- tidy_FragPipe_phospho_STYfile(styphosphopeptideTable = tri_methyl, locProbThreshold = LocProbThresh)

# add methylmass
peptide_14$index <- paste0(peptide_14$index,"~Mono")
peptide_28$index <- paste0(peptide_28$index,"~Di")
peptide_42$index <- paste0(peptide_42$index,"~Tri")

# combine
ptm_peptides <- rbind(peptide_14, peptide_28, peptide_42)
head(ptm_peptides)

# remove Homo sapiens here?
# multiSite_long <- multiSite_long |> filter(grepl("_TOXGM", x = fasta.id))
ptm_peptides <- ptm_peptides |> filter(grepl("_TOXG", x = protein))

# build up annotatation from raw file names
(allfiles <- unique(ptm_peptides$raw.file))
table(gsub(sapply(strsplit(allfiles, split = "_"), function(x)x[4]), pattern = "\\d$", replacement = ""))

rm(myAnnot)
# SampleName
table(sapply(strsplit(allfiles, split = "_"), function(x)x[4]))

# initalize myAnnot
myAnnot <- data.frame(raw.file = unique(ptm_peptides$raw.file))
myAnnot$SampleName <- sapply(strsplit(allfiles, split = "_"), function(x)x[4])
myAnnot$'Grouping.Var' <- gsub(sapply(strsplit(myAnnot$raw.file, split = "_"), function(x)x[4]), pattern = "\\d$", replacement = "")
table(myAnnot$'Grouping.Var')

# relabel grouping var
# https://fgcz-bfabric.uzh.ch/bfabric/order/show.html?id=37736&tab=details
myAnnot$'Grouping.Var'[myAnnot$'Grouping.Var' == "TIR"] <- "TIR1" # Control
myAnnot$'Grouping.Var'[myAnnot$'Grouping.Var' == "API"] <- "AKMTPCKMTIMC30"
myAnnot$'Grouping.Var'[myAnnot$'Grouping.Var' == "AP"] <- "AKMTPCKMT"
myAnnot$'Grouping.Var'[myAnnot$'Grouping.Var' == "A"] <- "AKMT"
myAnnot$'Grouping.Var'[myAnnot$'Grouping.Var' == "I"] <- "IMC30"
myAnnot$'Grouping.Var'[myAnnot$'Grouping.Var' == "P"] <- "PCKMT"
myAnnot$'Grouping.Var'[myAnnot$'Grouping.Var' == "AI"] <- "AKMTIMC30"

table(myAnnot$'Grouping.Var')

# table(myAnnot$'Grouping.Var')
# myAnnot$CONTROL <- "T"
#
# # what is control condition
# myAnnot$CONTROL[myAnnot$'Grouping.Var' == "Control"] <- "C"
# table(myAnnot$CONTROL)

# Here’s a summary of the groups we would like to compare :
# Each of the groups (singles, doubles, triple) vs TIR1
# AKMT vs AKMT IMC30
# AKMT vs AKMT PCKMT
# AKMT IMC30 vs AKMT PCKMT IMC30
# AKMT PCKMT vs AKMT PCKMT IMC30

# selfspecifiedContrasts
ContrastName <- c("AKMT_vs_TIR1","AKMTIMC30_vs_TIR1","AKMTPCKMT_vs_TIR1","AKMTPCKMTIMC30_vs_TIR1","IMC30_vs_TIR1","PCKMT_vs_TIR1","AKMT_vs_AKMTIMC30","AKMT_vs_AKMTPCKMT","AKMTIMC30_vs_AKMTPCKMTIMC30","AKMTPCKMT_vs_AKMTPCKMTIMC30")
nCont <- length(ContrastName)
ContrastName <- c(ContrastName, rep(NA, nrow(myAnnot) - nCont))
Contrast <- c("G_AKMT - G_TIR1","G_AKMTIMC30 - G_TIR1","G_AKMTPCKMT - G_TIR1","G_AKMTPCKMTIMC30 - G_TIR1","G_IMC30 - G_TIR1","G_PCKMT - G_TIR1","G_AKMT - G_AKMTIMC30","G_AKMT - G_AKMTPCKMT","G_AKMTIMC30 - G_AKMTPCKMTIMC30","G_AKMTPCKMT - G_AKMTPCKMTIMC30")
Contrast <- c(Contrast, rep(NA, nrow(myAnnot) - nCont))
myAnnot$ContrastName <- ContrastName
myAnnot$Contrast <- Contrast
myAnnot

# write out dsf to excel for prolfquapp
(xlsFN <- paste0(fgczProject,"_", WUID,"_",OIDfgcz,"Dataset.xlsx"))
write.xlsx(myAnnot, xlsFN)

# get annotation object -> this already contains all potential contrasts for specified control
annotation <- read_annotation(myAnnot, prefix = "G_")
annotation$contrasts

# working on prot annot
path = "../o37736_TTF_FPtables_MaxQuantLFQ/"
dir(path)

# get psm and fasta file
files <- get_BGS_files(path = path, bgs_pattern = "Experiment1_Report_BGS")
files$fasta
files$data # we are not gonna work with this report but rather use the PTM site centric report

annot <- annotation$annot
atable <- annotation$atable
annot <- dplyr::mutate(annot, raw.file = gsub("^x|.d.zip$|.raw$",
                                              "", (basename(annot[[atable$fileName]]))))

#filter BGS report for

# 2) BestLocalizationProbability > 0.75 # this is already handled upon reading the file
# siteReport <- siteReport |> filter(PTM.SiteProbability > best_loc_prob_threshold)

str(ptm_peptides)

# reshape column names
ptm_peptides <- dplyr::rename(ptm_peptides,  protein_Id = proteinid)
ptm_peptides <- dplyr::rename(ptm_peptides,  site = index)

# 4) reshape input before going into prolfqua
siteReport4proflqua <- ptm_peptides |> select(site, raw.file, protein_Id, protein, gene, intensity, peptide)
siteReport4proflqua$Probability <- 1
siteReport4proflqua$qValue <- 0.01
colnames(siteReport4proflqua)

# Prolfqua Configuration
atable$ident_Score = "Probability"
atable$ident_qValue = "qValue"
atable$hierarchy[["protein_Id"]] <- c("protein_Id")
atable$hierarchy[["site"]] <- c("protein_Id","site")
#atable$hierarchy[["mod_peptide_Id"]] <- c("peptide")
atable$set_response("intensity")
atable$hierarchyDepth <- 2
atable$get_response()

#drop contrast column from sf before joining
myAnnoForJoining <- annot
myAnnoForJoining$ContrastName <- NULL
myAnnoForJoining$Contrast <- NULL

# join
psma <- dplyr::inner_join(myAnnoForJoining, siteReport4proflqua, join_by("raw.file" == "raw.file"))
colnames(psma)

# Analysis Configuration
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(psma, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
colnames(lfqdata$data)

# Fasta file
fasta_file <- files$fasta

fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
# reshape fasta annot for ProteinAnnotation
fasta_annot <- dplyr::rename(fasta_annot, description = fasta.header)
fasta_annot <- dplyr::rename(fasta_annot, protein_Id = proteinname)

#colnames(psma)
#psma <- dplyr::rename(psma, protein_Id = proteinid)

# for protein annotation it is important to have Depth 1 (protein level) we want to annotation the proteins w/ descriptions
lfqdata$config$table$hierarchyDepth <- 1

colnames(fasta_annot)
colnames(lfqdata$data)

#ProteinAnnotation$debug("initialize")
prot_annot <- prolfquapp::ProteinAnnotation$new(
  lfqdata , fasta_annot,
  description = "description",
  cleaned_ids = "protein_Id", #this has to match one of the columns in fasta_annot that is the same as protein_Id in
  full_id = "fasta.id",
  exp_nr_children = "nr_tryptic_peptides", # not sure if nr_children should be used?
  pattern_contaminants = pattern_contaminants,
  pattern_decoys = pattern_decoys
)

# prolfqua specific
lfqdata$remove_small_intensities()

lfqdata$data #checks if all is here
lfqdata$response() #checks if all is here

lfqdata$config$table$hierarchyDepth <- 2

#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata$config$table$hierarchy_keys()

# gen Grp2
myPath <- getwd()
GRP2_ptm <- prolfquapp::make_DEA_config_R6(PATH = myPath,WORKUNITID = WUID,
                                            PROJECTID = fgczProject, ORDERID = OIDfgcz)

# important for phos
GRP2_ptm$pop$nr_peptdes <- 1
GRP2_ptm$pop$aggregate <- "sum_topN"
GRP2_ptm$pop$contrasts <- annotation$contrasts

logger::log_info("GENERATING DEA REPORTS")
logger::log_info("starting modelling")

# fix for old grps
GRP2_ptm$pop$aggregate
GRP2_ptm$processing_options$transform <- "robscale"
GRP2_ptm$processing_options$aggregate <- "topN"
GRP2_ptm$pop$transform <- GRP2_ptm$processing_options$transform
GRP2_ptm$pop$Diffthreshold <- GRP2_ptm$processing_options$diff_threshold
GRP2_ptm$pop$FDRthreshold <- GRP2_ptm$processing_options$FDR_threshold

logger::log_info("AGGREGATING PEPTIDE DATA!")
# nothing to aggregate here since site centric export is used!! (has been aggregated by Spectronaut)
colnames(lfqdata$data)
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method =
                                        GRP2_ptm$processing_options$aggregate,N = 1000)
logger::log_info("data aggregated: {GRP2_ptm$processing_options$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

lfqdata$hierarchy_counts()

# grp_ptm
grp_ptm <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2_ptm, prot_annot, GRP2_ptm$pop$contrasts)

# all fine and working as expected? do we see some significant hits?
grp_ptm$RES$lfqData$to_wide()
grp_ptm$RES$contrastsData_signif

logger::log_info("DONE WITH DEA REPORTS")
# result dir
length(unique(grp_ptm$RES$contrastsData$site))

(GRP2_ptm$zipdir_name <- paste0("/DEA_",datetoday, "_", fgczProject, "_", WUID, "_",OIDfgcz))
GRP2_ptm$get_zipdir()
dir.create(GRP2_ptm$get_zipdir())

# get markdown files here
# also copy the ptm specific Rmd files from prophosqua
source("~/GitHub/prophosqua/R/prophosqua_copy_helpers.R")
prolfquapp::copy_DEA_Files() # for bib file
# copy_phosDEA_Spectronaut()

# make sure some info in the reports is correct
grp_ptm$software <- "FragPipe-GUI"

GRP2 <- grp_ptm
# writing reports
descri <- "CustomizedMethylationSiteAnalysis"
print(descri)

# both htmls are rendered QC and Phospho
#outpath <- prolfquapp::write_DEA_all(grp_ptm, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")
outpath <- prolfquapp::write_DEA_all(grp_ptm, boxplot = FALSE, markdown = "_Grp2Analysis_PTMsite.Rmd")

# Save RData from enriched and total (only lfqdata is overwritten?) # keep lfqdata, grp, adata separate for phos and total!
(imageFN <- paste(fgczProject, descri, "MethylationDEA","Analyzed.RData", sep="_"))
save.image(imageFN)


# testing some
#grp2$RES$transformedlfqData$get_Plotter()

pl <- grp_ptm$RES$transformedlfqData$get_Plotter()
pl$heatmap()
