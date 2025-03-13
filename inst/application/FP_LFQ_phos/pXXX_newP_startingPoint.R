#!/usr/local/bin/Rscript
#
# 2025
#
# This script is used to run the R script for the analysis of the data


library(prolfquapp)
library(prophosqua)
library(readr)
library(dplyr)
library(openxlsx)
library(stringr)

# here we use Spectronaut "site centric" export from the PTM workflow
# Rmd are customized for phospho analysis and can be found in prophosqua (https://github.com/prolfqua/prophosqua/)

# variables
fgczProject <- "p29033"
WUID  <- "TotalProteome"
OIDfgcz <- "betterParsed"
datetoday <- format(Sys.Date(), "%Y%m%d")
best_loc_prob_threshold <- 0.75
global_qvalue_threshold <- 0.01
abundance_threshold <- 1
PTMtitle <- "Phospho"



pattern_decoys <- "rev_" # there are no decoys in this db?
pattern_contaminants <- "contam_sp"

# default BGS report for total
BGSreport_total <- read_tsv("../20250304_p29033_o37804_proteome_directDIA_mouse_correctDB_Report/Experiment1_Report_BGS Factory Report (Normal).tsv")

# create Dataset as before to run DEA
dsf_tot <- BGSreport_total |> select(R.FileName, R.Condition) |> distinct()
# we further need Grouping Var column and CONTROL
dsf_tot$R.Condition <- NULL
dsf_tot$raw.file <- dsf_tot$R.FileName
dsf_tot$SampleName <- gsub(x = dsf_tot$R.FileName, pattern = "\\d+_\\d+_S\\d+_", replacement = "")
dsf_tot$SampleName <- gsub(x = dsf_tot$SampleName, pattern = "_proteomeDIA", replacement = "")
dsf_tot$`Grouping Var` <- gsub(x = dsf_tot$SampleName, pattern = "_[I|C]\\d_", replacement = "_")
dsf_tot$R.FileName <- NULL
dsf_tot

unique(dsf_tot$`Grouping Var`)

#dsf_tot$CONTROL <- "T"
#dsf_tot$CONTROL[dsf_tot$`Grouping Var` == "siG_Untreated"] <- "C"

#unique(dsf_tot$`Grouping Var`)
#[1] "siG_Treated"   "siG_Untreated" "siC_Treated"   "siC_Untreated"

# selfspecifiedContrasts
ContrastName <- c("siG_vs_siC","siGIso_vs_siCIso", "siGIso_vs_siG")
nCont <- length(ContrastName)
ContrastName <- c(ContrastName, rep(NA, nrow(dsf_tot) - nCont))
Contrast <- c("G_siG_Untreated - G_siC_Untreated", "G_siG_Treated - G_siC_Treated", "G_siG_Treated - G_siG_Untreated")
Contrast <- c(Contrast, rep(NA, nrow(dsf_tot) - nCont))
dsf_tot$ContrastName <- ContrastName
dsf_tot$Contrast <- Contrast
dsf_tot

# write out dsf to excel for prolfquapp
(xlsFN <- paste0(fgczProject,"_", WUID,"_",OIDfgcz,".xlsx"))
write.xlsx(dsf_tot, xlsFN)

# go for prolfquapp
# copy shell scripts -> this only is need the first time
# system('R --vanilla -e "prolfquapp::copy_shell_script(workdir = \'.\')"')

# make executable
# system("chmod a+x prolfqua_*")

# BGS_DEFAULT_PROTEIN
# yaml to be generated with options
#GRP2_phos$zipdir_name <- paste0("/DEA_",datetoday, "_", fgczProject, "_", WUID, "_phos/")
(outdir <- paste0("DEA", "_",fgczProject,"_",datetoday,"_",WUID,"_",OIDfgcz))
(ymlF <- paste0("minimalYaml_robscale",".yaml"))
mkDirCMD <- paste("mkdir", outdir)
system(mkDirCMD)
softwareHere <- "BGS_DEFAULT_PROTEIN"
(ymlCMD <- paste0("bash prolfqua_yaml.sh --norm robscale --outdir ",outdir," --workunit ", WUID," -p ",fgczProject," -O ",OIDfgcz," -s ", softwareHere, " --yaml ",ymlF))
system(ymlCMD)

# run DEA
myInputFolder <- "../20250304_p29033_o37804_proteome_directDIA_mouse_correctDB_Report/"
firstPartStable <- paste0("bash prolfqua_dea.sh -i ",myInputFolder ," -d ")
middlePart <- paste0(" -s ",softwareHere," -o ")
(deaCMD <- paste0(firstPartStable, xlsFN, middlePart, outdir, " -y ", outdir,"/",ymlF))

# run DEA
# important tsv file has to have particular filename
# Experiment1_Report_BGS Factory Report (Normal)
system(deaCMD)


#
#
#
#     Phospho Enriched Part
#
#
#




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


