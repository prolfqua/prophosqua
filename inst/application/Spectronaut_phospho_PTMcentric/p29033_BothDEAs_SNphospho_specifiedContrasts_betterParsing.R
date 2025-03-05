#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#

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

# from Communication w/ Radhika
# - the conditions in your samples
# Samples are siC, siG, siC_Iso, siG_iso (3 replicates each)
#
# siC – Control Knockdown
# siG- Grb14 Knockdown
# siC_Iso - Control Knockdown + Iso treated
# siG_iso – Grb14 Knockdown + Iso treated
#
# - the comparisons:
#   1. siG vs siC
# 2. siG_Iso vs siC_Iso

# update!!
# 1. siG_Untreated vs siC_Untreated
# 2. siG_Treated vs siC_Treated
# 3. siG_Treated vs siG_Untreated

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
# go for phospho peptide part
# this is a next version -> here we use the site centric report! Each single site is individually reported!
#
#
WUID  <- "enriched"

#get Spectronaut SITE export from Antje in
fullSiteReport <- read_tsv("../20250304_p29033_o37804_phospho_DIA_wLib_mouse_correctDB_Report_Site.tsv")

# build dataset from report -> here we use self specified contrasts
# create Dataset as before to run custom DEA in prolfqua base
dsf <- fullSiteReport |> select(R.FileName, R.Condition) |> distinct()
# we further need Grouping Var column and CONTROL
dsf$`Grouping Var` <- dsf$R.Condition
dsf$R.Condition <- NULL
dsf$raw.file <- dsf$R.FileName
dsf$SampleName <- gsub(x = dsf$R.FileName, pattern = "\\d+_\\d+_S\\d+_", replacement = "")
dsf$SampleName <- gsub(x=dsf$SampleName, pattern = "_phoshoDIA", replacement = "")
dsf$R.FileName <- NULL
dsf

# watch out here .. conditions are named differently than before:
# "siG_Treated"   "siG_Untreated" "siC_Treated"   "siC_Untreated"
dsf$`Grouping Var` <- gsub(x = dsf$`Grouping Var`, pattern = "_I", replacement = "_Treated")
dsf$`Grouping Var` <- gsub(x = dsf$`Grouping Var`, pattern = "_C", replacement = "_Untreated")
# do we have the same conditions in both datasets?
dsf$`Grouping Var` %in% dsf_tot$`Grouping Var`


# double check
table(dsf$`Grouping Var`)

# extend self specified contrasts from dataset from total
dsf_tot
dsf$ContrastName <- dsf_tot$ContrastName
dsf$Contrast <- dsf_tot$Contrast
dsf

# write out dsf to excel for prolfquapp
(xlsFN <- paste0(fgczProject,"_", WUID,"_",OIDfgcz,".xlsx"))
write.xlsx(dsf, xlsFN)

# get annotation object -> this already contains all potential contrasts for specified control
annotation <- read_annotation(dsf, prefix = "G_")

# working on prot annot
path = "../20250304_p29033_o37804_phospho_DIA_wLib_mouse_correctDB_Report/"
dir(path)

# get psm and fasta file
files <- get_BGS_files(path = path, bgs_pattern = "Experiment1_Report_BGS")
files$fasta
files$data # we are not gonna work with this report but rather use the PTM site centric report

annot <- annotation$annot
atable <- annotation$atable
annot <- dplyr::mutate(annot, raw.file = gsub("^x|.d.zip$|.raw$",
                                              "", (basename(annot[[atable$fileName]]))))

#check max qvalue on ProteinGroup -> should be below 0.01 by default
fullSiteReport$PG.Qvalue |> max()

#filter BGS report for
# 1) keep only phospho
unique(fullSiteReport$PTM.ModificationTitle)
siteReport <- fullSiteReport |> filter(grepl(PTMtitle,x = PTM.ModificationTitle))

# 2) BestLocalizationProbability > 0.75
siteReport <- siteReport |> filter(PTM.SiteProbability > best_loc_prob_threshold)

# 3) filter on minimum threshold
siteReport <- siteReport |> filter(PTM.Quantity > abundance_threshold)

# 4) reshape input before going into prolfqua
siteReport4proflqua <- siteReport |> select(R.FileName, PG.ProteinAccessions, PTM.CollapseKey, PTM.Multiplicity, PTM.NrOfCollapsedPeptides, PTM.SiteAA, PTM.SiteLocation, PTM.Quantity)
siteReport4proflqua$Probability <- 1
siteReport4proflqua$qValue <- 0.01

# Prolfqua Configuration
atable$ident_Score = "Probability"
atable$ident_qValue = "qValue"
atable$hierarchy[["protein_Id"]] <- c("PG.ProteinAccessions")
atable$hierarchy[["site"]] <- c("PG.ProteinAccessions","PTM.CollapseKey")
atable$hierarchy[["mod_peptide_Id"]] <- c("PTM.CollapseKey")
atable$set_response("PTM.Quantity")
atable$hierarchyDepth <- 2
atable$get_response()

#drop contrast column from sf before joining
myAnnoForJoining <- annot
myAnnoForJoining$ContrastName <- NULL
myAnnoForJoining$Contrast <- NULL

# join
psma <- dplyr::inner_join(myAnnoForJoining, siteReport4proflqua, join_by("raw.file" == "R.FileName"))

# Analysis Configuration
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(psma, config)
lfqdata <- prolfqua::LFQData$new(adata, config)

# Fasta file
fasta_file <- files$fasta

fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
# reshape fasta annot for ProteinAnnotation
fasta_annot <- dplyr::rename(fasta_annot, description = fasta.header)
fasta_annot <- dplyr::rename(fasta_annot, protein_Id = proteinname)

# for protein annotation it is important to have Depth 1 (protein level) we want to annotation the proteins w/ descriptions
lfqdata$config$table$hierarchyDepth <- 1

#ProteinAnnotation$debug("initialize")
prot_annot <- prolfquapp::ProteinAnnotation$new(
  lfqdata , fasta_annot,
  description = "description",
  cleaned_ids = "protein_Id", #this has to match one of the columns in fasta_annot that is the same as protein_Id in
  full_id = "fasta.id",
  exp_nr_children = "nr_children",
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
GRP2_phos <- prolfquapp::make_DEA_config_R6(PATH = myPath,WORKUNITID = WUID,
                                            PROJECTID = fgczProject, ORDERID = OIDfgcz)

# important for phos
GRP2_phos$pop$nr_peptdes <- 1
GRP2_phos$pop$aggregate <- "sum_topN"
GRP2_phos$pop$contrasts <- annotation$contrasts

logger::log_info("GENERATING DEA REPORTS")
logger::log_info("starting modelling")

# fix for old grps
GRP2_phos$pop$aggregate
GRP2_phos$processing_options$transform <- "robscale"
GRP2_phos$processing_options$aggregate <- "topN"
GRP2_phos$pop$transform <- GRP2_phos$processing_options$transform
GRP2_phos$pop$Diffthreshold <- GRP2_phos$processing_options$diff_threshold
GRP2_phos$pop$FDRthreshold <- GRP2_phos$processing_options$FDR_threshold

logger::log_info("AGGREGATING PEPTIDE DATA!")
# nothing to aggregate here since site centric export is used!! (has been aggregated by Spectronaut)

lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method =
                                        GRP2_phos$processing_options$aggregate,N = 1000)
logger::log_info("data aggregated: {GRP2_phos$processing_options$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

lfqdata$hierarchy_counts()

# grp_phos
grp_phos <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2_phos, prot_annot, GRP2_phos$pop$contrasts)

# all fine and working as expected? do we see some significant hits?
grp_phos$RES$lfqData$to_wide()
grp_phos$RES$contrastsData_signif

logger::log_info("DONE WITH DEA REPORTS")
# result dir
length(unique(grp_phos$RES$contrastsData$site))

(GRP2_phos$zipdir_name <- paste0("/DEA_",datetoday, "_", fgczProject, "_", WUID, "_",OIDfgcz))
GRP2_phos$get_zipdir()
dir.create(GRP2_phos$get_zipdir())

# get markdown files here
# also copy the phos specific Rmd files from prophosqua
source("~/GitHub/prophosqua/R/prophosqua_copy_helpers.R")
prolfquapp::copy_DEA_Files() # for bib file
# copy_phosDEA_Spectronaut()

# make sure some info in the reports is correct
grp_phos$software <- "Spectronaut-GUI"

GRP2 <- grp_phos
# writing reports
descri <- "CustomPhosphoAnalysis"
print(descri)

# both htmls are rendered QC and Phospho
outpath <- prolfquapp::write_DEA_all(grp_phos, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")

# Save RData from enriched and total (only lfqdata is overwritten?) # keep lfqdata, grp, adata separate for phos and total!
(imageFN <- paste(fgczProject, descri, "BothDEA","Analyzed.RData", sep="_"))
save.image(imageFN)

