#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

library(prolfquapp)
library(prophosqua)
library(readr)
library(dplyr)
library(openxlsx)

# source FP ubi helper Functions to use with preprocess multiplex

# here we use Spectronaut "peptide centric" export
# Rmd are customized for UBI. We use the same Rmd as for Phos

# variables
fgczProject <- "p37127"
WUID  <- "specifiedContrasts"
OIDfgcz <- "o37127"
datetoday <- format(Sys.Date(), "%Y%m%d")
best_loc_prob_threshold <- 0.75
global_qvalue_threshold <- 0.01
abundance_threshold <- 1

pattern_decoys <- "REV_" # there are no decoys in this db?
pattern_contaminants <- "Cont"



# go for total
# https://github.com/prolfqua/prolfquapp

# still first the annotation has to be there
BGSreport_total <- read_tsv("../proteome/Experiment1_Report_BGS Factory Report (Normal).tsv")

# create Dataset as before to run DEA
dsf_tot <- BGSreport_total |> select(R.FileName, R.Condition) |> distinct()
# we further need Grouping Var column and CONTROL
dsf_tot$`Grouping Var` <- dsf_tot$R.Condition
dsf_tot$R.Condition <- NULL
dsf_tot$SampleName <- dsf_tot$R.FileName
dsf_tot$R.FileName <- NULL
dsf_tot$raw.file <- dsf_tot$SampleName
dsf_tot

# selfspecifiedContrasts
ContrastName <- c("UBE3BvsWT_D0", "UBE3BvsWT_D3","D3vsD0_WT", "D3vsD0_UBE3B", rep(NA, nrow(dsf_tot) - 4))
Contrast <- c("G_UBE3B_D0 - G_WT_D0", "G_UBE3B_D3 - G_WT_D3","G_WT_D3 - G_WT_D0", "G_UBE3B_D3 - G_UBE3B_D0", rep(NA, nrow(dsf_tot) - 4))

dsf_tot$ContrastName <- ContrastName
dsf_tot$Contrast <- Contrast
dsf_tot

# write out dsf to excel for prolfquapp
(xlsFN <- paste0(fgczProject,"_", WUID,"_dsf_Total.xlsx"))
write.xlsx(dsf_tot, xlsFN)

# go for prolfquapp
# copy shell scripts
system('R --vanilla -e "prolfquapp::copy_shell_script(workdir = \'.\')"')
# make executable
system("chmod a+x prolfqua_*")

# BGS_DEFAULT_PROTEIN
# dataset we already have
# yaml to be generated with options

(outdir <- paste0("DEA_",fgczProject,"_",WUID,"_totalProtein"))
(ymlF <- paste0("minimalYaml_robscale",".yaml"))
mkDirCMD <- paste("mkdir", outdir)
system(mkDirCMD)
softwareHere <- "BGS_DEFAULT_PROTEIN"
(ymlCMD <- paste0("bash prolfqua_yaml.sh --norm robscale --outdir ",outdir," --workunit ", WUID," -p ",fgczProject," -O ",OIDfgcz," -s ", softwareHere, " --yaml ",ymlF))
system(ymlCMD)

# run DEA
myInputFolder <- "../proteome/"
firstPartStable <- paste0("bash prolfqua_dea.sh -i ",myInputFolder ," -d ")
middlePart <- paste0(" -s ",softwareHere," -o ")
(deaCMD <- paste0(firstPartStable, xlsFN, middlePart, outdir, " -y ", outdir,"/",ymlF))

# run DEA
# important tsv file has to have particular filename
# Experiment1_Report_BGS Factory Report (Normal)
system(deaCMD)





#get Spectronaut peptide export from Antje in
fullBGSreport <- read_tsv("../enriched/20241212_142604_20241212_o37127_diGly_enriched_PTMReport.tsv")


# build dataset from report -> here we use self specified contrasts
# unique(fullBGSreport$R.FileName)
# allCond <- unique(fullBGSreport$R.Condition)
# allCond
# (controlCondition <- allCond[1])

# create Dataset as before to run DEA
dsf <- fullBGSreport |> select(R.FileName, R.Condition) |> distinct()
# we further need Grouping Var column and CONTROL
dsf$`Grouping Var` <- dsf$R.Condition
dsf$R.Condition <- NULL
#dsf$CONTROL <- "T"
#dsf$CONTROL[dsf$`Grouping Var` == controlCondition] <- "C"
dsf$SampleName <- dsf$R.FileName
dsf$R.FileName <- NULL
dsf$raw.file <- dsf$SampleName
dsf

# double check
table(dsf$`Grouping Var`)
#table(dsf$CONTROL, dsf$`Grouping Var`)

# add selfspecified contrasts
#dsf$CONTROL <- NULL

# do selfspecified contrasts as above
dsf$ContrastName <- ContrastName
dsf$Contrast <- Contrast

# write out dsf to excel for prolfquapp
write.xlsx(dsf, "Dataset_UBIenriched.xlsx")
dsf
#
# GO For GlyGly enriched data
#
# get annotation object -> this already contains all potential contrasts for specified control
annotation <- read_annotation(dsf, prefix = "G_")

# working on prot annot
path = "../enriched/"
dir(path)

# get psm and fasta file
files <- get_FP_PSM_files(path = path)
files$fasta
files$data

annot <- annotation$annot
atable <- annotation$atable
annot <- dplyr::mutate(annot, raw.file = gsub("^x|.d.zip$|.raw$",
                                              "", (basename(annot[[atable$fileName]]))))

#filter BGS report for
# 1) keep only GlyGly
# 2) BestLocalizationProbability > 0.75

#check max qvalue
fullBGSreport$EG.Qvalue |> max()
fullBGSreport |> filter(EG.Qvalue < global_qvalue_threshold) |> nrow()
fullBGSreport <- fullBGSreport |> filter(EG.Qvalue < global_qvalue_threshold)


# 1) keep only GlyGly
psmStart <- fullBGSreport |> filter(grepl("GlyGly",x = EG.ModifiedPeptide))

# 2) Localization threshold
psm_long <- psmStart |> filter(`EG.PTMProbabilities [GlyGly (K)]` > best_loc_prob_threshold | `EG.PTMProbabilities [LeuArgGlyGly]` > best_loc_prob_threshold)

# look at enrichment
(enrichment_ratio <- psm_long |> nrow() / fullBGSreport |> nrow())

# 3) prepare for peptidoform summarization
#colnames(psm_long)
psm_relevant <- psm_long|> dplyr::select(R.FileName, PG.ProteinAccessions, PEP.NrOfMissedCleavages, EG.ProteinPTMLocations, EG.ModifiedPeptide,EG.PrecursorId, `EG.TotalQuantity (Settings)`)

if (!is.null(abundance_threshold)) {
  psm_relevant <- dplyr::filter(psm_relevant, `EG.TotalQuantity (Settings)` > abundance_threshold)
}

# generate peptidoFormModSeq for roll-up from precursor to peptide sequence -> replace LeuArgGlyGly with GlyGly
psm_relevant$peptidoFormModSeq <- gsub(x = psm_relevant$EG.ModifiedPeptide, pattern = "LeuArgGlyGly", replacement = "GlyGly")

nrPeptides_exp <- dplyr::summarize(dplyr::group_by(dplyr::distinct(dplyr::select(psm_relevant,
                                                                                 PG.ProteinAccessions, EG.ProteinPTMLocations,  peptidoFormModSeq)), PG.ProteinAccessions), nrPeptides = dplyr::n())
#psm_relevant$site <- paste(psm_relevant$PG.ProteinAccessions, psm_relevant$peptidoFormModSeq, psm_relevant$EG.ProteinPTMLocations, sep = "_")
(colnames(psm_relevant) <- make.names(colnames(psm_relevant)))

# 4) peptidoform summarization
aggregate <- TRUE
colnames(psm_relevant)

if (aggregate) {
  peptidoForm <- psm_relevant |> select(R.FileName,  PG.ProteinAccessions, EG.ProteinPTMLocations,peptidoFormModSeq, EG.TotalQuantity..Settings.) |>
    dplyr::group_by(R.FileName, PG.ProteinAccessions, peptidoFormModSeq,EG.ProteinPTMLocations) |>
    dplyr::summarize(nr_psm = n(), abundance = sum(EG.TotalQuantity..Settings., na.rm = TRUE))
}

head(peptidoForm)

#psm <- peptidoForm
peptidoForm$qValue <- 0.01
peptidoForm$Probability <- 1
nr <- sum(annot[[annotation$atable$sampleName]] %in% sort(unique(peptidoForm$R.FileName)))
logger::log_info("nr : ", nr, " files annotated out of ",
                 length(unique(peptidoForm$R.FileName)))

stopifnot(nr > 0)
colnames(peptidoForm)

atable$ident_Score = "Probability"
atable$ident_qValue = "qValue"
#atable_phos$hierarchy[["protein_Id"]] <- c("ProteinID")
#atable_phos$hierarchy[["site"]] <- c("Index", "Peptide")
atable$hierarchy[["protein_Id"]] <- c("PG.ProteinAccessions")
atable$hierarchy[["site"]] <- c("PG.ProteinAccessions","EG.ProteinPTMLocations")
atable$hierarchy[["mod_peptide_Id"]] <- c("peptidoFormModSeq")
atable$set_response("abundance")
atable$hierarchyDepth <- 2
atable$get_response()

#drop contrast column from sf
annot$ContrastName <- NULL
annot$Contrast <- NULL

# join
psma <- dplyr::inner_join(annot, peptidoForm, join_by("raw.file" == "R.FileName"))

config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(psma, config)
lfqdata <- prolfqua::LFQData$new(adata, config)

fasta_file <- files$fasta

#fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
fasta_annot <- get_annot_from_fasta(fasta_file, pattern_decoys = pattern_decoys)
head(fasta_annot)

nrPeptides_exp$FirstProteinAccession <- sapply(nrPeptides_exp$PG.ProteinAccessions, function(x){ unlist(strsplit(x, "[ ;]"))[1]})

fasta_annot <- dplyr::left_join(nrPeptides_exp, fasta_annot,
                                by = c(FirstProteinAccession = "proteinname"))


fasta_annot <- dplyr::rename(fasta_annot, `:=`(!!lfqdata$config$table$hierarchy_keys_depth()[1],
                                               !!sym("PG.ProteinAccessions")))

fasta_annot <- dplyr::rename(fasta_annot, description = fasta.header)
colnames(fasta_annot)
colnames(lfqdata$data)


# for protein annotation it is important to have Depth 1 (protein level)
lfqdata$config$table$hierarchyDepth <- 1
#ProteinAnnotation$undebug("initialize")

prot_annot <- prolfquapp::ProteinAnnotation$new(
  lfqdata , fasta_annot,
  description = "description",
  cleaned_ids = "FirstProteinAccession",
  full_id = "fasta.id",
  exp_nr_children = "nrPeptides",
  pattern_contaminants = pattern_contaminants,
  pattern_decoys = pattern_decoys
)



lfqdata$remove_small_intensities()

lfqdata$data #checks if all is here
lfqdata$response() #checks if all is here

lfqdata$config$table$hierarchyDepth <- 2

#logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata$config$table$hierarchy_keys()


# gen Grp2
myPath <- getwd()
GRP2_ubi <- prolfquapp::make_DEA_config_R6(PATH = myPath,WORKUNITID = WUID, PROJECTID = fgczProject,
                                           ORDERID = OIDfgcz)


# important for ubi
GRP2_ubi$pop$nr_peptdes <- 1
GRP2_ubi$pop$aggregate <- "sum_topN"
GRP2_ubi$pop$contrasts <- annotation$contrasts


logger::log_info("GENERATING DEA REPORTS")
logger::log_info("starting modelling")

# fix for old grps
GRP2_ubi$pop$aggregate
GRP2_ubi$processing_options$transform <- "robscale"
GRP2_ubi$processing_options$aggregate <- "topN"
GRP2_ubi$pop$transform <- GRP2_ubi$processing_options$transform
GRP2_ubi$pop$Diffthreshold <- GRP2_ubi$processing_options$diff_threshold
GRP2_ubi$pop$FDRthreshold <- GRP2_ubi$processing_options$FDR_threshold



logger::log_info("AGGREGATING PEPTIDE DATA!")
lfqdata <- prolfquapp::aggregate_data(lfqdata, agg_method =
                                        GRP2_ubi$processing_options$aggregate,N = 1000)
logger::log_info("data aggregated: {GRP2_ubi$processing_options$aggregate}.")
logger::log_info("END OF PROTEIN AGGREGATION")

lfqdata$hierarchy_counts()




# grp_ubi
grp_ubi <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2_ubi, prot_annot, GRP2_ubi$pop$contrasts)

# all fine?
grp_ubi$RES$lfqData$to_wide()
grp_ubi$RES$contrastsData_signif

logger::log_info("DONE WITH DEA REPORTS")
# result dir


GRP2_ubi$zipdir_name <- paste0("/DEA_",datetoday, "_", fgczProject, "_", WUID, "_ubi/")
GRP2_ubi$get_zipdir()
dir.create(GRP2_ubi$get_zipdir())

# get markdown files here
# also copy the phospho specific Rmd files from prophosqua
prolfquapp::copy_DEA_Files() # for bib file
prophosqua::copy_phosphoDEA_FragPipe_TMT()

# make sure some info in the reports is correct
grp_ubi$software <- "Spectronaut-GUI"

GRP2 <- grp_ubi
# writing reports
#outpath <- prolfquapp::write_DEA_all(grp_ubi, boxplot = FALSE, markdown = "_Grp2Analysis_Phospho_V2.Rmd")
#outp2 <- prolfquapp::write_DEA_all(grp2 = grp_ubi, boxplot = FALSE, markdown = "_DiffExpQC_Phospho_V2.Rmd")

print(descri)
unique(grp_ubi$RES$contrastsData$contrast)

#outp2 <- prolfquapp::write_DEA_all(grp2 = grp_ubi, boxplot = FALSE, markdown = "_DiffExpQC_Ubi.Rmd")
outpath <- prolfquapp::write_DEA_all(grp_ubi, boxplot = FALSE, markdown = "_Grp2Analysis_Ubi.Rmd") # maybe only one is enough?


# Save RData from enriched and total (only lfqdata is overwritten?) # keep lfqdata, grp, adata separate for phos and total!
(imageFN <- paste(fgczProject, descri, "ubi","Analyzed.RData", sep="_"))
save.image(imageFN)








