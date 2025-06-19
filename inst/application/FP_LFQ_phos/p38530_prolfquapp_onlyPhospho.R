#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#
# remotes::install_github("prolfqua/prolfquappPTMreaders", dependencies = TRUE)

library(prolfquapp)
library(prolfquappPTMreaders)
library(prophosqua)
library(readr)
library(dplyr)
library(openxlsx)
library(stringr)

# variables
fgczProject <- "pX"
WUID  <- "Phospho"
OIDfgcz <- "o38530"
datetoday <- format(Sys.Date(), "%Y%m%d")
best_loc_prob_threshold <- 0.75
global_qvalue_threshold <- 0.01
abundance_threshold <- 1
PTMtitle <- "OnlyPhospho"


# https://github.com/prolfqua/prolfquappPTMreaders
#ext_reader:
#  extra_args: list(annotation_join_by = c("raw.file", "Name"))
#preprocess: prolfquappPTMreaders::preprocess_FP_combined_STY
#get_files: prolfquappPTMreaders::get_FP_combined_STY_files


# PTM DEA


WUID <- "enriched"

# generate dataset from experimental design
dsfPhosN <- "../phospho_enriched/experiment_annotation.tsv"

# make minimal dsfPhos
dsfPhos <- readr::read_tsv(dsfPhosN)
dsfPhos
dsfPhos$condition <- NULL
dsfPhos$replicate <- NULL
(dsfPhos$file <- gsub(x = dsfPhos$file,pattern = ".*(2025.*_\\d+_S\\d+_.*)\\.d", replacement = "\\1"))

# we further need CONTROL column and Grouping Var
# Grouping Var
((dsfPhos$`Grouping Var` <- gsub(dsfPhos$sample_name, pattern = "\\d+_\\d+_S\\d+_[A-Z]_", replacement = "")))
((dsfPhos$`Grouping Var` <- gsub(dsfPhos$`Grouping Var`, ,pattern = "_.*", replacement = "")))
table(dsfPhos$`Grouping Var`)

# Control
dsfPhos$CONTROL <- "T"
dsfPhos$CONTROL[dsfPhos$`Grouping Var` == "VRK1"] <- "C"
table(dsfPhos$CONTROL, dsfPhos$`Grouping Var`)

# write out annotation table
(fNxls <- paste0(OIDfgcz,"_",WUID,"_annotationTable.xlsx"))
openxlsx::write.xlsx(x = dsfPhos, file = fNxls)

getwd()
# run this in Terminal
# R --vanilla -e "prolfquapp::copy_shell_script(workdir = '.')"
system("chmod a+x prolfqua_*")


# generate yaml
# ./prolfqua_yaml.sh -n vsn -O 38530 -y PTMenriched.yaml -w 328298 -s FP -p 38530



# Yaml for PTM DEA
# extend external readers by:
# ext_reader:
#     extra_args: list()
# preprocess: []
# get_files: []

#   extra_args: list(annotation_join_by = "file")
# preprocess: prolfquappPTMreaders::preprocess_FP_combined_STY
# get_files: prolfquappPTMreaders::get_FP_combined_STY_files


# path: '.'
# zipdir_name: DEA_20250617_PI38530_O38530_WU328298_vsn
# prefix: DEA
# software: DIANN
# project_spec:
#     input_URL: https://fgcz-bfabric.uzh.ch/bfabric/
#     workunit_Id: '328298'
# order_Id: '38530'
# project_name: ''
# project_Id: '38530'
# processing_options:
#     model: prolfqua
# model_missing: yes
# interaction: no
# nr_peptides: 1.0
# pattern_contaminants: ^zz|^CON|Cont_
# pattern_decoys: ^REV_|^rev_
# remove_decoys: no
# remove_cont: no
# FDR_threshold: 0.1
# diff_threshold: 1.0
# aggregate: medpolish
# transform: vsn
# ext_reader:
#     extra_args: list(annotation_join_by = "file")
# preprocess: prolfquappPTMreaders::preprocess_FP_combined_STY
# get_files: prolfquappPTMreaders::get_FP_combined_STY_files
# group: G_
# RES: []
# pop: []


# run this in Terminal
# ./prolfqua_dea.sh -i ../phospho_enriched/ -d o38530_enriched_annotationTable.xlsx -y PTMenriched.yaml -w DEA_PTM -s FP




# 2025-06-18: comment from Antje:
#Hi Jonas,
#wir haben den Verdacht dass zwei Proben vertauscht wurden. Könntest du die Analyse noch mal laufen lassen und manuell die Proben S944602 und S944605 der jeweils anderen Gruppe zuorden? Also:
# 20250606_005_S944602_C_VRK1_treated_replicate_3_enriched zur "treated" Gruppe
# 20250606_006_S944605_F_deadVRK1_treated_replicate_3_enriched zur "control" Gruppe
# Zum einen würde der PCA dann besser passen aber wir sehen das auch in den TIC traces vom chromatogram.
#Danke und viele Grüsse
#
#Antje

# label switches were manually adapted in the excel:

# rerun DEA
# run this in Terminal
system("cd ~/FGCZ/Interdisziplinar/o38530_onlyPhos_LFQ/prophosqua")
# ./prolfqua_dea.sh -i ../phospho_enriched/ -d o38530_enriched_annotationTable_labelSwitch.xlsx -y PTMenriched.yaml -w DEA_PTM_labelSwitched -s FP

# take the switched samples out
# ./prolfqua_dea.sh -i ../phospho_enriched/ -d o38530_enriched_annotationTable_twoSamplesOut.xlsx -y PTMenriched.yaml -w DEA_PTM_twoSamplesOut -s FP


# after meeting with Laura and Antje
# request: rerun, labels shorter, contrast switched
# Yet once more
#./prolfqua_dea.sh -i ../phospho_enriched/ -d o38530_enriched_annotationTable_labelSwitch_v2.xlsx -y PTMenriched.yaml -w DEA_PTM_labelSwitchedv2 -s FP
#./prolfqua_dea.sh -i ../phospho_enriched/ -d o38530_enriched_annotationTable_labelSwitch_v2.xlsx -y PTMenriched_v2_robscale.yaml -w DEA_PTM_labelSwitchedv2 -s FP


# ./prolfqua_dea.sh -i ../phospho_enriched/ -d o38530_enriched_annotationTable_twoSamplesOut_v2.xlsx -y PTMenriched.yaml -w DEA_PTM_twoSamplesOutv2 -s FP
