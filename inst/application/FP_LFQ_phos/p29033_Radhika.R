library(tidyverse)

files <- prolfquapp::get_FP_combined_STY_files("TKOiWATLFQphospho/")
GRP2 <- prolfquapp::get_config(WORKUNITID = "TKOiWATLFQphospho", ORDERID = "29033")


cfg <- prolfquapp::R6_extract_values(GRP2)
ymlfile <- "Configuration.yml"
yaml::write_yaml(cfg, file = "Configuration.yml")
# if you want to change the parameters edit the Configuration.yml and then run
# GRP2 <- prolfquapp::get_config("Configuration.yml" ,WORKUNITID = "TKOiWATLFQphospho")


opt <- list()
opt$dataset <- "dataset.xlsx"

# DO NOT EDIT BELOW

annotation <- file.path(opt$dataset) |>
  prolfquapp::read_annotation_file() |> prolfquapp::read_annotation(prefix = GRP2$group)
logger::log_info("Contrasts: \n", paste(annotation$contrasts, collapse = "\n"))


prolfquapp::copy_DEA_Files()


xd <- prolfquapp::preprocess_FP_combined_STY(
  files$data,
  files$fasta,
  annotation,
  pattern_contaminants = "^zz|^CON",
  pattern_decoys = "REV_",
  annotation_join_by = "Name"
)

logger::log_info(paste(c("Protein Annotation :\n",capture.output( print(xd$protein_annotation$get_summary()))),collapse = "\n"))
logger::log_info("AGGREGATING PEPTIDE DATA: {GRP2$processing_options$aggregate}.")
lfqdata <- prolfquapp::aggregate_data(xd$lfqdata, agg_method = GRP2$processing_options$aggregate)
logger::log_info("END OF PROTEIN AGGREGATION")
logger::log_info("RUN ANALYSIS")
grp <- prolfquapp::generate_DEA_reports2(lfqdata, GRP2, xd$protein_annotation, annotation$contrasts)
logger::log_info("Writing results to: " ,  GRP2$get_zipdir())

outdir <- prolfquapp::write_DEA_all(
  grp, boxplot = FALSE, markdown = "_Grp2Analysis_V2.Rmd")

lfqdataIB <- xd$lfqdata$get_subset(xd$protein_annotation$clean(
  contaminants = GRP2$processing_options$remove_cont,
  decoys = GRP2$processing_options$remove_decoys))
# do not write when peptide level analysis.
if (length(xd$protein_annotation$pID) == 1) {
  ibaq <- prolfquapp::compute_IBAQ_values(lfqdataIB, xd$protein_annotation)
  writexl::write_xlsx(ibaq$to_wide()$data, path = file.path(grp$get_result_dir(), "IBAQ.xlsx"))
}

logger::log_info("Writing summarized experiment.")
SE <- prolfquapp::make_SummarizedExperiment(grp)
saveRDS(SE, file = file.path( grp$get_result_dir(), "SummarizedExperiment.rds"))


logger::log_info("Creating directory with input files :", GRP2$get_input_dir())
dir.create(GRP2$get_input_dir())
prolfquapp::copy_DEA_Files(workdir = GRP2$get_input_dir())
prolfquapp::copy_shell_script(workdir = GRP2$get_input_dir())

file.copy(c(files$data, files$fasta, ymlfile, opt$dataset), GRP2$get_input_dir())
logger::log_info("Wirte yaml with parameters: ", file.path(GRP2$get_input_dir(), "minimal.yaml"))


GRP2$RES <- NULL
GRP2$pop <- NULL
yaml::write_yaml(prolfquapp::R6_extract_values(GRP2), file = file.path(GRP2$get_input_dir(), "minimal.yaml"))


