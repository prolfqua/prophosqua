message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")

# annotation and comparison based on:
# https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details

# some functions
# write function to get spezificity and sensitivity
get_spezificity_sensitivity <- function(objectWithProteinIDs){
  # get the true positives
  TP <- sum(str_count(string = objectWithProteinIDs, pattern = "NoChange")==0)
  # get the false positives
  FP <- sum(str_count(string = objectWithProteinIDs, pattern = "NoChange")>0)
  # get the false negatives
  FN <- TP - 250
  # get the true negatives
  TN <- 750 - FP
  # length
  len <- length(objectWithProteinIDs)
  # get the spezificity
  spezificity <- TN / (TN + FP)
  # get the sensitivity
  sensitivity <- TP / (TP + FN)
  # get precision
  precision <- TP / (TP + FP)
  # get recall
  recall <- TP / (TP + FN)
  return(c(spezificity, sensitivity, precision, recall, len))
}

# write function to get spezificity and sensitivity
get_empirical_FDR <- function(objectWithProteinIDs){
  # get the true positives
  TP <- sum(str_count(string = objectWithProteinIDs, pattern = "NoChange")==0)
  # get the false positives
  FP <- sum(str_count(string = objectWithProteinIDs, pattern = "NoChange")>0)
  # length
  len <- length(objectWithProteinIDs)
  # get eFDR
  eFDR <- FP / (TP + FP)
  return(c(eFDR, TP, FP, len))
}


# look into prophosqua results to see differences in TP and FP
res_prophosqua <- read.xlsx(xlsxFile = "pXXXX_TMTphospho_integration_ownWeWSimulatedOne/Integration_ownWeWSimulatedOne.xlsx", sheet = "combinedStats")

# read in msstatsPTM results
load("adj_limma_models_sim1.rda")
res_MSstatsPTM <- adj_limma_sim1[[1]]
head(res_MSstatsPTM)
res_MSstatsPTM$adj.P.Val <- p.adjust(res_MSstatsPTM$pvalue, method = "BH")


sigThreshold <- 0.1
# get the significant results
sig_prophosqua <- res_prophosqua[res_prophosqua$MSstatsPTMadj_FDR < sigThreshold,]
sig_MSstatsPTM <- res_MSstatsPTM[res_MSstatsPTM$adj.P.Val < sigThreshold,]

get_spezificity_sensitivity(objectWithProteinIDs = sig_MSstatsPTM$PTM)
get_spezificity_sensitivity(objectWithProteinIDs = sig_prophosqua$IDcolumn.x)

get_empirical_FDR(objectWithProteinIDs = sig_MSstatsPTM$PTM)
get_empirical_FDR(objectWithProteinIDs = sig_prophosqua$IDcolumn.x)

# get the overlap
overlap <- intersect(sig_prophosqua$IDcolumn.x, sig_MSstatsPTM$PTM)
length(overlap)

