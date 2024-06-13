#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# here we want to try prophosqua with simlated data from the MSstatsPTM paper
# https://www.sciencedirect.com/science/article/pii/S1535947622002857#tbl1
# https://github.com/devonjkohler/MSstatsPTM_simulations/tree/main/data

library(tidyverse)
library(prolfqua)
library(prolfquapp)
library(readr)
library(openxlsx)

# func
# write function to get eFDR
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




#load(file = "simulation1_data.rda") # downloaded from github page -> this one is flawed .. 10 features also in the PTM data
load(file="simulation1_data_newByWeW.rda") # this one is fixed")

# source script for different idx 1:8 -> all with 2 grps only but different number of files
# rep(4, 2) rep(6, 2) rep(10, 2) rep(20, 2)

# this takes a while
# for (j in 1:8) {
#   idxOfInterest <- j
#   source("pXXXX_prophosquaForSimData_4sourcing_DEA_and_Integration.R")
# }


# read back results
resultProphosqua <- list()
for (j in 1:8) {
  # folder
  fo <- paste0("Simulation_integration_id_",j)
  fi <- paste0("Integration_id_", j, ".xlsx")
  f <- paste0(fo, "/", fi)
  print(f)
  resultProphosqua[[j]] <- readxl::read_xlsx(f)
}

# get MSstatsPTM results
load("adj_limma_models_sim1.rda")

res_MSstatsPTM <- list()
for (j in 1:8) {
  res_MSstatsPTM[[j]] <- adj_limma_sim1[[j]]
  res_MSstatsPTM[[j]]$adj.P.Val <- p.adjust(res_MSstatsPTM[[j]]$pvalue, method = "BH")
}

# get significant proteins and evaluate eFDR
eFDR_msStats <- list()
eFDR_prophosqua <- list()
sigThreshold <- 0.05

# prolfqua
for (j in 1:8) {
  # get the significant results
  sigProteins <- resultProphosqua[[j]] %>% filter(MSstatsPTMadj_FDR < sigThreshold) %>% pull(protein_Id)
  eFDR_prophosqua[[j]] <- get_empirical_FDR(sigProteins)
}
eFDR_prophosqua


# MSstatsPTM
for (j in 1:8) {
  # get the significant results
  sigProteins <- res_MSstatsPTM[[j]] %>% filter(adj.P.Val < sigThreshold) %>% pull(PTM)
  eFDR_msStats[[j]] <- get_empirical_FDR(sigProteins)
}

# compare performance of eFDRs from msStats and Prophosqua
# create a data frame with the columns eFDR, method, simulation
msSts <- as.data.frame(matrix(unlist(eFDR_msStats), nrow = 8, byrow = TRUE))
msSts$method <- "MSstatsPTM"
msSts$GrpSize <- c(rep(2, 2), rep(3, 2), rep(5, 2), rep(10, 2))

prophosqua_df <- as.data.frame(matrix(unlist(eFDR_prophosqua), nrow = 8, byrow = TRUE))
prophosqua_df$method <- "Prophosqua"
prophosqua_df$GrpSize <- c(rep(2, 2), rep(3, 2), rep(5, 2), rep(10, 2))

# combine the data frames
df <- rbind(msSts, prophosqua_df)
colnames(df) <- c("eFDR", "TP", "FP", "len", "method", "GrpSize")

# plot using ggplot2
ggplot(df, aes(x = GrpSize, y = eFDR, color = method)) +
  geom_point() +
  geom_line() +
  labs(title = "Comparison of eFDRs from MSstatsPTM and Prophosqua",
       x = "Group Size",
       y = "eFDR") +
  theme_minimal()

# show a barplot for each grpsize of the len for each method with respect to the group size with error bars for the len
ggplot(df, aes(x = GrpSize, y = len, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = len - sqrt(len), ymax = len + sqrt(len)), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Comparison of the number of significant proteins from MSstatsPTM and Prophosqua",
       x = "Group Size",
       y = "Number of significant proteins") +
  theme_minimal()

# show a multiple boxplot for each grpsize of the len for each method with respect to the group size
ggplot(df, aes(x = GrpSize, y = len, fill = method)) +
  geom_boxplot() +
  labs(title = "Comparison of the number of significant proteins from MSstatsPTM and Prophosqua",
       x = "Group Size",
       y = "Number of significant proteins") +
  theme_minimal()


library(lattice)
# show a multiple boxplot for each grpsize of the len for each method with respect to the group size
bwplot(len ~ GrpSize | method, data = df, layout = c(1, 2), scales = list(x = list(relation = "free")))


# compare overlap of significant proteins for prolfqua and msStatsPTM

sp <- resultProphosqua[[8]] %>% filter(MSstatsPTMadj_FDR < sigThreshold) %>% pull(protein_Id)
sm <- res_MSstatsPTM[[8]] %>% filter(adj.P.Val < sigThreshold) %>% pull(Protein)
length(intersect(sp, sm))


# get the significant proteins in lists
sigList_prolfqua <- list()
sigList_msStatsPTM <- list()

# prolfqua
for (j in 1:8) {
  # get the significant results
  sigList_prolfqua[[j]] <- resultProphosqua[[j]] %>% filter(MSstatsPTMadj_FDR < sigThreshold) %>% pull(protein_Id)
  sigList_msStatsPTM[[j]] <- res_MSstatsPTM[[j]] %>% filter(adj.P.Val < sigThreshold) %>% pull(Protein)
}

# get intersection of both significant proteins
intersection <- list()
for (j in 1:8) {
  intersection[[j]] <- intersect(sigList_prolfqua[[j]], sigList_msStatsPTM[[j]])
}

# get the length of the intersection
lengths <- sapply(intersection, length)
myIntersection <- rep(lengths, 2)

df <- cbind(df, myIntersection)

# plot using ggplot2
# stacked barplots for each group size with respect to the number of TP, FP and intersection
ggplot(df, aes(x = GrpSize, y = len, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = len - sqrt(len), ymax = len + sqrt(len)), position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Comparison of the number of significant proteins from MSstatsPTM and Prophosqua",
       x = "Group Size",
       y = "Number of significant proteins") +
  theme_minimal()


write_tsv(df, "sim_x_results.tsv")

# preprocess ala witold
# Convert data to long format for plotting
# here instead of  2 cols w TP and FP  we add a column for Type and take counts appart
data_long <- df %>%
  gather(key = "Type", value = "Count", TP, FP)

# ad id per type n method
data_long <- data_long |>
  group_by(method, Type) |>
  mutate(ID = row_number()) |> ungroup()


# visualize w/ stacked barplots
# Plot
pdf("StackedBars_comparison_sim1.pdf",3,3)
ggplot(data_long, aes(x = factor(ID), y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Group Size", y = "Count", fill = "Type") +
  scale_x_discrete(labels = data_long$GrpSize) +
  theme_minimal() +
  facet_wrap(~method, ncol = 1)
dev.off()





