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


#load(file = "simulation1_data.rda") # downloaded from github page -> this one is flawed .. 10 features also in the PTM data
load(file="simulation1_data_newByWeW.rda") # this one is fixed")


for (j in 1:8) {
  idxOfInterest <- j
  source("pXXXX_prophosquaForSimData_sourcing_DEA_and_Integration_simByWeW.R")
}


# read back results
resultProlfqua <- list()

idx <- 1:8

idx[2]

# folder
ff <- paste0("Simulation_integration_id_", idx[2])
fi <- paste0("Integration_id_", idx[2], "xlsx")

f <- paste0(ff, "/", fi)

resultProlfqua[[idx[2]]] readxl::read_xlsx(f)
