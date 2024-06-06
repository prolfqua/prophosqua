message("prolfqua Version :", packageVersion("prolfqua"), "\n")
message("prolfquaapp Version :", packageVersion("prolfquapp"), "\n")

# annotation and comparison based on:
# https://fgcz-bfabric.uzh.ch/bfabric/dataset/show.html?id=47399&tab=details


# libs
library(prolfqua)
library(prolfquapp)
library(prophosqua)
library(readr)
library(dplyr)
library(stringr)
library(openxlsx)
library(seqinr)


# parameters and thresholds
# params ideally taken from yaml
fgczProject <- "pXXXX"
descri <- "TMTphospho_integration"
compari <- "_SimulatedOne"
WUID <- "WUxx"

# read back in results
#paste0(resDir,"Results_DEA_WU",WUID,"/DE_Groups_vs_Controls.xlsx" )
totRes <- read.xlsx(xlsxFile = "SimulationONETMTphospho_TotalProteome_PI_pXXXX_OI_oYYYY_WU__none/Results_DEA_WU/DE_Groups_vs_Controls_WU.xlsx", sheet = "diff_exp_analysis")
phosRes <- read.xlsx(xlsxFile = "SimulationONETMTphospho_PhosphoEnriched_WU__none/Results_DEA_WU/DE_Groups_vs_Controls_WU.xlsx", sheet = "diff_exp_analysis")

(resDir <- paste0(fgczProject, "_",descri, compari))

# there is no fasta for simlulated data
# myFasta <- seqinr::read.fasta("../fgcz_9606_reviewed_cnl_20230330.fasta", seqtype = "AA", as.string = TRUE)

# work on phosRes
head(phosRes$site)

# Filter out FGCZContaminants -> problematic at sites
phosRes <- phosRes |> filter(!grepl("FGCZCont", protein_Id))


#
phosRes$originalSite <- phosRes$site
#phosRes$peptideSequence <- sapply(strsplit((phosRes$site), split = "~"), function(x)x[2])
phosRes$site <- sapply(strsplit((phosRes$site), split = "~"), function(x)x[1])
phosRes$AccFromSite <- sapply(strsplit((phosRes$site), split = "_"), function(x)x[1])
phosRes$startModSite <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[2]))
phosRes$endModSite <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[3]))
#phosRes$NumPhos <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[4]))
#phosRes$LocalizedNumPhos <- as.numeric(sapply(strsplit((phosRes$site), split = "_"), function(x)x[5]))
phosRes$NumPhos <- 1
phosRes$LocalizedNumPhos <- 1
phosRes$PhosSites <- gsub(x = phosRes$site, pattern = "_",replacement = "")
phosRes$SinglePhos_bool <- phosRes$NumPhos == 1
phosRes$AllLocalized <- phosRes$NumPhos == phosRes$LocalizedNumPhos
phosRes$SinglePhosLocalized_bool <- phosRes$NumPhos == 1 & phosRes$LocalizedNumPhos == 1
phosRes$peptideSequence <- "PEPTIDEK"

# work on phospho to
phosRes$AA <- gsub(x = phosRes$PhosSites, pattern = "\\d+", replacement = "")
# parse position in protein
phosRes$posInProtein <- as.numeric(gsub(x = phosRes$PhosSites, pattern = "[STY]", replacement = ""))

# add info if found in more than one unique peptide
phosRes$SiteFoundInManyPeptides <- 1

table(phosRes$AA[phosRes$SinglePhosLocalized_bool])
head(phosRes$posInProtein[phosRes$SinglePhosLocalized_bool])

# find mapping for sequence window
phosRes |> select(protein_Id, peptideSequence, posInProtein) |> distinct() |> dim_desc()
phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct() |> dim_desc()

uniqueProtPepSeq <- phosRes |> filter(SinglePhosLocalized_bool == TRUE) |> select(protein_Id, peptideSequence, posInProtein) |> distinct()

# we do not have a fasta
# for (i in 1:nrow(uniqueProtPepSeq)) {
#   uniqueProtPepSeq$SequenceWindows[i] <- getSequenceWindowForLogo(poi = uniqueProtPepSeq$protein_Id[i], soi = uniqueProtPepSeq$posInProtein[i], fastaObject = myFasta)
# }

uniqueProtPepSeq$SequenceWindows <- "MYPEPTIDEKWINDWWFAKE"
tail(uniqueProtPepSeq)

# join sequence windows back!
phosRes <- left_join(x = phosRes, y = uniqueProtPepSeq)

# join sheets
# get rid  of decoys
phosRes  <- phosRes |> filter(!grepl("REV_", protein_Id))
totRes  <- totRes |> filter(!grepl("REV_", protein_Id))

# parse  phosRes from "Protein_128|NoChange1_S_1" back to totRes$protein_Id[2] -> "Protein_1|NoChange1"

phosRes$fProt <- sapply(strsplit((phosRes$protein_Id), split = "\\|"), function(x)x[1])
# cut away after _.* -> "Protein_1"
# q: how can I do gsub with back reference in R? -> \\1
phosRes$fProt <- gsub(x = phosRes$fProt, pattern = "(Protein_\\d+).*", replacement = "\\1")
phosRes$Accpart <- sapply(strsplit((phosRes$protein_Id), split = "\\|"), function(x)x[2])
# replace NA with ""
phosRes$Accpart[is.na(phosRes$Accpart)] <- ""
phosRes$Accpart <- gsub(x = phosRes$Accpart, pattern = "(NoChange\\d+).*", replacement = "\\1")

phosRes$cleanPhosProtein <- paste(phosRes$fProt, phosRes$Accpart, sep = "|")
#problematic Protein_1| for no Accpart
phosRes$cleanPhosProtein <- gsub(x = phosRes$cleanPhosProtein[phosRes$Accpart == ""], pattern = "\\|", replacement = "")

# combine
combo <- left_join(x = phosRes, y = totRes, join_by("cleanPhosProtein" == "protein_Id", "contrast" == "contrast"))


# MS stats like adjustment for protein change
comboWithAdj <- doMSstatsLikeSiteNormalizationUsingProteinStatsOnComboObject(combo)
colnames(comboWithAdj)
head(comboWithAdj)

# save.image("Integration_o33038_allIn.RData")
# load("Integration_vs5hpiNoGSK.RData")

# for RMD report
GRP2 <- prolfquapp::make_DEA_config(PROJECTID = fgczProject, ORDERID = fgczProject, WORKUNITID = "WUxxx")

# Idea for reporting and downstream processing

# Number crunching
# number of sites split STY
# all of them
#barplot(table(comboWithAdj$NumPhos[comboWithAdj$AllLocalized]), main = "Number of identified Phosphorylations per Peptide")

# only single phos-localized -> in RMD
# table(comboWithAdj$AA[comboWithAdj$SinglePhosLocalized_bool])
#barplot(table(comboWithAdj$AA[comboWithAdj$SinglePhosLocalized_bool]))


# Abundances and STY -> in RMD
# p <- ggplot2::ggplot(na.omit(comboWithAdj[comboWithAdj$SinglePhosLocalized_bool,]), ggplot2::aes(x=AA, y=avgAbd.x)) +
#   ggplot2::geom_boxplot()
# p


# separate for each constrast
(NumCoi <- length(unique(comboWithAdj$contrast)))


# N-to-C plotting per contrast
# fdrThreshold <- 0.1 # for plotting N-to-C
#
# for (i in 1:length(unique(comboWithAdj$contrast))) {
#   comboMat <- comboWithAdj[comboWithAdj$contrast == unique(comboWithAdj$contrast)[i],]
#   # Towards N-to-C plots
#   sum(comboWithAdj$FDR.x < fdrThreshold)
#   candidateMat <- comboWithAdj[comboWithAdj$FDR.x < fdrThreshold,]
#   generateNtoCProteinPDFsWithPhosphoPeptides_FragPipeTMT(globalNphosphoCombinedNResultMatrix = comboMat, candidateMatrix = candidateMat, expName = unique(comboWithAdj$contrast)[i])
# }


# # render integration html
# (resultPath <- resDir)
# (htmlFN <- paste0("Integration",compari))
# comboWithAdj$protein_Id <- comboWithAdj$protein_Id.x
# prolfquapp::render_DEA(GRP2 = comboWithAdj, outpath = resultPath, htmlname = htmlFN, word = FALSE, markdown = "_Overview_PhosphoAndIntegration.Rmd")

# write to excel
#prolfquapp::write_DEA(GRP2 = comboWithAdj, outpath = "IntegrationResults", xlsxname = "IntegrationPhosphoCentric", write = TRUE)
excelResultList <- list()
excelResultList$combinedStats <- comboWithAdj
writexl::write_xlsx(excelResultList, path = paste0(resultPath, "/",htmlFN,".xlsx"))





# look into ms-stats pmt results --> wtf these results are not feature centric!
load("adj_limma_models_sim1.rda")
# look into adj_limma_models_sim1
adj_limma_sim1[[1]] |> dim()
tt <-  adj_limma_sim1[[1]]
load("ptm_models_sim1.rda")
# look into ptm_models_sim1
ptm_models_sim1[[1]] |> dim()
adj_limma_sim1[[1]]


# look into prophosqua results to see differences in TP and FP
res_prophosqua <- read.xlsx(xlsxFile = "pXXXX_TMTphospho_integration_SimulatedOne/Integration_SimulatedOne.xlsx", sheet = "combinedStats")

sigThreshold <- 0.05

res_prophosqua$potentialTP_bool <- str_count(res_prophosqua$site ,"NoChange") == 0
res_prophosqua$is_sig <- res_prophosqua$MSstatsPTMadj_FDR < sigThreshold
res_prophosqua$sureFP <- str_count(res_prophosqua$site ,"NoChange") > 0 & res_prophosqua$MSstatsPTMadj_FDR < sigThreshold
res_prophosqua$likelyTP <- res_prophosqua$potentialTP_bool &  res_prophosqua$is_sig

# q: draw two histogram for res_prophosqua$avgAbd.x where likelyTP and sureFPlibrary(ggplot2)
library(ggplot2)
# Filter the data based on conditions
likelyTP_hist <- res_prophosqua[res_prophosqua$likelyTP == TRUE, ]
sureFP_hist <- res_prophosqua[res_prophosqua$sureFP == TRUE, ]

# Create histograms
ggplot() +
  # Histogram for likelyTP
  geom_histogram(data = likelyTP_hist, aes(x = avgAbd.x), fill = "blue", alpha = 0.5, binwidth = 0.5) +
  # Histogram for sureFP
  geom_histogram(data = sureFP_hist, aes(x = avgAbd.x), fill = "red", alpha = 0.5, binwidth = 0.5) +
  # Axis and plot labels
  labs(x = "avgAbd.x", y = "Frequency", title = "Histograms for likelyTP sureFP") +
  # Legend
  scale_fill_manual(values = c("blue", "red"), labels = c("likelyTP", "sureFP"))


# q: how can I have more granularity on my histogram (more breaks)?
ggplot() +
  # Histogram for likelyTP
  geom_histogram(data = likelyTP_hist, aes(x = avgAbd.x), fill = "blue", alpha = 0.5, binwidth = 0.05) +
  # Histogram for sureFP
  geom_histogram(data = sureFP_hist, aes(x = avgAbd.x), fill = "red", alpha = 0.5, binwidth = 0.05) +
  # Axis and plot labels
  labs(x = "avgAbd.x", y = "Frequency", title = "Histograms for likelyTP sureFP") +
  # Legend
  scale_fill_manual(values = c("blue", "red"), labels = c("likelyTP", "sureFP")) +
  # More breaks
  scale_x_continuous(breaks = seq(0, 10, by = 0.01))


# now I want the same but look at the diff.x values
ggplot() +
  # Histogram for likelyTP
  geom_histogram(data = likelyTP_hist, aes(x = diff.x), fill = "blue", alpha = 0.5, binwidth = 0.05) +
  # Histogram for sureFP
  geom_histogram(data = sureFP_hist, aes(x = diff.x), fill = "red", alpha = 0.5, binwidth = 0.05) +
  # Axis and plot labels
  labs(x = "diff.x", y = "Frequency", title = "Histograms for likelyTP sureFP") +
  # Legend
  scale_fill_manual(values = c("blue", "red"), labels = c("likelyTP", "sureFP")) +
  # More breaks
  scale_x_continuous(breaks = seq(-10, 10, by = 0.5))

# q: I want to draw xy-plot in ggplot2 where x-axis is diff.x and y-axis is diff.y
# q: how can I use continous color of -log10(MSstatsPTMadj_FDR) for the points in the scatter plot?
# q: how can I also use -log10(MSstatsPTMadj_FDR) size for the dots in the scatter plot?
# q: how can I add alpha blending to the points in the scatter plot or make them transparent or invert the color scale?

# q: how can I change color range for the points in the scatter plot?
# q: how can I invert the color range for the points in the scatter plot?
# q: how can I change the shape of the points in the scatter plot?  https://ggplot2.tidyverse.org/reference/geom_point.html
# q: how can I add a red circle around plots where MSstatsPTMadj_FDR < 0.05?
# a: https://ggplot2.tidyverse.org/reference/geom_point.html

ggplot() +
  # Scatter plot
  geom_point(alpha = 1/0.1, data = res_prophosqua, aes(x = diff.x, y = diff.y,
                                                        color = -log10(MSstatsPTMadj_FDR),
                                                       size = -log10(MSstatsPTMadj_FDR),
                                                       shape = factor(is_sig))) +
  scale_colour_gradient2() +
  # Axis and plot labels
  labs(x = "diff.x", y = "diff.y", title = "Scatter plot for diff.x and diff.y")


ggplot() +
  # Scatter plot
  geom_point(alpha = 1/10, data = res_prophosqua, aes(x = diff.x, y = diff.y, color = -log10(MSstatsPTMadj_FDR), size = 0.01)) +
  scale_colour_gradient2() +
  # Axis and plot labels
  labs(x = "diff.x", y = "diff.y", title = "Scatter plot for diff.x and diff.y")


ggplot() +
  # Scatter plot
  geom_point(alpha = 1/0.1, data = res_prophosqua, aes(x = diff.x, y = diff.y, color = -log10(MSstatsPTMadj_FDR), size = -log10(MSstatsPTMadj_FDR))) +
  scale_colour_gradient(low = "lightblue", high = "black") +
  # Axis and plot labels
  labs(x = "diff.x", y = "diff.y", title = "Scatter plot for diff.x and diff.y")

ggplot() +
  # Scatter plot
  geom_point(data = res_prophosqua, aes(x = diff.x, y = diff.y, color = -log10(MSstatsPTMadj_FDR))) +
  # Axis and plot labels
  labs(x = "diff.x", y = "diff.y", title = "Scatter plot for diff.x and diff.y")

ggplot() +
  # Scatter plot
  geom_point(data = res_prophosqua, aes(x = diff.x, y = diff.y), color = "blue") +
  # Axis and plot labels
  labs(x = "diff.x", y = "diff.y", title = "Scatter plot for diff.x and diff.y")




