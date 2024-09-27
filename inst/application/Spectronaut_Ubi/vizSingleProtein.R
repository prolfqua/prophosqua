#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#
library(ggplot2)
library(ggrepel)

myDat <- GRP2$RES$contrastsData

# what is the context
print(descri)

# preprocessing for all proteins
myDat$PTMpos <- sapply(myDat$site, function(x){ unlist(strsplit(x, "~"))[2]})

# clean it
myDat$PTMpos <- gsub(gsub(myDat$PTMpos, pattern = "\\(", replacement = ""), pattern = "\\)", replacement = "")
# get rid of Methionine and number plus C and number (-> keep only Ks)
myDat$PTMpos <- gsub(gsub(gsub(myDat$PTMpos, pattern = "M\\d+", replacement = ""), pattern = "C\\d+", replacement = ""), pattern = ",", replacement = "")
myDat$PTMpos2 <- as.numeric(gsub(myDat$PTMpos, pattern = "K", replacement = ""))

# identify protein of interest
poi <- "PCNA"

# focus on poi only
poiDat <- myDat[myDat$gene_name == poi,]


# loop over all contrasts
for (i in 1:length(unique(poiDat$contrast))){
  COI <- unique(poiDat$contrast)[i]

  cpoi <- poiDat[poiDat$contrast == COI,]

  # draw it
  # even better
  NCyPos <- -0.10
  pp <- ggplot(cpoi, aes(x = PTMpos2, y = diff, color = FDR)) +
    geom_point() +
    geom_text_repel(aes(label = PTMpos), size = 3) +

    # Dashed horizontal line from 0 to protein_length at y = 0
    geom_segment(aes(x = 0, y = 0, xend = protein_length, yend = 0),
                 linetype = "dashed", color = "black") +

    # Vertical lines from each point to the y = 0 line
    geom_segment(aes(x = PTMpos2, xend = PTMpos2, y = diff, yend = 0),
                 color = "gray", linetype = "solid") +

    # Annotate "N" at x = 0 and "C" at protein_length
    annotate("text", x = 0, y = NCyPos, label = "N", size = 5, color = "black") +
    annotate("text", x = cpoi$protein_length, y = NCyPos, label = "C", size = 5, color = "black") +

    theme_minimal() +
    labs(title = paste0(poi, " for contrast <-> ", COI), x = "Position", y = "log2FC")

  fN <- paste(poi,"_",COI, "_NtoC_Ubi.pdf", sep = "")
  pdf(fN, width = 10, height = 10)
  print(pp)
  dev.off()
}

# cpoi <- poiDat |> filter(contrast == COI)
# cpoi$PTMpos <- sapply(cpoi$site, function(x){ unlist(strsplit(x, "~"))[2]})
#
# # clean it
# cpoi$PTMpos <- gsub(gsub(cpoi$PTMpos, pattern = "\\(", replacement = ""), pattern = "\\)", replacement = "")
#
# # get rid of Methionine and number plus C and number (-> keep only Ks)
# cpoi$PTMpos <- gsub(gsub(gsub(cpoi$PTMpos, pattern = "M\\d+", replacement = ""), pattern = "C\\d+", replacement = ""), pattern = ",", replacement = "")
# cpoi$PTMpos2 <- as.numeric(gsub(cpoi$PTMpos, pattern = "K", replacement = ""))

# draw it
library(ggplot2)
library(ggrepel)

# # add a line from 0 to protein_length
# ggplot(cpoi, aes(x = PTMpos2, y = diff, color = FDR)) +
#   geom_point() +
#   geom_text_repel(aes(label = PTMpos), size = 3) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(title = paste0(poi, " ", COI), x = "Position", y = "log2FC") +
#   scale_color_gradient(low = "blue", high = "red")
#
#
# ggplot(cpoi, aes(x = PTMpos2, y = diff, color = FDR)) +
#   geom_point() +
#   geom_text_repel(aes(label = PTMpos), size = 3) +
#   geom_segment(aes(x = 0, y = 0, xend = protein_length, yend = 0),
#                linetype = "dashed", color = "black") +  # Dashed line from 0 to protein_length
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(title = paste0(poi, " ", COI), x = "Position", y = "log2FC") +
#   scale_color_gradient(low = "blue", high = "red")
#
#
# # with vertical lines
# ggplot(cpoi, aes(x = PTMpos2, y = diff, color = FDR)) +
#   geom_point() +
#   geom_text_repel(aes(label = PTMpos), size = 3) +
#
#   # Add dashed horizontal line from 0 to protein_length at y = 0
#   geom_segment(aes(x = 0, y = 0, xend = protein_length, yend = 0),
#                linetype = "dashed", color = "black") +
#
#   # Add vertical lines from each point to the y=0 line
#   geom_segment(aes(x = PTMpos2, xend = PTMpos2, y = diff, yend = 0),
#                color = "gray", linetype = "solid") +
#
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(title = paste0(poi, " for contrast <-> ", COI), x = "Position", y = "log2FC") +
#   scale_color_gradient(low = "blue", high = "red")


# # even better
# NCyPos <- -0.10
# ggplot(cpoi, aes(x = PTMpos2, y = diff, color = FDR)) +
#   geom_point() +
#   geom_text_repel(aes(label = PTMpos), size = 3) +
#
#   # Dashed horizontal line from 0 to protein_length at y = 0
#   geom_segment(aes(x = 0, y = 0, xend = protein_length, yend = 0),
#                linetype = "dashed", color = "black") +
#
#   # Vertical lines from each point to the y = 0 line
#   geom_segment(aes(x = PTMpos2, xend = PTMpos2, y = diff, yend = 0),
#                color = "gray", linetype = "solid") +
#
#   # Annotate "N" at x = 0 and "C" at protein_length
#   annotate("text", x = 0, y = NCyPos, label = "N", size = 5, color = "black") +
#   annotate("text", x = cpoi$protein_length, y = NCyPos, label = "C", size = 5, color = "black") +
#
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(title = paste0(poi, " ", COI), x = "Position", y = "log2FC") +
#   scale_color_gradient(low = "blue", high = "red")
#
#







