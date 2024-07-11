#' N to C plot using ggplot2
#' @export
#'
N_to_C_plot <- function(POI_matrixMin, protein_name, protLength, contrast,thrA = 0.05, thrB = 0.2) {
  get_significance <- function(fdr, thrA = 0.05, thrB = 0.2) {
    if (fdr < thrA) {
      return("**")
    } else if (fdr < thrB) {
      return("*")
    } else {
      return("")
    }
  }


  POI_matrixMin$linetype <- ifelse(grepl("imputed",POI_matrixMin$model_site), "imputed", "observed")
  class(POI_matrixMin[["startModSite"]])  <- "numeric"
  class(POI_matrixMin[["endModSite"]])  <- "numeric"

  POI_matrixMin <- POI_matrixMin |> mutate(posInProtein = ifelse(AllLocalized, posInProtein, as.integer(startModSite + endModSite)/2))
  POI_matrixMin <- POI_matrixMin |> mutate(modAA = ifelse(AllLocalized, modAA, "NotLoc"))

  mean_diff_prot <- mean(POI_matrixMin$diff.prot, na.rm = TRUE)
  POI_matrixMin$significance <- sapply(POI_matrixMin$FDR.site, get_significance, thrA ,thrB)

  plot_title <- paste0("Prot : ",protein_name, "; length: ", protLength ,"; # sites:", nrow(POI_matrixMin), "; # not localized sites:", sum(!POI_matrixMin$AllLocalized))


  p <- ggplot(data = POI_matrixMin) +
    geom_segment(aes(x = posInProtein, xend = posInProtein, y = 0, yend = diff.site, color = modAA,
                     linetype = linetype)) +
    scale_linetype_manual(values = c("imputed" = "dashed", "observed" = "solid")) +
    annotate("segment", x = 0, xend = protLength, y = 0, yend = 0, color = "black") +
    scale_color_manual(values = c("S" = "blue", "T" = "green", Y = "brown", NotLoc = "pink")) +
    scale_x_continuous(limits = c(0, protLength)) +
    # protein box
    # annotate("rect", xmin = 0, xmax = protLength, ymin = 0, ymax = mean_diff_prot, alpha = 0.3, fill = "yellow") +
    annotate("text", x = 0 , y = mean_diff_prot, label = "N", vjust = 0, hjust = 0) +
    annotate("text", x = protLength, y = mean_diff_prot, label = "C", vjust = 0, hjust = 0) +
    geom_text(aes(x = posInProtein, y = diff.site, label = significance), vjust = 0 , size = 7, color = "red") +
    labs(y = paste0("diff : ",contrast),title = plot_title) +
    theme_minimal()

  legend_data <- data.frame(
    xmin = 0,
    xmax = protLength,
    ymin = 0,
    ymax = ifelse(is.na(mean_diff_prot), 0, mean_diff_prot),
    fill = ifelse(is.na(mean_diff_prot), NA, "diff of protein")
  )

  if (!is.na(mean_diff_prot)) {
    #p <- p + annotate("rect", xmin = 0, xmax = protLength, ymin = 0, ymax = mean_diff_prot, alpha = 0.3, fill = "yellow")
    p <- p + geom_rect(data = legend_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.3) +
      scale_fill_manual(values = c("diff of protein" = "yellow")) +
      guides(fill = guide_legend(title = "Rectangle"))
  } else {
    yext <- max(POI_matrixMin$diff.site, na.rm = TRUE)
    p <- p + annotate("rect", xmin = 0, xmax = protLength, ymin = -yext/2, ymax = +yext/2, alpha = 0.3, fill = "white", color = "red", linetype = "dashed") +
      annotate("text", x = protLength / 2, y = 0, label = "No estimate for diff of protein", color = "red", angle = 45, size = 6)
  }
  return(p)
}

