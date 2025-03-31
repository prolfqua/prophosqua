#' N to C plot using ggplot2
#' @param POI_matrixMin data.frame
#' @param protein_name name of protein
#' @param protLength protein length
#' @param contrast name of contrast
#' @param thrA significance threshold small default 0.05
#' @param thrB significance threshold small large default 0.20
#'
#' @export
#' @examples
#' data(exampleN_C_dat)
#' N_to_C_plot(exampleN_C_dat,"A0A1I9LPZ1",2160,"H1FC")
#' N_to_C_plot(exampleN_C_dat_no_prot, "A0A178US29",806,"H1FC")
N_to_C_plot <- function(
    POI_matrixMin,
    protein_name,
    protLength,
    contrast,
    thrA = 0.05,
    thrB = 0.2,
    color_protein = "yellow"
) {
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

  POI_matrixMin <- POI_matrixMin |> dplyr::mutate(posInProtein = ifelse(AllLocalized, posInProtein, as.integer(startModSite + endModSite)/2))
  POI_matrixMin <- POI_matrixMin |> dplyr::mutate(modAA = ifelse(AllLocalized, modAA, "NotLoc"))

  mean_diff_prot <- mean(POI_matrixMin$diff.protein, na.rm = TRUE)
  POI_matrixMin$significance <- sapply(POI_matrixMin$FDR.site, get_significance, thrA ,thrB)

  plot_title <- paste0("Prot : ",protein_name, "; length: ", protLength ,"; # sites:", nrow(POI_matrixMin), "; # not localized sites:", sum(!POI_matrixMin$AllLocalized))


  p <- ggplot(data = POI_matrixMin) +
    geom_segment(aes(x = posInProtein, xend = posInProtein, y = 0, yend = diff.site, color = modAA,
                     linetype = linetype)) +
    scale_linetype_manual(values = c("imputed" = "dashed", "observed" = "solid")) +
    annotate("segment", x = 0, xend = protLength, y = 0, yend = 0, color = "black") +
    scale_color_manual(values = c("S" = "blue", "T" = "green", Y = "brown", NotLoc = "pink")) +
    scale_x_continuous(limits = c(0, protLength)) +
    geom_text(aes(x = posInProtein, y = diff.site, label = significance), vjust = 0.4 , size = 7, color = "red") +
    labs(y = paste0("diff : ",contrast),title = plot_title) +
    theme_minimal()



  if (!is.na(mean_diff_prot)) {
    p <- p + annotate("text", x = 0 , y = mean_diff_prot, label = "N", vjust = 0, hjust = 0) +
      annotate("text", x = protLength, y = mean_diff_prot, label = "C", vjust = 0, hjust = 0)

    legend_data <- data.frame(
      xmin = 0,
      xmax = protLength,
      ymin = 0,
      ymax = ifelse(is.na(mean_diff_prot), 0, mean_diff_prot),
      fill = ifelse(is.na(mean_diff_prot), NA, "diff of protein"))

      p <- p + geom_rect(data = legend_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.3) +
        scale_fill_manual(values = c("diff of protein" = color_protein)) +
        guides(fill = guide_legend(title = "Rectangle"))
  } else {
    yext <- max(POI_matrixMin$diff.site, na.rm = TRUE)
    p <- p + annotate("rect", xmin = 0, xmax = protLength, ymin = -yext/2, ymax = +yext/2, alpha = 0.3, fill = "white", color = "red", linetype = "dashed") +
      annotate("text", x = protLength / 2, y = 0, label = "No estimate for diff of protein", color = "red", angle = 45, size = 6)
  }
  return(p)
}


#' N to C for integrated results
#' @param POI_matrixMin
#' @export
#' @examples
#' data(n_c_integrated_df)
#' N_to_C_plot_integrated(n_c_integrated_df,"A0A1I9LT44",539,"WTFC")
N_to_C_plot_integrated <- function(
    POI_matrixMin,
    protein_name,
    protLength,
    contrast,
    thrA = 0.05,
    thrB = 0.2,
    color_protein = "yellow"
) {
  get_significance <- function(fdr, thrA = 0.05, thrB = 0.2) {
    if (fdr < thrA) {
      return("**")
    } else if (fdr < thrB) {
      return("*")
    } else {
      return("")
    }
  }
  POI_matrixMin$significance <- sapply(POI_matrixMin$FDR_I, get_significance, thrA ,thrB)

  plot_title <- paste0("Prot : ",protein_name, "; length: ", protLength ,"; # sites:", nrow(POI_matrixMin))
  mean_diff_prot <- 0
  p <- ggplot(data = POI_matrixMin) +
    geom_segment(aes(x = posInProtein, xend = posInProtein, y = 0, yend = diff_diff, color = modAA)) +
    annotate("segment", x = 0, xend = protLength, y = 0, yend = 0, color = "black") +
    scale_color_manual(values = c("S" = "blue", "T" = "green", Y = "brown", NotLoc = "pink")) +
    scale_x_continuous(limits = c(0, protLength)) +
    annotate("text", x = 0 , y = mean_diff_prot, label = "N", vjust = 0, hjust = 0) +
    annotate("text", x = protLength, y = mean_diff_prot, label = "C", vjust = 0, hjust = 0) +
    geom_text(aes(x = posInProtein, y = diff_diff, label = significance), vjust = 0.4 , size = 7, color = "red") +
    labs(y = paste0("diff : ",contrast),title = plot_title) +
    theme_minimal()

  legend_data <- data.frame(
    xmin = 0,
    xmax = protLength,
    ymin = 0,
    ymax = ifelse(is.na(mean_diff_prot), 0, mean_diff_prot),
    fill = ifelse(is.na(mean_diff_prot), NA, "diff of protein")
  )
  return(p)
}


