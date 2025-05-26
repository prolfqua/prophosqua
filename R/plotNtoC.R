#' @importFrom rlang .data
NULL

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
#' n_to_c_plot(exampleN_C_dat, "A0A1I9LPZ1", 2160, "H1FC")
#' n_to_c_plot(exampleN_C_dat_no_prot, "A0A178US29", 806, "H1FC")
n_to_c_plot <- function(
    poi_matrix_min,
    protein_name,
    prot_length,
    contrast,
    thr_a = 0.05,
    thr_b = 0.2,
    color_protein = "yellow") {
  get_significance <- function(fdr, thr_a = 0.05, thr_b = 0.2) {
    if (fdr < thr_a) {
      return("**")
    } else if (fdr < thr_b) {
      return("*")
    } else {
      return("")
    }
  }
  poi_matrix_min$linetype <- ifelse(grepl("imputed", poi_matrix_min$model_site), "imputed", "observed")
  class(poi_matrix_min[["startModSite"]]) <- "numeric"
  class(poi_matrix_min[["endModSite"]]) <- "numeric"

  poi_matrix_min <- poi_matrix_min |>
    dplyr::mutate(
      posInProtein =
        ifelse(.data$AllLocalized, .data$posInProtein, as.integer(.data$startModSite + .data$endModSite) / 2)
    )
  poi_matrix_min <- poi_matrix_min |> dplyr::mutate(modAA = ifelse(.data$AllLocalized, .data$modAA, "NotLoc"))

  mean_diff_prot <- mean(poi_matrix_min$diff.protein, na.rm = TRUE)
  poi_matrix_min$significance <- sapply(poi_matrix_min$FDR.site, get_significance, thr_a, thr_b)

  plot_title <- paste0(
    "Prot : ", protein_name,
    "; length: ", prot_length, "; # sites:",
    nrow(poi_matrix_min), "; # not localized sites:", sum(!poi_matrix_min$AllLocalized)
  )


  p <- ggplot(data = poi_matrix_min) +
    geom_segment(aes(
      x = .data$posInProtein, xend = .data$posInProtein, y = 0, yend = .data$diff.site, color = .data$modAA,
      linetype = .data$linetype
    )) +
    scale_linetype_manual(values = c("imputed" = "dashed", "observed" = "solid")) +
    annotate("segment", x = 0, xend = prot_length, y = 0, yend = 0, color = "black") +
    scale_color_manual(values = c("S" = "blue", "T" = "green", Y = "brown", NotLoc = "pink")) +
    scale_x_continuous(limits = c(0, prot_length)) +
    geom_text(aes(x = .data$posInProtein, y = .data$diff.site, label = .data$significance),
      vjust = 0.4, size = 7, color = "red"
    ) +
    labs(y = paste0("diff : ", contrast), title = plot_title) +
    theme_minimal()



  if (!is.na(mean_diff_prot)) {
    p <- p + annotate("text", x = 0, y = mean_diff_prot, label = "N", vjust = 0, hjust = 0) +
      annotate("text", x = prot_length, y = mean_diff_prot, label = "C", vjust = 0, hjust = 0)

    legend_data <- data.frame(
      xmin = 0,
      xmax = prot_length,
      ymin = 0,
      ymax = ifelse(is.na(mean_diff_prot), 0, mean_diff_prot),
      fill = ifelse(is.na(mean_diff_prot), NA, "diff of protein")
    )

    p <- p + geom_rect(data = legend_data, aes(
      xmin = .data$xmin,
      xmax = .data$ xmax, ymin = .data$ymin,
      ymax = .data$ymax, fill = .data$fill
    ), alpha = 0.3) +
      scale_fill_manual(values = c("diff of protein" = color_protein)) +
      guides(fill = guide_legend(title = "Rectangle"))
  } else {
    yext <- max(poi_matrix_min$diff.site, na.rm = TRUE)
    p <- p + annotate("rect",
      xmin = 0, xmax = prot_length,
      ymin = -yext / 2, ymax = +yext / 2, alpha = 0.3,
      fill = "white", color = "red", linetype = "dashed"
    ) +
      annotate("text",
        x = prot_length / 2,
        y = 0, label = "No estimate for diff of protein", color = "red", angle = 45, size = 6
      )
  }
  return(p)
}


#' N to C for integrated results
#' @param POI_matrixMin data.frame
#' @export
#' @examples
#' data(n_c_integrated_df)
#' n_to_c_plot_integrated(n_c_integrated_df, "A0A1I9LT44", 539, "WTFC")
n_to_c_plot_integrated <- function(
    poi_matrix_min,
    protein_name,
    prot_length,
    contrast,
    thr_a = 0.05,
    thr_b = 0.2,
    color_protein = "yellow") {
  get_significance <- function(fdr, thr_a = 0.05, thr_b = 0.2) {
    if (fdr < thr_a) {
      return("**")
    } else if (fdr < thr_b) {
      return("*")
    } else {
      return("")
    }
  }
  poi_matrix_min$significance <- sapply(poi_matrix_min$FDR_I, get_significance, thr_a, thr_b)

  plot_title <- paste0("Prot : ", protein_name, "; length: ", prot_length, "; # sites:", nrow(poi_matrix_min))
  mean_diff_prot <- 0
  p <- ggplot(data = poi_matrix_min) +
    geom_segment(aes(
      x = .data$posInProtein,
      xend = .data$posInProtein, y = 0, yend = .data$diff_diff,
      color = .data$modAA
    )) +
    annotate("segment", x = 0, xend = prot_length, y = 0, yend = 0, color = "black") +
    scale_color_manual(values = c("S" = "blue", "T" = "green", Y = "brown", NotLoc = "pink")) +
    scale_x_continuous(limits = c(0, prot_length)) +
    annotate("text", x = 0, y = mean_diff_prot, label = "N", vjust = 0, hjust = 0) +
    annotate("text", x = prot_length, y = mean_diff_prot, label = "C", vjust = 0, hjust = 0) +
    geom_text(aes(x = .data$posInProtein, y = .data$diff_diff, label = .data$significance),
      vjust = 0.4, size = 7, color = "red"
    ) +
    labs(y = paste0("diff : ", contrast), title = plot_title) +
    theme_minimal()

  return(p)
}
