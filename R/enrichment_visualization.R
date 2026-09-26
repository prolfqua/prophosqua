#' Enrichment Visualization Functions
#'
#' Plots shared by the PTM-SEA, Kinase GSEA and MEA sections of the enrichment
#' report.
#'
#' @importFrom dplyr mutate filter arrange group_by summarize case_when pull slice_min
#' @importFrom ggplot2 ggplot aes geom_point geom_tile geom_text geom_vline
#'   geom_hline scale_color_gradient2 scale_fill_gradient2
#'   scale_alpha_manual scale_color_manual facet_wrap theme_bw theme_minimal
#'   theme element_text labs
#' @importFrom stats reorder
#' @importFrom rlang .data
#' @importFrom purrr map_dfr
#' @name enrichment_visualization
NULL

# Adds -log10(FDR), the direction of the NES and whether the FDR passes the
# threshold.
prepare_enrichment_data <- function(data, fdr_col = "FDR", fdr_threshold = 0.1) {
  data |>
    dplyr::mutate(
      neg_log_fdr = -log10(pmax(.data[[fdr_col]], 1e-10)),
      direction = dplyr::case_when(
        .data$NES > 0 ~ "Up",
        .data$NES < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      significant = .data[[fdr_col]] < fdr_threshold
    )
}

#' Create volcano plot for enrichment results
#'
#' @param data Data frame with columns: NES, p.adjust/FDR, contrast, item (kinase/pathway)
#' @param item_col Name of item column for labels (default: "kinase")
#' @param fdr_col Name of FDR column (default: "FDR")
#' @param fdr_threshold FDR threshold for significance line (default: 0.1)
#' @param label_fdr_threshold FDR threshold for labeling points (default: 0.05)
#' @param n_labels Number of top labels per contrast (default: 5)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return ggplot object
#' @export
#' @examples
#' kinase_results <- data.frame(
#'   kinase = rep(c("PKACA", "AKT1", "MAPK1", "CDK1"), 2),
#'   NES = c(2.1, 1.5, -1.8, 0.5, 1.2, -0.8, 1.9, -1.1),
#'   FDR = c(0.001, 0.02, 0.01, 0.3, 0.05, 0.2, 0.008, 0.06),
#'   contrast = rep(c("Treatment_vs_Control", "Timepoint_vs_Baseline"), each = 4)
#' )
#' plot_enrichment_volcano(kinase_results, title = "Kinase Enrichment Volcano")
plot_enrichment_volcano <- function(
  data,
  item_col = "kinase",
  fdr_col = "FDR",
  fdr_threshold = 0.1,
  label_fdr_threshold = 0.05,
  n_labels = 5,
  title = NULL,
  subtitle = NULL
) {
  volcano_data <- prepare_enrichment_data(data, fdr_col, fdr_threshold)
  if (is.null(subtitle)) {
    subtitle <- paste0("Dashed line: FDR = ", fdr_threshold, "; Labels: top ", n_labels, " by FDR per contrast")
  }
  label_data <- volcano_data |>
    dplyr::filter(.data[[fdr_col]] < label_fdr_threshold) |>
    dplyr::group_by(.data$contrast) |>
    dplyr::slice_min(.data[[fdr_col]], n = n_labels) |>
    dplyr::ungroup() |>
    dplyr::mutate(label_hjust = ifelse(.data$NES < 0, 1.1, -0.1))

  ggplot2::ggplot(volcano_data, ggplot2::aes(x = .data$NES, y = .data$neg_log_fdr)) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$direction, alpha = .data$significant),
      size = 2
    ) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey30") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(label = .data[[item_col]], hjust = .data$label_hjust),
      size = 2.5,
      check_overlap = TRUE
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.16, 0.16))) +
    ggplot2::scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey50")) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    ggplot2::facet_wrap(~ .data$contrast, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 18, 5.5, 18)) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "NES",
      y = "-log10(FDR)",
      color = "Direction"
    )
}

#' Create NES heatmap for top items across contrasts
#'
#' @param data Data frame with columns: item (ID/kinase), NES, p.adjust/FDR, contrast
#' @param item_col Name of item column (default: "ID")
#' @param fdr_col Name of FDR column (default: "p.adjust")
#' @param fdr_filter FDR threshold for selecting top items (default: 0.15)
#' @param n_top Number of top items to show (default: 25)
#' @param item_label_col Optional column for shorter labels (default: NULL, uses item_col)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return ggplot object
#' @export
#' @examples
#' ptmsea_results <- data.frame(
#'   ID = rep(c("KINASE-PSP_PKACA", "KINASE-PSP_AKT1", "KINASE-PSP_MAPK1"), 3),
#'   NES = c(2.1, 1.5, -1.2, 1.8, 0.9, -0.8, 1.2, 1.1, -1.5),
#'   p.adjust = c(0.001, 0.02, 0.05, 0.005, 0.12, 0.18, 0.08, 0.09, 0.03),
#'   contrast = rep(c("Early", "Mid", "Late"), each = 3)
#' )
#' plot_enrichment_heatmap(ptmsea_results, n_top = 3, title = "PTMSEA Kinase Activity")
plot_enrichment_heatmap <- function(
  data,
  item_col = "ID",
  fdr_col = "p.adjust",
  fdr_filter = 0.15,
  n_top = 25,
  item_label_col = NULL,
  title = NULL,
  subtitle = NULL
) {
  top_items <- data |>
    dplyr::filter(.data[[fdr_col]] < fdr_filter) |>
    dplyr::group_by(.data[[item_col]]) |>
    dplyr::summarize(min_padj = min(.data[[fdr_col]]), .groups = "drop") |>
    dplyr::arrange(.data$min_padj) |>
    utils::head(n_top) |>
    dplyr::pull(.data[[item_col]])
  if (is.null(item_label_col)) {
    item_label_col <- item_col
  }
  heatmap_data <- data |>
    dplyr::filter(.data[[item_col]] %in% top_items) |>
    dplyr::mutate(
      item_label = .data[[item_label_col]],
      sig_label = dplyr::case_when(
        .data[[fdr_col]] < 0.01 ~ "***",
        .data[[fdr_col]] < 0.05 ~ "**",
        .data[[fdr_col]] < 0.1 ~ "*",
        TRUE ~ ""
      )
    )
  if (is.null(subtitle)) {
    subtitle <- "Top items (* p<0.1, ** p<0.05, *** p<0.01)"
  }

  ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(x = .data$contrast, y = reorder(.data$item_label, .data$NES), fill = .data$NES)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), color = "black", size = 4) +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "",
      y = "",
      fill = "NES"
    )
}

# One tidy table of clusterProfiler GSEA results, with the right columns and
# zero rows when nothing was tested.
extract_gsea_results <- function(results) {
  if (length(results) == 0) {
    return(dplyr::tibble(
      ID = character(),
      NES = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      setSize = integer(),
      contrast = character()
    ))
  }
  purrr::map_dfr(names(results), function(contrast) {
    res <- results[[contrast]]@result
    dplyr::tibble(
      ID = as.character(res[["ID"]]),
      NES = as.numeric(res[["NES"]]),
      pvalue = as.numeric(res[["pvalue"]]),
      p.adjust = as.numeric(res[["p.adjust"]]),
      setSize = as.integer(res[["setSize"]]),
      contrast = contrast
    )
  })
}
