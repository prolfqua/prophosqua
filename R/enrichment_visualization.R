#' Enrichment Visualization Functions
#'
#' Visualization functions for PTMSEA, KinaseLibrary GSEA, and MEA analyses.
#' These functions provide consistent plotting styles across enrichment analyses.
#'
#' @importFrom dplyr mutate filter arrange group_by summarize case_when pull slice_min
#' @importFrom ggplot2 ggplot aes geom_point geom_tile geom_text geom_vline
#'   geom_hline scale_color_gradient2 scale_fill_gradient2 scale_size_continuous
#'   scale_alpha_manual scale_color_manual facet_wrap theme_bw theme_minimal
#'   theme element_text labs
#' @importFrom stats reorder
#' @importFrom forcats fct_reorder
#' @importFrom rlang .data
#' @importFrom purrr map_dfr map_int
#' @name enrichment_visualization
NULL


#' Prepare enrichment data with common computed columns
#'
#' Adds neg_log_fdr, direction, and significant columns to enrichment results.
#' This is a shared helper used by plot functions and analyses.
#'
#' @param data Data frame with NES and FDR columns
#' @param fdr_col Name of FDR column (default: "FDR")
#' @param fdr_threshold FDR threshold for significance (default: 0.1)
#' @return Data frame with added neg_log_fdr, direction, significant columns
#' @export
#' @examples
#' # Basic usage with default FDR column
#' df <- data.frame(NES = c(1.5, -1.2, 0.5, -0.3), FDR = c(0.01, 0.08, 0.05, 0.25))
#' result <- prepare_enrichment_data(df)
#' result$direction   # "Up", "Down", "Up", "Down"
#' result$significant # TRUE, TRUE, TRUE, FALSE
#'
#' # With custom FDR column name (e.g., from clusterProfiler)
#' df2 <- data.frame(NES = c(2.0, -1.5), p.adjust = c(0.02, 0.12))
#' result2 <- prepare_enrichment_data(df2, fdr_col = "p.adjust")
#'
#' # With stricter FDR threshold
#' result3 <- prepare_enrichment_data(df, fdr_threshold = 0.05)
#' result3$significant # TRUE, FALSE, TRUE, FALSE
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


#' Create enrichment dotplot for top items by FDR
#'
#' @param data Data frame with columns: item (kinase/pathway), NES, p.adjust/FDR, contrast
#' @param item_col Name of item column (default: "kinase")
#' @param fdr_col Name of FDR column (default: "FDR")
#' @param n_top Number of top items to show (default: 30)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return ggplot object
#' @export
#' @examples
#' # Kinase enrichment results (e.g., from KinaseLibrary)
#' kinase_results <- data.frame(
#'   kinase = c("PKACA", "AKT1", "MAPK1", "CDK1", "CSNK2A1", "GSK3B"),
#'   NES = c(2.1, 1.8, -1.5, 1.2, -0.8, 0.3),
#'   FDR = c(0.001, 0.01, 0.02, 0.08, 0.15, 0.45)
#' )
#' plot_enrichment_dotplot(kinase_results, n_top = 6, title = "Kinase Enrichment")
#'
#' # Pathway enrichment with custom column names
#' pathway_results <- data.frame(
#'   pathway = c("Cell cycle", "Apoptosis", "MAPK signaling"),
#'   NES = c(1.9, -1.4, 1.1),
#'   p.adjust = c(0.005, 0.03, 0.09)
#' )
#' plot_enrichment_dotplot(pathway_results, item_col = "pathway",
#'                         fdr_col = "p.adjust", n_top = 3)
plot_enrichment_dotplot <- function(data, item_col = "kinase", fdr_col = "FDR",
                                    n_top = 30, title = NULL, subtitle = "Top 30 by FDR") {
  # Prepare data using shared helper
  plot_data <- prepare_enrichment_data(data, fdr_col, 0.1) |>
    dplyr::arrange(.data[[fdr_col]]) |>
    utils::head(n_top) |>
    dplyr::mutate(item = forcats::fct_reorder(.data[[item_col]], .data$NES))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$NES, y = .data$item)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$neg_log_fdr, color = .data$NES, alpha = .data$significant)) +
    ggplot2::scale_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0) +
    ggplot2::scale_size_continuous(name = "-log10(FDR)", range = c(2, 8)) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9)) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Normalized Enrichment Score (NES)",
      y = ""
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
#' # Multi-contrast kinase enrichment results
#' kinase_results <- data.frame(
#'   kinase = rep(c("PKACA", "AKT1", "MAPK1", "CDK1"), 2),
#'   NES = c(2.1, 1.5, -1.8, 0.5, 1.2, -0.8, 1.9, -1.1),
#'   FDR = c(0.001, 0.02, 0.01, 0.3, 0.05, 0.2, 0.008, 0.06),
#'   contrast = rep(c("Treatment_vs_Control", "Timepoint_vs_Baseline"), each = 4)
#' )
#' plot_enrichment_volcano(kinase_results, title = "Kinase Enrichment Volcano")
#'
#' # Customize labeling thresholds
#' plot_enrichment_volcano(kinase_results, fdr_threshold = 0.05,
#'                         label_fdr_threshold = 0.02, n_labels = 3)
plot_enrichment_volcano <- function(data, item_col = "kinase", fdr_col = "FDR",
                                    fdr_threshold = 0.1, label_fdr_threshold = 0.05,
                                    n_labels = 5, title = NULL, subtitle = NULL) {
  # Prepare data using shared helper
  volcano_data <- prepare_enrichment_data(data, fdr_col, fdr_threshold)

  # Default subtitle
  if (is.null(subtitle)) {
    subtitle <- paste0("Dashed line: FDR = ", fdr_threshold, "; Labels: top ", n_labels, " by FDR per contrast")
  }

  ggplot2::ggplot(volcano_data, ggplot2::aes(x = .data$NES, y = .data$neg_log_fdr)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$direction, alpha = .data$significant), size = 2) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey30") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::geom_text(
      data = volcano_data |>
        dplyr::filter(.data[[fdr_col]] < label_fdr_threshold) |>
        dplyr::group_by(.data$contrast) |>
        dplyr::slice_min(.data[[fdr_col]], n = n_labels),
      ggplot2::aes(label = .data[[item_col]]),
      size = 2.5,
      hjust = -0.1,
      check_overlap = TRUE
    ) +
    ggplot2::scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey50")) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    ggplot2::facet_wrap(~.data$contrast, scales = "free") +
    ggplot2::theme_bw() +
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
#' # PTMSEA results across multiple contrasts
#' ptmsea_results <- data.frame(
#'   ID = rep(c("KINASE-PSP_PKACA", "KINASE-PSP_AKT1", "KINASE-PSP_MAPK1"), 3),
#'   NES = c(2.1, 1.5, -1.2, 1.8, 0.9, -0.8, 1.2, 1.1, -1.5),
#'   p.adjust = c(0.001, 0.02, 0.05, 0.005, 0.12, 0.18, 0.08, 0.09, 0.03),
#'   contrast = rep(c("Early", "Mid", "Late"), each = 3)
#' )
#' plot_enrichment_heatmap(ptmsea_results, n_top = 3,
#'                         title = "PTMSEA Kinase Activity")
#'
#' # With custom item labels (shorter names)
#' ptmsea_results$short_name <- gsub("KINASE-PSP_", "", ptmsea_results$ID)
#' plot_enrichment_heatmap(ptmsea_results, item_label_col = "short_name")
plot_enrichment_heatmap <- function(data, item_col = "ID", fdr_col = "p.adjust",
                                    fdr_filter = 0.15, n_top = 25,
                                    item_label_col = NULL, title = NULL, subtitle = NULL) {
  # Find top items across all contrasts
  top_items <- data |>
    dplyr::filter(.data[[fdr_col]] < fdr_filter) |>
    dplyr::group_by(.data[[item_col]]) |>
    dplyr::summarize(min_padj = min(.data[[fdr_col]]), .groups = "drop") |>
    dplyr::arrange(.data$min_padj) |>
    utils::head(n_top) |>
    dplyr::pull(.data[[item_col]])

  # Use item_col for labels if no separate label column provided
  if (is.null(item_label_col)) {
    item_label_col <- item_col
  }

  # Prepare heatmap data
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

  # Default subtitle
  if (is.null(subtitle)) {
    subtitle <- "Top items (* p<0.1, ** p<0.05, *** p<0.01)"
  }

  ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data$contrast, y = reorder(.data$item_label, .data$NES), fill = .data$NES)) +
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
      x = "", y = "", fill = "NES"
    )
}


#' Export GSEA enrichment plots to PDF
#'
#' Exports clusterProfiler GSEA plots for significant gene sets to a PDF file.
#' Requires the enrichplot package.
#'
#' @param results Named list of clusterProfiler GSEA result objects
#' @param output_file Path to output PDF file
#' @param fdr_threshold FDR threshold for including plots (default: 0.25)
#' @param width PDF width in inches (default: 10)
#' @param height PDF height in inches (default: 6)
#' @param prefix_pattern Regex pattern to strip from pathway names for titles (default: NULL)
#' @return Number of plots exported
#' @export
#' @examples
#' \dontrun{
#' # Requires clusterProfiler GSEA results
#' export_gsea_plots_pdf(gsea_results, "enrichment_plots.pdf")
#' }
export_gsea_plots_pdf <- function(results, output_file, fdr_threshold = 0.25,
                                  width = 10, height = 6, prefix_pattern = NULL) {
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required for export_gsea_plots_pdf()")
  }

  grDevices::pdf(output_file, width = width, height = height)

  n_plots <- 0
  for (ct in names(results)) {
    res <- results[[ct]]
    genesets <- res@result |>
      dplyr::as_tibble() |>
      dplyr::filter(.data$p.adjust < fdr_threshold) |>
      dplyr::arrange(.data$pvalue) |>
      dplyr::pull(.data$ID)

    for (geneset in genesets) {
      row <- res@result |> dplyr::as_tibble() |> dplyr::filter(.data$ID == geneset)

      # Clean pathway name for title
      pathway_label <- geneset
      if (!is.null(prefix_pattern)) {
        pathway_label <- gsub(prefix_pattern, "", geneset)
      }

      nes_val <- round(row$NES, 2)
      fdr <- signif(row$p.adjust, 2)

      p <- enrichplot::gseaplot2(res, geneSetID = geneset,
                                 title = paste0(ct, ": ", pathway_label, " (NES=", nes_val, ", FDR=", fdr, ")"))
      print(p)
      n_plots <- n_plots + 1
    }
  }

  grDevices::dev.off()
  return(n_plots)
}


#' Create summary table for enrichment results
#'
#' Generates a summary table with counts of significant results at different FDR thresholds.
#'
#' @param results Named list of clusterProfiler GSEA result objects OR data frame
#' @param fdr_thresholds Vector of FDR thresholds to count (default: c(0.1, 0.05))
#' @return tibble with summary statistics
#' @export
#' @examples
#' # Summarize MEA/kinase enrichment data frame
#' mea_results <- data.frame(
#'   contrast = c(rep("Treatment_vs_Control", 50), rep("Drug_vs_Vehicle", 50)),
#'   kinase = paste0("Kinase", 1:100),
#'   NES = rnorm(100),
#'   FDR = runif(100, 0, 0.3)
#' )
#' summarize_enrichment_results(mea_results)
#' # Returns: contrast, total, FDR < 0.1, FDR < 0.05
summarize_enrichment_results <- function(results, fdr_thresholds = c(0.1, 0.05)) {
  if (is.data.frame(results)) {
    # Handle data frame input (MEA style)
    results |>
      dplyr::group_by(.data$contrast) |>
      dplyr::summarize(
        total = dplyr::n(),
        `FDR < 0.1` = sum(.data$FDR < 0.1, na.rm = TRUE),
        `FDR < 0.05` = sum(.data$FDR < 0.05, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    # Handle clusterProfiler results list
    dplyr::tibble(
      Contrast = names(results),
      `Total` = purrr::map_int(results, ~nrow(.x@result)),
      `FDR < 0.1` = purrr::map_int(results, ~sum(.x@result$p.adjust < 0.1, na.rm = TRUE)),
      `FDR < 0.05` = purrr::map_int(results, ~sum(.x@result$p.adjust < 0.05, na.rm = TRUE))
    )
  }
}


#' Extract results from clusterProfiler GSEA objects
#'
#' Converts a list of clusterProfiler GSEA result objects into a tidy tibble.
#'
#' @param results Named list of clusterProfiler GSEA result objects
#' @return tibble with columns: ID, NES, pvalue, p.adjust, setSize, contrast
#' @export
#' @examples
#' \dontrun{
#' # Requires clusterProfiler GSEA results
#' df <- extract_gsea_results(gsea_results)
#' }
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

  names(results) |>
    purrr::map_dfr(function(ct) {
      res <- results[[ct]]@result
      dplyr::tibble(
        ID = as.character(res[["ID"]]),
        NES = as.numeric(res[["NES"]]),
        pvalue = as.numeric(res[["pvalue"]]),
        p.adjust = as.numeric(res[["p.adjust"]]),
        setSize = as.integer(res[["setSize"]]),
        contrast = ct
      )
    })
}
