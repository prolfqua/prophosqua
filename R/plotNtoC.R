#' @importFrom rlang .data
NULL

get_significance <- function(fdr, thr_a = 0.05, thr_b = 0.2) {
  if (is.na(fdr)) {
    return("")
  } else if (fdr < thr_a) {
    return("**")
  } else if (fdr < thr_b) {
    return("*")
  } else {
    return("")
  }
}

#' Prepare data for N-to-C plotting
#' @param poi_matrix_min data.frame with phosphorylation data
#' @export
#' @return data.frame with prepared data for plotting
#' @keywords internal
#' @examples
#' # example code
#' data(exampleN_C_dat)
#' poi_matrix_min <- prepare_n_to_c_data(exampleN_C_dat)
#'
prepare_n_to_c_data <- function(poi_matrix_min, model_site = "model_site") {
  # Add imputation status based on model name
  poi_matrix_min$imputation_status <- ifelse(grepl("imputed", poi_matrix_min[[model_site]]), "imputed", "observed")

  # Ensure numeric types for position columns
  class(poi_matrix_min[["startModSite"]]) <- "numeric"
  class(poi_matrix_min[["endModSite"]]) <- "numeric"

  # Calculate positions and handle non-localized sites
  poi_matrix_min <- poi_matrix_min |>
    dplyr::mutate(
      posInProtein = ifelse(
        .data$AllLocalized,
        .data$posInProtein,
        as.integer(.data$startModSite + .data$endModSite) / 2
      )
    ) |>
    dplyr::mutate(
      modAA = ifelse(.data$AllLocalized, .data$modAA, "NotLoc")
    )

  return(poi_matrix_min)
}

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
#' # Prepare data for plotting
#' poi_matrix_min <- prepare_n_to_c_data(exampleN_C_dat)
#'
#' n_to_c_plot(subset(poi_matrix_min, protein_Id == "A0A1I9LPZ1"), "A0A1I9LPZ1", 2160, "H1FC")
#' n_to_c_plot(subset(poi_matrix_min, protein_Id == "A0A178US29"), "A0A178US29", 806, "H1FC")
n_to_c_plot <- function(
  poi_matrix_min,
  protein_name,
  prot_length,
  contrast,
  thr_a = 0.05,
  thr_b = 0.2,
  color_protein = "yellow"
) {
  # Validate required columns
  required_cols <- c(
    "diff.protein", "diff.site", "FDR.site", "posInProtein",
    "modAA", "imputation_status"
  )
  missing_cols <- setdiff(required_cols, colnames(poi_matrix_min))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in poi_matrix_min: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  mean_diff_prot <- mean(poi_matrix_min$diff.protein, na.rm = TRUE)
  poi_matrix_min$significance <- sapply(poi_matrix_min$FDR.site, get_significance, thr_a, thr_b)

  plot_title <- paste0(
    "Prot : ", protein_name,
    "; length: ", prot_length, "; # sites:",
    nrow(poi_matrix_min), "; # not localized sites:", sum(poi_matrix_min$modAA == "NotLoc")
  )


  p <- ggplot(data = poi_matrix_min) +
    geom_segment(aes(
      x = .data$posInProtein, xend = .data$posInProtein, y = 0, yend = .data$diff.site, color = .data$modAA,
      linetype = .data$imputation_status
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
      xmax = .data$xmax, ymin = .data$ymin,
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
#' n_c_integrated_df$imputation_status <- "observed"
#' n_to_c_plot_integrated(n_c_integrated_df, "A0A1I9LT44", 539, "WTFC")
n_to_c_plot_integrated <- function(
  poi_matrix_min,
  protein_name,
  prot_length,
  contrast,
  thr_a = 0.05,
  thr_b = 0.2,
  color_protein = "yellow"
) {
  # Validate required columns
  required_cols <- c("diff_diff", "FDR_I", "posInProtein", "modAA")
  missing_cols <- setdiff(required_cols, colnames(poi_matrix_min))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in poi_matrix_min: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  poi_matrix_min$significance <- sapply(poi_matrix_min$FDR_I, get_significance, thr_a, thr_b)

  plot_title <- paste0(
    "Prot : ", protein_name,
    "; length: ",
    prot_length,
    "; # sites:", nrow(poi_matrix_min),
    "; # not localized sites:",
    sum(poi_matrix_min$modAA == "NotLoc")
  )


  mean_diff_prot <- 0
  p <- ggplot(data = poi_matrix_min) +
    geom_segment(aes(
      x = .data$posInProtein,
      xend = .data$posInProtein, y = 0, yend = .data$diff_diff,
      color = .data$modAA,
      linetype = .data$imputation_status
    )) +
    scale_linetype_manual(values = c("imputed" = "dashed", "observed" = "solid")) +
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

#' Calculate optimal grid layout for multi-panel plots
#' @param n_panels Number of panels to arrange
#' @return List with ncol and nrow for grid layout
#' @keywords internal
#' @keywords internal
calculate_optimal_layout <- function(n_panels) {
  if (n_panels == 1) {
    return(list(ncol = 1, nrow = 1))
  }

  # Calculate square root and round up for columns
  # Prefer wider layouts (more columns than rows) for better readability
  ncol <- ceiling(sqrt(n_panels))
  nrow <- ceiling(n_panels / ncol)

  return(list(ncol = ncol, nrow = nrow))
}

#' Plot protein and PTM expression
#' @param combined_site_prot_long data frame with combined site and protein data
#' @param contrast_name name of the contrast to plot
#' @param FDR_threshold FDR threshold for filtering significant sites (default 0.05)
#' @param fc_threshold Fold change threshold for filtering (default 0)
#' @param impute_flag Flag for imputed values (default "Imputed_Mean_moderated")
#' @return data frame with protein_Id, protein_length, contrast, data, and plot
#' @export
n_to_c_expression <- function(
  combined_site_prot_long,
  contrast_name,
  FDR_threshold = 0.05,
  fc_threshold = 0,
  impute_flag = "Imputed_Mean_moderated"
) {
  if (!contrast_name %in% unique(combined_site_prot_long$contrast)) {
    stop(contrast_name, " not in ", paste0(unique(combined_site_prot_long$contrast)))
  }
  combined_site_prot_long <- combined_site_prot_long |>
    dplyr::filter(contrast == contrast_name)
  combined_site_prot_long <- combined_site_prot_long |>
    dplyr::group_by(protein_Id) |>
    dplyr::filter(any(FDR.site < FDR_threshold & abs(diff.site) > fc_threshold)) |>
    dplyr::ungroup()

  required_cols <- c(
    "protein_Id",
    "contrast",
    "protein_length",
    "site", "diff.protein", "diff.site", "FDR.site",
    "posInProtein", "modAA", "imputation_status"
  )

  plot_data <- combined_site_prot_long |>
    dplyr::mutate(imputation_status = dplyr::case_when(modelName.site == impute_flag ~ "imputed", TRUE ~ "observed")) |>
    dplyr::select(dplyr::all_of(required_cols)) |>
    dplyr::group_by(.data$protein_Id, .data$contrast, .data$protein_length) |>
    tidyr::nest()

  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))

  pb <- txtProgressBar(min = 0, max = nrow(plot_data), style = 3)
  for (i in seq_len(nrow(plot_data))) {
    setTxtProgressBar(pb, i)
    plot_data$plot[[i]] <- n_to_c_plot(
      plot_data$data[[i]],
      plot_data$protein_Id[[i]],
      plot_data$protein_length[[i]],
      plot_data$contrast[[i]]
    )
  }
  close(pb)

  return(plot_data = plot_data)
}


#' Plot PTM occupation, use after adjusting for total proteome
#' @param data_combined_diff data frame with combined differential analysis results
#' @param contrast_name name of the contrast to plot
#' @param FDR_threshold FDR threshold for filtering significant sites (default 0.05)
#' @param fc_threshold Fold change threshold for filtering (default 0)
#' @param impute_flag Flag for imputed values (default "Imputed_Mean_moderated")
#' @param protein_Id Column name for protein identifier (default "protein_Id")
#' @return data frame with protein_Id, protein_length, contrast, data, and plot
#' @export
n_to_c_usage <- function(
  data_combined_diff,
  contrast_name,
  FDR_threshold = 0.05,
  fc_threshold = 0,
  impute_flag = "Imputed_Mean_moderated",
  protein_Id = "protein_Id"
) {
  data_combined_diff <- data_combined_diff |>
    dplyr::filter(contrast == contrast_name)
  data_combined_diff <- data_combined_diff |>
    dplyr::group_by(!!dplyr::sym(protein_Id)) |>
    dplyr::filter(any(FDR_I < FDR_threshold & abs(diff_diff) > fc_threshold)) |>
    dplyr::ungroup()

  required_columns <- c(
    "protein_Id",
    "contrast",
    "protein_length",
    "site",
    "diff_diff",
    "FDR_I",
    "posInProtein",
    "modAA",
    "imputation_status"
  )
  data_combined_diff$imputation_status <- ifelse(data_combined_diff$modelName.site == "Imputed_Mean_moderated" | data_combined_diff$modelName.protein == "Imputed_Mean_moderated", "imputed", "observed")

  # Create integrated N-to-C plots
  # Prepare data for plotting
  plot_data <- data_combined_diff |>
    dplyr::select(dplyr::all_of(required_columns)) |>
    dplyr::group_by(.data$protein_Id, .data$contrast, .data$protein_length) |>
    tidyr::nest()

  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))


  # Create plots
  message("Creating plots...")
  pb <- txtProgressBar(min = 0, max = nrow(plot_data), style = 3)
  for (i in seq_len(nrow(plot_data))) {
    setTxtProgressBar(pb, i)
    plot_data$plot[[i]] <- n_to_c_plot_integrated(
      plot_data$data[[i]],
      plot_data$protein_Id[[i]],
      plot_data$protein_length[[i]],
      plot_data$contrast[[i]]
    )
  }
  close(pb)
  return(plot_data)
}

#' Plot protein and PTM expression across all contrasts in multi-panel layout
#' @param combined_site_prot_long data frame with combined site and protein data
#' @param FDR_threshold FDR threshold for filtering significant sites (default 0.05)
#' @param fc_threshold Fold change threshold for filtering (default 0)
#' @param impute_flag Flag for imputed values (default "Imputed_Mean_moderated")
#' @param max_plots Maximum number of plots to generate (default NULL = no limit).
#'   If specified, only the first max_plots proteins will be plotted.
#' @param include_proteins Character vector of protein IDs to always include in the output,
#'   regardless of max_plots limit (default NULL). These proteins will be added even if
#'   they exceed the max_plots limit.
#' @return data frame with protein_Id, protein_length, n_contrasts, and multi-panel plot
#' @export
#' @examples
#' # data(combined_site_prot_data)
#' # result <- n_to_c_expression_multicontrast(combined_site_prot_data,
#' #   max_plots = 50, include_proteins = c("Q64337"))
#' # print(result$plot[[1]])  # Display first protein's multi-contrast plot
n_to_c_expression_multicontrast <- function(
  combined_site_prot_long,
  FDR_threshold = 0.05,
  fc_threshold = 0,
  impute_flag = "Imputed_Mean_moderated",
  max_plots = NULL,
  include_proteins = NULL
) {
  # Get all unique contrasts
  all_contrasts <- unique(combined_site_prot_long$contrast)
  n_contrasts <- length(all_contrasts)

  message("Found ", n_contrasts, " contrasts: ", paste(all_contrasts, collapse = ", "))

  # Filter proteins that have at least one significant site in ANY contrast
  significant_proteins <- combined_site_prot_long |>
    dplyr::group_by(protein_Id) |>
    dplyr::filter(any(FDR.site < FDR_threshold & abs(diff.site) > fc_threshold)) |>
    dplyr::ungroup()

  # Get unique proteins with their lengths
  proteins_to_plot <- significant_proteins |>
    dplyr::select(protein_Id, protein_length) |>
    dplyr::distinct()

  # Limit number of plots if max_plots is specified
  n_total <- nrow(proteins_to_plot)
  if (!is.null(max_plots) && n_total > max_plots) {
    # Keep first max_plots proteins
    proteins_limited <- proteins_to_plot[seq_len(max_plots), ]

    # Add include_proteins if specified and not already in the set
    if (!is.null(include_proteins)) {
      missing_proteins <- setdiff(include_proteins, proteins_limited$protein_Id)
      if (length(missing_proteins) > 0) {
        extra_proteins <- proteins_to_plot |>
          dplyr::filter(protein_Id %in% missing_proteins)
        proteins_limited <- dplyr::bind_rows(proteins_limited, extra_proteins)
        message(
          "Added ", nrow(extra_proteins), " requested protein(s): ",
          paste(missing_proteins, collapse = ", ")
        )
      }
    }
    proteins_to_plot <- proteins_limited
    message("Limiting to ", nrow(proteins_to_plot), " plots (out of ", n_total, " significant proteins)")
  }

  message("Creating multi-contrast plots for ", nrow(proteins_to_plot), " proteins...")

  # Prepare required columns
  required_cols <- c(
    "protein_Id", "contrast", "protein_length", "site",
    "diff.protein", "diff.site", "FDR.site",
    "posInProtein", "modAA", "imputation_status"
  )

  # Add imputation status
  combined_site_prot_long <- combined_site_prot_long |>
    dplyr::mutate(imputation_status = dplyr::case_when(
      modelName.site == impute_flag ~ "imputed",
      TRUE ~ "observed"
    ))

  # Calculate optimal layout
  layout <- calculate_optimal_layout(n_contrasts)

  # Initialize results
  plot_data <- proteins_to_plot
  plot_data$n_contrasts <- n_contrasts
  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))

  # Create multi-panel plots for each protein
  pb <- txtProgressBar(min = 0, max = nrow(plot_data), style = 3)

  for (i in seq_len(nrow(plot_data))) {
    setTxtProgressBar(pb, i)

    current_protein <- plot_data$protein_Id[[i]]
    current_length <- plot_data$protein_length[[i]]

    # Get data for this protein across all contrasts
    protein_data <- combined_site_prot_long |>
      dplyr::filter(protein_Id == current_protein)

    # Create individual plots for each contrast
    contrast_plots <- list()

    for (j in seq_along(all_contrasts)) {
      contrast <- all_contrasts[j]

      # Get data for this specific contrast
      contrast_data <- protein_data |>
        dplyr::filter(contrast == !!contrast) |>
        dplyr::select(dplyr::all_of(required_cols))

      # Create plot if data exists
      if (nrow(contrast_data) > 0) {
        contrast_plots[[j]] <- n_to_c_plot(
          contrast_data,
          current_protein,
          current_length,
          contrast
        )
      } else {
        # Create empty placeholder plot
        contrast_plots[[j]] <- ggplot2::ggplot() +
          ggplot2::annotate("text",
            x = 0.5, y = 0.5,
            label = paste0("No data for\n", contrast),
            size = 5, color = "gray50"
          ) +
          ggplot2::theme_void() +
          ggplot2::labs(title = contrast)
      }
    }

    # Combine plots using patchwork
    combined_plot <- patchwork::wrap_plots(
      contrast_plots,
      ncol = layout$ncol,
      nrow = layout$nrow
    ) +
      patchwork::plot_annotation(
        title = paste0("Protein: ", current_protein, " (length: ", current_length, ")"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
      )

    plot_data$plot[[i]] <- combined_plot
  }

  close(pb)
  message("Created ", nrow(plot_data), " multi-contrast plots")

  return(plot_data)
}

#' Plot PTM usage across all contrasts in multi-panel layout
#' @param data_combined_diff data frame with combined differential analysis results
#' @param FDR_threshold FDR threshold for filtering significant sites (default 0.05)
#' @param fc_threshold Fold change threshold for filtering (default 0)
#' @param impute_flag Flag for imputed values (default "Imputed_Mean_moderated")
#' @param protein_Id Column name for protein identifier (default "protein_Id")
#' @param max_plots Maximum number of plots to generate (default NULL = no limit).
#'   If specified, only the first max_plots proteins will be plotted.
#' @param include_proteins Character vector of protein IDs to always include in the output,
#'   regardless of max_plots limit (default NULL). These proteins will be added even if
#'   they exceed the max_plots limit.
#' @return data frame with protein_Id, protein_length, n_contrasts, and multi-panel plot
#' @export
#' @examples
#' # data(combined_diff_data)
#' # result <- n_to_c_usage_multicontrast(combined_diff_data,
#' #   max_plots = 50, include_proteins = c("Q64337"))
#' # print(result$plot[[1]])  # Display first protein's multi-contrast plot
n_to_c_usage_multicontrast <- function(
  data_combined_diff,
  FDR_threshold = 0.05,
  fc_threshold = 0,
  impute_flag = "Imputed_Mean_moderated",
  protein_Id = "protein_Id",
  max_plots = NULL,
  include_proteins = NULL
) {
  # Get all unique contrasts
  all_contrasts <- unique(data_combined_diff$contrast)
  n_contrasts <- length(all_contrasts)

  message("Found ", n_contrasts, " contrasts: ", paste(all_contrasts, collapse = ", "))

  # Filter proteins that have at least one significant site in ANY contrast
  significant_proteins <- data_combined_diff |>
    dplyr::group_by(!!dplyr::sym(protein_Id)) |>
    dplyr::filter(any(FDR_I < FDR_threshold & abs(diff_diff) > fc_threshold)) |>
    dplyr::ungroup()

  # Get unique proteins with their lengths
  proteins_to_plot <- significant_proteins |>
    dplyr::select(dplyr::all_of(c(protein_Id, "protein_length"))) |>
    dplyr::distinct()

  # Limit number of plots if max_plots is specified
  n_total <- nrow(proteins_to_plot)
  if (!is.null(max_plots) && n_total > max_plots) {
    # Keep first max_plots proteins
    proteins_limited <- proteins_to_plot[seq_len(max_plots), ]

    # Add include_proteins if specified and not already in the set
    if (!is.null(include_proteins)) {
      missing_proteins <- setdiff(include_proteins, proteins_limited[[protein_Id]])
      if (length(missing_proteins) > 0) {
        extra_proteins <- proteins_to_plot |>
          dplyr::filter(!!dplyr::sym(protein_Id) %in% missing_proteins)
        proteins_limited <- dplyr::bind_rows(proteins_limited, extra_proteins)
        message(
          "Added ", nrow(extra_proteins), " requested protein(s): ",
          paste(missing_proteins, collapse = ", ")
        )
      }
    }
    proteins_to_plot <- proteins_limited
    message("Limiting to ", nrow(proteins_to_plot), " plots (out of ", n_total, " significant proteins)")
  }

  message("Creating multi-contrast plots for ", nrow(proteins_to_plot), " proteins...")

  # Prepare required columns
  required_columns <- c(
    "protein_Id", "contrast", "protein_length", "site",
    "diff_diff", "FDR_I", "posInProtein", "modAA",
    "imputation_status"
  )

  # Add imputation status
  data_combined_diff$imputation_status <- ifelse(
    data_combined_diff$modelName.site == impute_flag |
      data_combined_diff$modelName.protein == impute_flag,
    "imputed", "observed"
  )

  # Calculate optimal layout
  layout <- calculate_optimal_layout(n_contrasts)

  # Initialize results
  plot_data <- proteins_to_plot
  plot_data$n_contrasts <- n_contrasts
  plot_data$plot <- vector(mode = "list", length = nrow(plot_data))

  # Create multi-panel plots for each protein
  pb <- txtProgressBar(min = 0, max = nrow(plot_data), style = 3)

  for (i in seq_len(nrow(plot_data))) {
    setTxtProgressBar(pb, i)

    current_protein <- plot_data[[protein_Id]][[i]]
    current_length <- plot_data$protein_length[[i]]

    # Get data for this protein across all contrasts
    protein_data <- data_combined_diff |>
      dplyr::filter(!!dplyr::sym(protein_Id) == current_protein)

    # Create individual plots for each contrast
    contrast_plots <- list()

    for (j in seq_along(all_contrasts)) {
      contrast <- all_contrasts[j]

      # Get data for this specific contrast
      contrast_data <- protein_data |>
        dplyr::filter(contrast == !!contrast) |>
        dplyr::select(dplyr::all_of(required_columns))

      # Create plot if data exists
      if (nrow(contrast_data) > 0) {
        contrast_plots[[j]] <- n_to_c_plot_integrated(
          contrast_data,
          current_protein,
          current_length,
          contrast
        )
      } else {
        # Create empty placeholder plot
        contrast_plots[[j]] <- ggplot2::ggplot() +
          ggplot2::annotate("text",
            x = 0.5, y = 0.5,
            label = paste0("No data for\n", contrast),
            size = 5, color = "gray50"
          ) +
          ggplot2::theme_void() +
          ggplot2::labs(title = contrast)
      }
    }

    # Combine plots using patchwork
    combined_plot <- patchwork::wrap_plots(
      contrast_plots,
      ncol = layout$ncol,
      nrow = layout$nrow
    ) +
      patchwork::plot_annotation(
        title = paste0("Protein: ", current_protein, " (length: ", current_length, ")"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
      )

    plot_data$plot[[i]] <- combined_plot
  }

  close(pb)
  message("Created ", nrow(plot_data), " multi-contrast plots")

  return(plot_data)
}
