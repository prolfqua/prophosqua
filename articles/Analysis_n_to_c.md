# PTM Analysis - N-to-C Plots

## Introduction

This document generates N-to-C plots for **DPA** analysis.

``` r

desc <- switch(params$analysis_type,
  "dpa" = "**DPA (Differential PTM Abundance)**: Raw PTM signal changes without protein abundance correction. Shows both site-level and protein-level fold changes.",
  "dpu" = "**DPU (Differential PTM Usage)**: Protein-normalized PTM changes (model-first approach). Reflects genuine modification stoichiometry changes.",
  "cf" = "**CF (CorrectFirst)**: Alternative protein-correction approach where abundances are normalized before statistical modeling."
)
cat(desc)
```

**DPA (Differential PTM Abundance)**: Raw PTM signal changes without
protein abundance correction. Shows both site-level and protein-level
fold changes.

## Data Loading

``` r

if (pipeline_mode) {
  message("Loading data from: ", params$xlsx_file, " (sheet: ", params$sheet, ")")
  data <- readxl::read_xlsx(params$xlsx_file, sheet = params$sheet)
  output_dir <- if (!is.null(params$output_dir)) params$output_dir else dirname(params$xlsx_file)
} else {
  # Vignette mode: use example data
  data("combined_test_diff_example", package = "prophosqua")
  data <- combined_test_diff_example
  output_dir <- tempdir()
  max_fig <- 5
  message("Using example data from prophosqua package")
}

data_info <- tibble(
  Property = c("Mode", "Sheet", "Analysis Type", "Rows", "Contrasts"),
  Value = c(
    if (pipeline_mode) basename(params$xlsx_file) else "Example data",
    params$sheet, toupper(params$analysis_type),
    nrow(data), paste(unique(data$contrast), collapse = ", ")
  )
)
knitr::kable(
  data_info,
  caption = paste0(
    "Input of this N-to-C report: source workbook, sheet, analysis type, number of ",
    "site-by-contrast rows read and the contrasts they cover. Every row is one phosphosite ",
    "in one contrast, before the significance filter that selects the plotted proteins."
  )
)
```

| Property | Value |
|:---|:---|
| Mode | Example data |
| Sheet | DPA |
| Analysis Type | DPA |
| Rows | 105824 |
| Contrasts | KO_vs_WT, KO_vs_WT_at_Early, KO_vs_WT_at_Late, KO_vs_WT_at_Uninfect |

Input of this N-to-C report: source workbook, sheet, analysis type,
number of site-by-contrast rows read and the contrasts they cover. Every
row is one phosphosite in one contrast, before the significance filter
that selects the plotted proteins. {.table}

``` r

# Prepare data for plotting using shared function
plot_data <- prepare_ntoc_data(data, params$analysis_type)

all_contrasts <- unique(plot_data$contrast)
cat("Found", length(all_contrasts), "contrasts:", paste(all_contrasts, collapse = ", "), "\n")
```

    ## Found 4 contrasts: KO_vs_WT, KO_vs_WT_at_Early, KO_vs_WT_at_Late, KO_vs_WT_at_Uninfect

## Generate N-to-C Plots

**Filtering criteria:** Proteins are included if they contain at least
one PTM site with FDR \< 0.05 and \|log2 fold change\| \> 0.6.

``` r

if (params$analysis_type == "dpa") {
  # DPA: Use expression plotting (shows protein + site)
  plot_result <- n_to_c_expression_multicontrast(
    plot_data,
    FDR_threshold = params$fdr,
    fc_threshold = params$fc,
    max_plots = max_fig
  )
} else {
  # DPU/CF: Use usage plotting (shows site only)
  plot_result <- n_to_c_usage_multicontrast(
    plot_data,
    FDR_threshold = params$fdr,
    fc_threshold = params$fc,
    max_plots = max_fig
  )
}

n_plots <- nrow(plot_result)
stopifnot("No significant sites found. Adjust fdr/log2fc thresholds in config." = n_plots > 0)
cat("Generated", n_plots, "multi-contrast plots\n")
```

### Example Plots

Showing example plots in the HTML report: one protein WITH protein-level
data and one WITHOUT (if available).

``` r

# Pick the examples first, so that each one can be given its own figure caption.
matched_proteins <- character(0)
protein_note <- NULL

if (params$analysis_type == "dpa") {
  # DPA: show one protein with and one without matched protein-level data
  matched_proteins <- plot_data |>
    dplyr::filter(!is.na(diff.protein)) |>
    dplyr::pull(protein_Id) |>
    unique()

  unmatched_proteins <- plot_data |>
    dplyr::filter(is.na(diff.protein)) |>
    dplyr::pull(protein_Id) |>
    unique()

  matched_idx <- which(plot_result$protein_Id %in% matched_proteins)
  unmatched_idx <- which(plot_result$protein_Id %in% unmatched_proteins)

  example_idx <- c(
    if (length(matched_idx) > 0) matched_idx[1],
    if (length(unmatched_idx) > 0) unmatched_idx[1]
  )
  if (length(example_idx) == 0) {
    example_idx <- seq_len(min(2, n_plots))
  }

  protein_note <- paste0(
    "**Note:** ", length(matched_idx), " proteins have protein-level data, ",
    length(unmatched_idx), " proteins do not."
  )
} else {
  # DPU/CF: no protein-level track to compare against, show the first two
  example_idx <- seq_len(min(2, n_plots))
}

example_has_protein <- plot_result$protein_Id[example_idx] %in% matched_proteins
example_labels <- ifelse(example_has_protein, " (with protein data)", "")
if (params$analysis_type == "dpa") {
  example_labels <- ifelse(example_has_protein, " (with protein data)", " (NO protein data)")
}

example_caps <- paste0(
  "Position-resolved ", toupper(params$analysis_type), " changes along protein ",
  plot_result$protein_Id[example_idx], example_labels,
  ", one panel per contrast. The x-axis is the residue position from the N- to the ",
  "C-terminus of the protein and the y-axis is the site-level log2 fold change of that ",
  "contrast; each vertical stick is one PTM site, coloured by modified residue (S blue, ",
  "T green, Y brown, not localized pink) and dashed when the value was imputed. Asterisks ",
  "give the site FDR (** < 0.05, * < 0.2), and the yellow rectangle spans the protein from ",
  "N to C at the protein-level log2 fold change, so sticks reaching beyond it change more ",
  "than the protein itself. Proteins are shown when at least one site has FDR < ",
  params$fdr, " and |log2FC| > ", params$fc, "."
)
```

**Note:** 5 proteins have protein-level data, 0 proteins do not.

``` r

if (n_plots > 0) {
  for (i in example_idx) {
    cat("\n\n### ", plot_result$protein_Id[i], example_labels[match(i, example_idx)], "\n\n", sep = "")
    print(plot_result$plot[[i]])
    cat("\n\n")
  }
} else {
  cat("No significant proteins found with the current filtering criteria.\n")
}
```

#### sp\|A0A140LIF8\|IRGM2_MOUSE (with protein data)

![Position-resolved DPA changes along protein
sp\|A0A140LIF8\|IRGM2_MOUSE (with protein data), one panel per contrast.
The x-axis is the residue position from the N- to the C-terminus of the
protein and the y-axis is the site-level log2 fold change of that
contrast; each vertical stick is one PTM site, coloured by modified
residue (S blue, T green, Y brown, not localized pink) and dashed when
the value was imputed. Asterisks give the site FDR (\*\* \< 0.05, \* \<
0.2), and the yellow rectangle spans the protein from N to C at the
protein-level log2 fold change, so sticks reaching beyond it change more
than the protein itself. Proteins are shown when at least one site has
FDR \< 0.05 and \|log2FC\| \>
0.6.](Analysis_n_to_c_files/figure-html/displayExamples-1.png)

Position-resolved DPA changes along protein sp\|A0A140LIF8\|IRGM2_MOUSE
(with protein data), one panel per contrast. The x-axis is the residue
position from the N- to the C-terminus of the protein and the y-axis is
the site-level log2 fold change of that contrast; each vertical stick is
one PTM site, coloured by modified residue (S blue, T green, Y brown,
not localized pink) and dashed when the value was imputed. Asterisks
give the site FDR (\*\* \< 0.05, \* \< 0.2), and the yellow rectangle
spans the protein from N to C at the protein-level log2 fold change, so
sticks reaching beyond it change more than the protein itself. Proteins
are shown when at least one site has FDR \< 0.05 and \|log2FC\| \> 0.6.

## Export Plots

``` r

# Use explicit output_dir if provided, otherwise use xlsx directory
pdf_dir <- output_dir

if (!dir.exists(pdf_dir)) {
  dir.create(pdf_dir, recursive = TRUE)
}

# Determine output filename based on analysis type
pdf_name <- switch(params$analysis_type,
  "dpa" = "Site_differential_Expression_MultiContrast.pdf",
  "dpu" = "Site_differential_UsageChange_MultiContrast.pdf",
  "cf"  = "Site_differential_CorrectFirst_MultiContrast.pdf"
)

pdf_path <- file.path(pdf_dir, pdf_name)

if (n_plots > 0) {
  pdf(pdf_path, width = 14, height = 10)
  for (i in seq_len(n_plots)) {
    print(plot_result$plot[[i]])
  }
  dev.off()
  cat("Exported", n_plots, "plots to:", pdf_path, "\n")
} else {
  cat("No plots to export.\n")
}
```

``` r

message("Vignette mode: PDF export skipped.")
pdf_path <- NULL
```

## Results Summary

``` r

summary_info <- tibble(
  Metric = c("Analysis Type", "Total Proteins", "Plots Generated",
             "Shown in HTML", "FDR Threshold", "FC Threshold"),
  Value = c(toupper(params$analysis_type), n_distinct(plot_data$protein_Id),
            n_plots, min(2, n_plots), params$fdr, params$fc)
)
knitr::kable(
  summary_info,
  caption = paste0(
    "Coverage of this report: how many proteins carry data, how many passed the significance ",
    "filter and were plotted, and how many of those plots are embedded here. Total Proteins ",
    "counts all proteins in the sheet, Plots Generated the proteins with at least one site ",
    "below the FDR threshold and above the fold-change threshold (capped by the max_fig ",
    "setting), and the remaining plots are available in the exported PDF."
  )
)
```

| Metric          | Value |
|:----------------|:------|
| Analysis Type   | DPA   |
| Total Proteins  | 9581  |
| Plots Generated | 5     |
| Shown in HTML   | 2     |
| FDR Threshold   | 0.05  |
| FC Threshold    | 0.6   |

Coverage of this report: how many proteins carry data, how many passed
the significance filter and were plotted, and how many of those plots
are embedded here. Total Proteins counts all proteins in the sheet,
Plots Generated the proteins with at least one site below the FDR
threshold and above the fold-change threshold (capped by the max_fig
setting), and the remaining plots are available in the exported PDF.
{.table}

``` r

if (pipeline_mode && n_plots > 0) {
  cat("\nPDF exported to:", pdf_path, "\n")
}
```

## Session Info

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] patchwork_1.3.2  readxl_1.5.0     dplyr_1.2.1      prophosqua_0.3.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       jsonlite_2.0.0     compiler_4.6.1     tidyselect_1.2.1  
    ##  [5] optparse_1.8.2     jquerylib_0.1.4    systemfonts_1.3.2  scales_1.4.0      
    ##  [9] textshaping_1.0.5  yaml_2.3.12        fastmap_1.2.0      ggplot2_4.0.3     
    ## [13] R6_2.6.1           labeling_0.4.3     generics_0.1.4     knitr_1.51        
    ## [17] htmlwidgets_1.6.4  forcats_1.0.1      tibble_3.3.1       bookdown_0.47     
    ## [21] desc_1.4.3         RColorBrewer_1.1-3 bslib_0.12.0       pillar_1.11.1     
    ## [25] rlang_1.3.0        cachem_1.1.0       xfun_0.60          S7_0.2.2          
    ## [29] fs_2.1.0           sass_0.4.10        otel_0.2.0         cli_3.6.6         
    ## [33] withr_3.0.3        pkgdown_2.2.1      magrittr_2.0.5     digest_0.6.39     
    ## [37] grid_4.6.1         lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5    
    ## [41] glue_1.8.1         cellranger_1.1.0   farver_2.1.2       ggseqlogo_0.2.2   
    ## [45] ragg_1.5.2         purrr_1.2.2        rmarkdown_2.31     tools_4.6.1       
    ## [49] pkgconfig_2.0.3    htmltools_0.5.9
