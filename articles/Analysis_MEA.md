# Kinase Activity (Kinase Library + MEA)

## Overview

**Motif Enrichment Analysis (MEA)** asks, for every kinase, whether the
phosphosites that carry that kinase’s preferred sequence motif change
*coordinately* in a contrast. It therefore reports inferred kinase
*activity* changes in the **DPA** results — not kinase abundance, which
a phosphoproteomics experiment does not measure.

### How the numbers in this report are produced

| Step | What happens | Software |
|:---|:---|:---|
| 1\. Rank the sites | All quantified sites of a contrast are sorted by their site-level moderated t-statistic (`statistic.site`: most increased first, most decreased last) and written to a `.rnk` file, one per contrast. The sites, the sheet and the ranking statistic come from the combined `PTM_results.xlsx` that PTM-SEA and the KinaseLib GSEA also read, so the three enrichment analyses rank identically. For DPU that statistic is the site-versus-protein interaction, which the combined workbook stores under the same name | `prep_kinaselib.R` of the PTM pipeline |
| 2\. Score the motifs | The ±7-residue sequence window of every site is scored against the position-specific scoring matrix of every kinase, and each score is converted into a percentile relative to a reference phosphoproteome | [kinase-library](https://github.com/TheKinaseLibrary/kinase-library) (Python) |
| 3\. Build substrate sets | Per kinase, all sites scoring above the configured percentile threshold (95th percentile by default) form its predicted substrate set — one set per kinase (311 Ser/Thr, 78 Tyr kinases in Kinase Library v1.2.0) | kinase-library |
| 4\. Test for enrichment | Pre-ranked GSEA of those sets against the ranked site list, giving ES, NES, a permutation p-value (1000 permutations by default) and a Benjamini-Hochberg FDR per kinase | [GSEApy](https://gseapy.readthedocs.io) `prerank`, called by kinase-library |
| 5\. This report | Reads the per-contrast `mea_*.csv` files and produces the summaries, figures and Excel export below | R: `prophosqua`, `dplyr`, `ggplot2`, `DT` |

Steps 2-4 run outside R, through the `run-mea` command line entry point
of the kinase-library package. The kinase type, percentile threshold and
permutation count actually used are recorded in the pipeline log written
next to the results (`logs/mea_<contrast>.log`).

### What “enriched at the extremes of the ranked data” means

The ranked list from step 1 is the whole quantified phosphoproteome of
one contrast, ordered by the moderated t-statistic: strongly
up-regulated sites at the top, unchanged sites in the middle, strongly
down-regulated sites at the bottom. Ranking on the t-statistic rather
than on the fold change alone means a site needs both a sizeable change
and a reproducible one to reach an extreme of the list. Pre-ranked GSEA
walks this list from top to bottom with a running sum that increases at
every site belonging to the kinase’s substrate set and decreases at
every site outside it. The largest deviation of that walk from zero is
the **enrichment score (ES)**:

- substrates spread evenly over the list → the walk stays near zero → ES
  ≈ 0, no signal;
- substrates piling up at the **top** of the list → large **positive**
  ES: the motif is enriched at the up-regulated extreme;
- substrates piling up at the **bottom** → large **negative** ES:
  enriched at the down-regulated extreme.

Randomly permuting the ranks calibrates ES into the **normalized
enrichment score (NES)**, which is comparable between kinases and
contrasts, plus a p-value and an FDR. The **leading-edge substrates**
are the sites in front of the running-sum peak: they are the sites that
actually produce the score, and they are listed per kinase in the
`mea_*.csv` files.

So a kinase with NES \> 0 and a small FDR has its predicted substrates
concentrated among the sites that go up in this contrast, which is read
as increased activity of that kinase; NES \< 0 is read as decreased
activity.

Two limitations to keep in mind when interpreting the tables below:

- Substrate sets are **motif predictions**, not experimentally curated
  kinase-substrate pairs. Kinases within a family share near-identical
  motifs, so their sets overlap heavily and their NES values move
  together — MEA points at a motif class rather than at one specific
  kinase. The PTM-SEA report is the complementary view, based on curated
  PhosphoSitePlus kinase-substrate annotations.
- Only sites with a usable sequence window enter the ranking, and a
  site’s window is scored regardless of whether that kinase is expressed
  in the studied system.

### References and further reading

- The Kinase Library web tool: <https://kinase-library.phosphosite.org>
- kinase-library Python package:
  <https://github.com/TheKinaseLibrary/kinase-library>
  ([PyPI](https://pypi.org/project/kinase-library/)); this pipeline
  calls the fork <https://github.com/wolski/kinase-library>, which adds
  the `scan-motifs` and `run-mea` command line entry points
- Ser/Thr kinome atlas: Johnson JL *et al.* (2023) An atlas of substrate
  specificities for the human serine/threonine kinome. *Nature*
  613:759-766. <https://doi.org/10.1038/s41586-022-05575-3>
- Tyrosine kinome atlas: Yaron-Barir TM *et al.* (2024) The intrinsic
  substrate specificity of the human tyrosine kinome. *Nature*
  629:1174-1181. <https://doi.org/10.1038/s41586-024-07407-y>
- GSEA method: Subramanian A *et al.* (2005) Gene set enrichment
  analysis. *PNAS* 102:15545-15550.
  <https://doi.org/10.1073/pnas.0506580102>
- GSEApy implementation: Fang Z, Liu X, Peltz G (2023) GSEApy: a
  comprehensive package for performing gene set enrichment analysis in
  Python. *Bioinformatics* 39:btac757.
  <https://doi.org/10.1093/bioinformatics/btac757>

## Load MEA Results

``` r

if (pipeline_mode) {
  # Pipeline mode: load from kinaselib directory
  mea_files <- list.files(params$kinaselib_dir, pattern = "^mea_.*\\.csv$", full.names = TRUE)

  if (length(mea_files) == 0) {
    stop("No MEA result files found in: ", params$kinaselib_dir)
  }

  # Load all MEA results
  mea_results <- mea_files |>
    set_names(gsub("^mea_|\\.csv$", "", basename(mea_files))) |>
    map_dfr(~ read.csv(.x, stringsAsFactors = FALSE), .id = "contrast")

  output_dir <- params$kinaselib_dir
} else {
  # Vignette mode: use bundled example data
  bundled_zip <- system.file("extdata", "mea_results.zip", package = "prophosqua")

  if (!file.exists(bundled_zip)) {
    stop("Bundled MEA example data not found. Run in pipeline mode with actual data.")
  }

  temp_dir <- tempdir()
  unzip(bundled_zip, exdir = temp_dir)
  mea_files <- list.files(temp_dir, pattern = "^mea_.*\\.csv$", full.names = TRUE)

  mea_results <- mea_files |>
    set_names(gsub("^mea_|\\.csv$", "", basename(mea_files))) |>
    map_dfr(~ read.csv(.x, stringsAsFactors = FALSE), .id = "contrast")

  output_dir <- tempdir()
  message("Using bundled MEA example data from prophosqua package")
}

# Standardize column names (MEA output uses Kinase, we need kinase)
if ("Kinase" %in% names(mea_results)) {
  mea_results <- mea_results |> rename(kinase = Kinase)
}
if ("p.value" %in% names(mea_results)) {
  mea_results <- mea_results |> rename(pvalue = `p.value`)
} else if ("p-value" %in% names(mea_results)) {
  mea_results <- mea_results |> rename(pvalue = `p-value`)
}
# Split the "Subs fraction" column ("85/471") into leading-edge count and set size
subs_col <- intersect(c("Subs.fraction", "Subs fraction"), names(mea_results))
if (length(subs_col) > 0) {
  mea_results <- mea_results |>
    mutate(
      n_leading = as.numeric(sub("/.*", "", .data[[subs_col[1]]])),
      set_size = as.numeric(sub(".*/", "", .data[[subs_col[1]]]))
    )
}

cat("Loaded", nrow(mea_results), "results from", length(mea_files), "contrasts\n")
```

    ## Loaded 622 results from 2 contrasts

``` r

cat("Contrasts:", paste(unique(mea_results$contrast), collapse = ", "), "\n")
```

    ## Contrasts: KO_vs_WT_at_Early, KO_vs_WT

``` r

# Clean and prepare results using shared helper
mea_clean <- prepare_enrichment_data(mea_results, "FDR", 0.1)

# Summary
summary_df <- mea_clean |>
  group_by(contrast) |>
  summarize(
    total_kinases = n(),
    sig_up = sum(FDR < 0.1 & NES > 0, na.rm = TRUE),
    sig_down = sum(FDR < 0.1 & NES < 0, na.rm = TRUE),
    .groups = "drop"
  )

knitr::kable(
  summary_df,
  caption = paste0(
    "Number of kinase motif sets tested per contrast and how many of them are ",
    "significantly enriched (", params$analysis_type, " results, GSEA on sites ranked by ",
    "the site-level moderated t-statistic). total_kinases counts the motif sets tested, ",
    "sig_up the kinases with FDR < 0.1 and NES > 0 (predicted substrates enriched among ",
    "up-regulated sites, i.e. inferred activation) and sig_down those with FDR < 0.1 and ",
    "NES < 0 (inferred inactivation)."
  )
)
```

| contrast          | total_kinases | sig_up | sig_down |
|:------------------|--------------:|-------:|---------:|
| KO_vs_WT          |           311 |     79 |       58 |
| KO_vs_WT_at_Early |           311 |     55 |      118 |

Number of kinase motif sets tested per contrast and how many of them are
significantly enriched (DPA results, GSEA on sites ranked by the
site-level moderated t-statistic). total_kinases counts the motif sets
tested, sig_up the kinases with FDR \< 0.1 and NES \> 0 (predicted
substrates enriched among up-regulated sites, i.e. inferred activation)
and sig_down those with FDR \< 0.1 and NES \< 0 (inferred inactivation).
{.table}

``` r

# Contrast order and per-contrast figure captions, used by the loop below
contrast_ids <- unique(mea_clean$contrast)
dotplot_caps <- paste0(
  "Kinase motif enrichment in contrast ", contrast_ids, " (", params$analysis_type,
  "): normalized enrichment score of the 30 kinases with the smallest FDR. Each point is ",
  "one kinase motif set; the x-axis is the NES from pre-ranked GSEA over sites ranked by ",
  "the site-level moderated t-statistic, so NES > 0 means the predicted substrates ",
  "accumulate among ",
  "up-regulated sites and NES < 0 among down-regulated sites. Point size is -log10(FDR), ",
  "point colour follows the NES (blue negative, red positive), faded points are kinases ",
  "with FDR >= 0.1, and the dashed line marks NES = 0. Kinases are ordered by NES."
)

heatmap_cap <- paste0(
  "Normalized enrichment score per kinase motif set and contrast (", params$analysis_type,
  "). Rows are the 30 kinases with the smallest FDR in any contrast (only kinases reaching ",
  "FDR < 0.15 somewhere are eligible), ordered by NES; columns are the contrasts. Tile ",
  "colour is the NES (blue: substrates enriched among down-regulated sites, red: among ",
  "up-regulated sites, white: no enrichment) and the asterisks give the FDR of that ",
  "kinase in that contrast (* < 0.1, ** < 0.05, *** < 0.01). Empty tiles mean the kinase ",
  "was not tested in that contrast."
)

volcano_cap <- paste0(
  "Enrichment strength versus significance for every tested kinase motif set, one panel ",
  "per contrast (", params$analysis_type, "). Each point is one kinase; the x-axis is the ",
  "NES from pre-ranked GSEA (negative: substrates enriched among down-regulated sites, ",
  "positive: among up-regulated sites) and the y-axis -log10(FDR), with the horizontal ",
  "dashed line at FDR = 0.1 and the vertical dashed line at NES = 0. Colour encodes the ",
  "direction of enrichment, faded points are kinases above the FDR cutoff, and the five ",
  "kinases with the smallest FDR per contrast are labelled. Axes are free per panel."
)
```

## Results by Contrast

``` r

for (ctr in contrast_ids) {
  cat("\n\n## ", ctr, "\n\n")

  ctr_data <- mea_clean |> filter(contrast == ctr)
  n_sig <- sum(ctr_data$FDR < 0.1, na.rm = TRUE)
  cat("**Significant kinases (FDR < 0.1):** ", n_sig, "\n\n")

  # Top kinases dotplot using shared function
  p <- plot_enrichment_dotplot(
    ctr_data,
    item_col = "kinase",
    fdr_col = "FDR",
    title = paste0(params$analysis_type, " - ", ctr),
    subtitle = "Top 30 kinases by FDR"
  )
  print(p)
  cat("\n\n")

  # Significant kinases table
  cat("### Significant Kinases\n\n")
  sig_table <- ctr_data |>
    filter(FDR < 0.1) |>
    select(kinase, NES, pvalue, FDR, n_leading, set_size) |>
    arrange(FDR) |>
    mutate(across(where(is.numeric), ~round(.x, 4)))
  print(htmltools::tagList(
    DT::datatable(sig_table,
                  extensions = 'Buttons',
                  options = list(pageLength = 15, scrollX = TRUE,
                                 dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
                  caption = paste0(
                    "Kinases enriched at FDR < 0.1 in contrast ", ctr, ", sorted by FDR. ",
                    "NES is the normalized enrichment score (> 0 substrates enriched among ",
                    "up-regulated sites, < 0 among down-regulated sites), pvalue the GSEA ",
                    "permutation p-value, FDR its Benjamini-Hochberg adjustment, ",
                    "set_size the number of ranked sites whose sequence window matches the ",
                    "kinase motif, and n_leading how many of them lie in the leading edge ",
                    "and thus drive the score."
                  ))
  ))
  cat("\n\n")
}
```

### KO_vs_WT_at_Early

**Significant kinases (FDR \< 0.1):** 173

![Kinase motif enrichment in contrast KO_vs_WT_at_Early (DPA):
normalized enrichment score of the 30 kinases with the smallest FDR.
Each point is one kinase motif set; the x-axis is the NES from
pre-ranked GSEA over sites ranked by the site-level moderated
t-statistic, so NES \> 0 means the predicted substrates accumulate among
up-regulated sites and NES \< 0 among down-regulated sites. Point size
is -log10(FDR), point colour follows the NES (blue negative, red
positive), faded points are kinases with FDR \>= 0.1, and the dashed
line marks NES = 0. Kinases are ordered by
NES.](Analysis_MEA_files/figure-html/contrastPlots-1.png)

Kinase motif enrichment in contrast KO_vs_WT_at_Early (DPA): normalized
enrichment score of the 30 kinases with the smallest FDR. Each point is
one kinase motif set; the x-axis is the NES from pre-ranked GSEA over
sites ranked by the site-level moderated t-statistic, so NES \> 0 means
the predicted substrates accumulate among up-regulated sites and NES \<
0 among down-regulated sites. Point size is -log10(FDR), point colour
follows the NES (blue negative, red positive), faded points are kinases
with FDR \>= 0.1, and the dashed line marks NES = 0. Kinases are ordered
by NES.

#### Significant Kinases

### KO_vs_WT

**Significant kinases (FDR \< 0.1):** 137

![Kinase motif enrichment in contrast KO_vs_WT (DPA): normalized
enrichment score of the 30 kinases with the smallest FDR. Each point is
one kinase motif set; the x-axis is the NES from pre-ranked GSEA over
sites ranked by the site-level moderated t-statistic, so NES \> 0 means
the predicted substrates accumulate among up-regulated sites and NES \<
0 among down-regulated sites. Point size is -log10(FDR), point colour
follows the NES (blue negative, red positive), faded points are kinases
with FDR \>= 0.1, and the dashed line marks NES = 0. Kinases are ordered
by NES.](Analysis_MEA_files/figure-html/contrastPlots-2.png)

Kinase motif enrichment in contrast KO_vs_WT (DPA): normalized
enrichment score of the 30 kinases with the smallest FDR. Each point is
one kinase motif set; the x-axis is the NES from pre-ranked GSEA over
sites ranked by the site-level moderated t-statistic, so NES \> 0 means
the predicted substrates accumulate among up-regulated sites and NES \<
0 among down-regulated sites. Point size is -log10(FDR), point colour
follows the NES (blue negative, red positive), faded points are kinases
with FDR \>= 0.1, and the dashed line marks NES = 0. Kinases are ordered
by NES.

#### Significant Kinases

## Summary Heatmap

``` r

# Using shared heatmap function
plot_enrichment_heatmap(
  mea_clean,
  item_col = "kinase",
  fdr_col = "FDR",
  n_top = 30,
  title = paste0("MEA Summary Heatmap - ", params$analysis_type)
)
```

![Normalized enrichment score per kinase motif set and contrast (DPA).
Rows are the 30 kinases with the smallest FDR in any contrast (only
kinases reaching FDR \< 0.15 somewhere are eligible), ordered by NES;
columns are the contrasts. Tile colour is the NES (blue: substrates
enriched among down-regulated sites, red: among up-regulated sites,
white: no enrichment) and the asterisks give the FDR of that kinase in
that contrast (\* \< 0.1, \*\* \< 0.05, \*\*\* \< 0.01). Empty tiles
mean the kinase was not tested in that
contrast.](Analysis_MEA_files/figure-html/heatmap-1.png)

Normalized enrichment score per kinase motif set and contrast (DPA).
Rows are the 30 kinases with the smallest FDR in any contrast (only
kinases reaching FDR \< 0.15 somewhere are eligible), ordered by NES;
columns are the contrasts. Tile colour is the NES (blue: substrates
enriched among down-regulated sites, red: among up-regulated sites,
white: no enrichment) and the asterisks give the FDR of that kinase in
that contrast (\* \< 0.1, \*\* \< 0.05, \*\*\* \< 0.01). Empty tiles
mean the kinase was not tested in that contrast.

## Volcano Plot

``` r

# Using shared volcano function
plot_enrichment_volcano(
  mea_clean,
  item_col = "kinase",
  fdr_col = "FDR",
  title = paste0("MEA Volcano Plots - ", params$analysis_type)
)
```

![Enrichment strength versus significance for every tested kinase motif
set, one panel per contrast (DPA). Each point is one kinase; the x-axis
is the NES from pre-ranked GSEA (negative: substrates enriched among
down-regulated sites, positive: among up-regulated sites) and the y-axis
-log10(FDR), with the horizontal dashed line at FDR = 0.1 and the
vertical dashed line at NES = 0. Colour encodes the direction of
enrichment, faded points are kinases above the FDR cutoff, and the five
kinases with the smallest FDR per contrast are labelled. Axes are free
per panel.](Analysis_MEA_files/figure-html/volcano-1.png)

Enrichment strength versus significance for every tested kinase motif
set, one panel per contrast (DPA). Each point is one kinase; the x-axis
is the NES from pre-ranked GSEA (negative: substrates enriched among
down-regulated sites, positive: among up-regulated sites) and the y-axis
-log10(FDR), with the horizontal dashed line at FDR = 0.1 and the
vertical dashed line at NES = 0. Colour encodes the direction of
enrichment, faded points are kinases above the FDR cutoff, and the five
kinases with the smallest FDR per contrast are labelled. Axes are free
per panel.

## Diagnostics

``` r

pval_diag <- mea_clean |>
  group_by(contrast) |>
  summarise(
    `Min p-value` = signif(min(pvalue, na.rm = TRUE), 3),
    `p < 0.05` = sum(pvalue < 0.05, na.rm = TRUE),
    `p < 0.01` = sum(pvalue < 0.01, na.rm = TRUE),
    Total = n(),
    .groups = "drop"
  ) |>
  rename(Contrast = contrast)
knitr::kable(
  pval_diag,
  caption = paste0(
    "Distribution of the unadjusted GSEA permutation p-values per contrast (",
    params$analysis_type, "), before FDR adjustment. Total is the number of kinase motif ",
    "sets tested; a contrast with a real kinase signal shows more sets below the ",
    "thresholds than the roughly 5% and 1% of Total expected by chance. The smallest ",
    "attainable p-value is 1/(number of permutations)."
  )
)
```

| Contrast          | Min p-value | p \< 0.05 | p \< 0.01 | Total |
|:------------------|------------:|----------:|----------:|------:|
| KO_vs_WT          |       0.001 |       145 |       102 |   311 |
| KO_vs_WT_at_Early |       0.001 |       174 |       123 |   311 |

Distribution of the unadjusted GSEA permutation p-values per contrast
(DPA), before FDR adjustment. Total is the number of kinase motif sets
tested; a contrast with a real kinase signal shows more sets below the
thresholds than the roughly 5% and 1% of Total expected by chance. The
smallest attainable p-value is 1/(number of permutations). {.table}

## All Results

``` r

# All kinases across all contrasts
all_results_dt <- mea_clean |>
  select(contrast, kinase, NES, pvalue, FDR, n_leading, set_size) |>
  arrange(contrast, FDR) |>
  mutate(across(where(is.numeric), ~round(.x, 4)))

DT::datatable(all_results_dt,
  filter = "top",
  extensions = 'Buttons',
  options = list(pageLength = 15, scrollX = TRUE,
                 dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
  caption = paste0(
    "Complete MEA result table: every tested kinase motif set in every contrast (",
    params$analysis_type, "), sorted by contrast and FDR, including the ",
    "non-significant ones. NES is the normalized enrichment score (> 0 substrates ",
    "enriched among up-regulated sites, < 0 among down-regulated sites), pvalue the ",
    "permutation p-value, FDR its Benjamini-Hochberg adjustment, set_size the number of ",
    "ranked sites matching the kinase motif and n_leading the leading-edge subset that ",
    "drives the score. Use the column filters to search a kinase."
  ))
```

## Export Results

``` r

# Export to Excel
export_list <- list(
  all_results = mea_clean |>
    select(contrast, kinase, NES, pvalue, FDR, n_leading, set_size) |>
    arrange(contrast, FDR),
  significant = mea_clean |>
    filter(FDR < 0.1) |>
    select(contrast, kinase, NES, pvalue, FDR, n_leading, set_size) |>
    arrange(contrast, FDR),
  summary = summary_df
)

xlsx_file <- file.path(output_dir, paste0("MEA_", params$analysis_type, "_results.xlsx"))
writexl::write_xlsx(export_list, xlsx_file)

rds_file <- file.path(output_dir, paste0("MEA_", params$analysis_type, "_results.rds"))
saveRDS(mea_clean, rds_file)

cat("Exported:\n -", xlsx_file, "\n -", rds_file, "\n")
```

``` r

message("Vignette mode: File export skipped.")
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
    ## [1] prophosqua_0.3.0 writexl_2.0.1    forcats_1.0.1    purrr_1.2.2     
    ## [5] ggplot2_4.0.3    tidyr_1.3.2      DT_0.34.0        dplyr_1.2.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       jsonlite_2.0.0     compiler_4.6.1     tidyselect_1.2.1  
    ##  [5] optparse_1.8.2     jquerylib_0.1.4    systemfonts_1.3.2  scales_1.4.0      
    ##  [9] textshaping_1.0.5  yaml_2.3.12        fastmap_1.2.0      R6_2.6.1          
    ## [13] labeling_0.4.3     patchwork_1.3.2    generics_0.1.4     knitr_1.51        
    ## [17] htmlwidgets_1.6.4  tibble_3.3.1       bookdown_0.47      desc_1.4.3        
    ## [21] RColorBrewer_1.1-3 bslib_0.12.0       pillar_1.11.1      rlang_1.3.0       
    ## [25] cachem_1.1.0       xfun_0.60          S7_0.2.2           fs_2.1.0          
    ## [29] sass_0.4.10        otel_0.2.0         cli_3.6.6          withr_3.0.3       
    ## [33] pkgdown_2.2.1      magrittr_2.0.5     crosstalk_1.2.2    digest_0.6.39     
    ## [37] grid_4.6.1         lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5    
    ## [41] glue_1.8.1         farver_2.1.2       ggseqlogo_0.2.2    ragg_1.5.2        
    ## [45] rmarkdown_2.31     tools_4.6.1        pkgconfig_2.0.3    htmltools_0.5.9
