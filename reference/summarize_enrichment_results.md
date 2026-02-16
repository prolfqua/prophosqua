# Create summary table for enrichment results

Generates a summary table with counts of significant results at
different FDR thresholds.

## Usage

``` r
summarize_enrichment_results(results, fdr_thresholds = c(0.1, 0.05))
```

## Arguments

- results:

  Named list of clusterProfiler GSEA result objects OR data frame

- fdr_thresholds:

  Vector of FDR thresholds to count (default: c(0.1, 0.05))

## Value

tibble with summary statistics

## Examples

``` r
# Summarize MEA/kinase enrichment data frame
mea_results <- data.frame(
  contrast = c(rep("Treatment_vs_Control", 50), rep("Drug_vs_Vehicle", 50)),
  kinase = paste0("Kinase", 1:100),
  NES = rnorm(100),
  FDR = runif(100, 0, 0.3)
)
summarize_enrichment_results(mea_results)
#> # A tibble: 2 × 4
#>   contrast             total `FDR < 0.1` `FDR < 0.05`
#>   <chr>                <int>       <int>        <int>
#> 1 Drug_vs_Vehicle         50          16           12
#> 2 Treatment_vs_Control    50          19           12
# Returns: contrast, total, FDR < 0.1, FDR < 0.05
```
