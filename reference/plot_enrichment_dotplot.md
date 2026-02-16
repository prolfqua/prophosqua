# Create enrichment dotplot for top items by FDR

Create enrichment dotplot for top items by FDR

## Usage

``` r
plot_enrichment_dotplot(
  data,
  item_col = "kinase",
  fdr_col = "FDR",
  n_top = 30,
  title = NULL,
  subtitle = "Top 30 by FDR"
)
```

## Arguments

- data:

  Data frame with columns: item (kinase/pathway), NES, p.adjust/FDR,
  contrast

- item_col:

  Name of item column (default: "kinase")

- fdr_col:

  Name of FDR column (default: "FDR")

- n_top:

  Number of top items to show (default: 30)

- title:

  Plot title

- subtitle:

  Plot subtitle

## Value

ggplot object

## Examples

``` r
# Kinase enrichment results (e.g., from KinaseLibrary)
kinase_results <- data.frame(
  kinase = c("PKACA", "AKT1", "MAPK1", "CDK1", "CSNK2A1", "GSK3B"),
  NES = c(2.1, 1.8, -1.5, 1.2, -0.8, 0.3),
  FDR = c(0.001, 0.01, 0.02, 0.08, 0.15, 0.45)
)
plot_enrichment_dotplot(kinase_results, n_top = 6, title = "Kinase Enrichment")


# Pathway enrichment with custom column names
pathway_results <- data.frame(
  pathway = c("Cell cycle", "Apoptosis", "MAPK signaling"),
  NES = c(1.9, -1.4, 1.1),
  p.adjust = c(0.005, 0.03, 0.09)
)
plot_enrichment_dotplot(pathway_results, item_col = "pathway",
                        fdr_col = "p.adjust", n_top = 3)
```
