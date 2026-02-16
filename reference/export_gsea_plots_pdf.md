# Export GSEA enrichment plots to PDF

Exports clusterProfiler GSEA plots for significant gene sets to a PDF
file. Requires the enrichplot package.

## Usage

``` r
export_gsea_plots_pdf(
  results,
  output_file,
  fdr_threshold = 0.25,
  width = 10,
  height = 6,
  prefix_pattern = NULL
)
```

## Arguments

- results:

  Named list of clusterProfiler GSEA result objects

- output_file:

  Path to output PDF file

- fdr_threshold:

  FDR threshold for including plots (default: 0.25)

- width:

  PDF width in inches (default: 10)

- height:

  PDF height in inches (default: 6)

- prefix_pattern:

  Regex pattern to strip from pathway names for titles (default: NULL)

## Value

Number of plots exported

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires clusterProfiler GSEA results
export_gsea_plots_pdf(gsea_results, "enrichment_plots.pdf")
} # }
```
