# Extract results from clusterProfiler GSEA objects

Converts a list of clusterProfiler GSEA result objects into a tidy
tibble.

## Usage

``` r
extract_gsea_results(results)
```

## Arguments

- results:

  Named list of clusterProfiler GSEA result objects

## Value

tibble with columns: ID, NES, pvalue, p.adjust, setSize, contrast

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires clusterProfiler GSEA results
df <- extract_gsea_results(gsea_results)
} # }
```
