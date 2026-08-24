# Collect the Motif Enrichment Results of One Analysis

Reads the per-contrast \`mea\_\*.csv\` files the kinase-library tool
wrote and brings them onto the column names the enrichment reports
share. Unlike the two GSEA analyses this does no testing of its own —
the enrichment already happened, one contrast per file.

## Usage

``` r
compute_mea(kinaselib_dir)
```

## Arguments

- kinaselib_dir:

  Directory holding the \`mea\_\*.csv\` files.

## Value

A list with \`mea_clean\`, its per-contrast \`summary_df\`, and
\`has_results\`.

## Examples

``` r
# Needs the mea_*.csv files of a pipeline run.
if (FALSE) { # \dontrun{
res <- compute_mea("PTM_DPA/KinaseLib")
res$summary_df
} # }
```
