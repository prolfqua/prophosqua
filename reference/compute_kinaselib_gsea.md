# Compute Kinase-Library GSEA

Ranks the sites of each contrast by their site-level statistic and runs
pre-ranked GSEA against the kinase-substrate sets a motif scan assigned.

## Usage

``` r
compute_kinaselib_gsea(
  xlsx_file,
  sheet,
  term2gene_file,
  min_size = 15,
  max_size = 5000,
  n_perm = 1000,
  max_kinase_sets = NULL
)
```

## Arguments

- xlsx_file:

  Combined \`PTM_results.xlsx\`.

- sheet:

  Sheet to read.

- term2gene_file:

  Kinase-to-site assignment written by the motif scan.

- min_size, max_size:

  Set-size filter for testable kinase sets.

- n_perm:

  Permutations of the enrichment test.

- max_kinase_sets:

  Subsample the kinase sets to at most this many, for a quick standalone
  run. \`NULL\`, the default, tests all of them.

## Value

A list with \`gsea_results\`, \`all_results\`, the summary tables the
report shows, and \`has_results\`.

## Examples

``` r
# Needs a combined workbook and a term2gene assignment.
if (FALSE) { # \dontrun{
res <- compute_kinaselib_gsea(
  "PTM_results.xlsx", "DPA", "PTM_DPA/KinaseLib/term2gene.csv"
)
res$gsea_info
} # }
```
