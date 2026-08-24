# Compute PTM-SEA against PTMsigDB

Ranks the sites of each contrast by \`stat_column\`, matches them to
PTMsigDB signatures on their trimmed flanking sequence, and runs
pre-ranked GSEA.

## Usage

``` r
compute_ptmsea(
  xlsx_file,
  sheet,
  stat_column,
  ptmsigdb_file,
  trim_to = 15,
  min_size = 10,
  max_size = 500,
  n_perm = 1000
)
```

## Arguments

- xlsx_file:

  Combined \`PTM_results.xlsx\`.

- sheet:

  Sheet to read.

- stat_column:

  Per-site statistic the sites are ranked on.

- ptmsigdb_file:

  Filtered PTMsigDB, as \`.rds\` or \`.gmt\`.

- trim_to:

  Flanking width the database uses.

- min_size, max_size:

  Set-size filter for testable signatures.

- n_perm:

  Permutations of the enrichment test.

## Value

A list with \`results\` (one clusterProfiler object per contrast),
\`all_clean\` (all results as one table), \`pathways\`, the summary
tables the report shows, and \`has_results\`.

## Details

\`trim_to\` must match the width the signature database was trimmed to,
or nothing matches: the site identifier PTMsigDB keys on \*is\* the
flanking sequence.

## Examples

``` r
# Needs a combined workbook and a filtered PTMsigDB.
if (FALSE) { # \dontrun{
res <- compute_ptmsea(
  "PTM_results.xlsx", "DPA", "statistic.site",
  "data/ptmsigdb/ptmsigdb_filtered_KINASE-PSP_15mer.rds"
)
res$results_info
} # }
```
