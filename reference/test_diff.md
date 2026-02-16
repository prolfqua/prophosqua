# Compute MSstats-like test statistics for differential PTM usage

Joins site-level and protein-level results, computes
difference-of-differences, and performs t-tests for differential PTM
usage.

## Usage

``` r
test_diff(
  phos_res,
  tot_res,
  join_column = c("protein_Id", "contrast", "description", "protein_length",
    "nr_tryptic_peptides")
)
```

## Arguments

- phos_res:

  Data frame with phospho site-level results

- tot_res:

  Data frame with protein-level results

- join_column:

  Character vector of columns to join by

## Value

Data frame with combined results including diff_diff test statistics
