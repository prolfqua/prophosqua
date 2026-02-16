# Split PTMsigDB pathways by direction

Splits pathway gene sets into up- and down-regulated subsets based on
direction suffixes (;u and ;d) used in PTMsigDB.

## Usage

``` r
split_ptmsigdb_pathways(pathways)
```

## Arguments

- pathways:

  Named list of pathway gene vectors with ;u/;d suffixes

## Value

Named list with \_up and \_down pathway subsets
