# Run PTM-SEA without direction

Runs GSEA on PTMsigDB pathways, ignoring direction suffixes (;u/;d).

## Usage

``` r
run_ptmsea(
  ranks_list,
  pathways,
  min_size = 3,
  max_size = 500,
  n_perm = 1000,
  pvalueCutoff = 0.1
)
```

## Arguments

- ranks_list:

  Named list of ranked vectors from ptmsea_data_prep().

- pathways:

  List of pathways from fgsea::gmtPathways().

- min_size:

  Integer. Minimum pathway size. Default 3.

- max_size:

  Integer. Maximum pathway size. Default 500.

- n_perm:

  Integer. Number of permutations. Default 1000.

- pvalueCutoff:

  Numeric. Adjusted p-value cutoff. Default 0.1.

## Value

Named list of clusterProfiler gseaResult objects (one per contrast).
