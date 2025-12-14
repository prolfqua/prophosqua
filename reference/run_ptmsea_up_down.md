# Run directional PTM-SEA analysis

Splits PTMsigDB pathways by direction (up/down) and runs fgsea twice per
contrast (ascending and descending rank order).

## Usage

``` r
run_ptmsea_up_down(
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

  Named list of named numeric vectors from ptmsea_data_prep().

- pathways:

  Named list of pathway definitions (e.g., from fgsea::gmtPathways()).

- min_size:

  Integer. Minimum pathway size. Default 3.

- max_size:

  Integer. Maximum pathway size. Default 500.

- n_perm:

  Integer. Number of permutations. Default 1000.

- pvalueCutoff:

  Numeric. Adjusted p-value cutoff. Default 0.1.

## Value

Named list of clusterProfiler gseaResult objects, named
\<contrast\>\_\<ascending/descending\>.
