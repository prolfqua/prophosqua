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

## Examples

``` r
# Mock ranked lists
seqs <- paste0("SEQ", 1:15, "-p")
rank_vec <- setNames(rnorm(15), seqs)
ranks_list <- list(contrast1 = sort(rank_vec, decreasing = TRUE))

# Mock pathways
pathways <- list(
  PathwayA = seqs[1:10],
  PathwayB = seqs[6:15]
)

# Run PTM-SEA (adjusting min_size for small mock data)
results <- run_ptmsea(ranks_list, pathways, min_size = 1, pvalueCutoff = 1.0)
#> Running PTM-SEA for 'contrast1' (15 sites, 15 overlap)
#> 
```
