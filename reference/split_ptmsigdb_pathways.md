# Run PTM-SEA analysis

Applies fgsea to each contrast's ranked site list against PTMsigDB
pathways.

## Usage

``` r
split_ptmsigdb_pathways(pathways)
```

## Arguments

- pathways:

  Named list of pathway definitions (e.g., from fgsea::gmtPathways()).

- ranks_list:

  Named list of named numeric vectors from ptmsea_data_prep().

- min_size:

  Integer. Minimum pathway size. Default 3.

- max_size:

  Integer. Maximum pathway size. Default 500.

- n_perm:

  Integer. Number of permutations. Default 1000.

## Value

Named list of fgsea result data.frames, one per contrast.

## Examples

``` r
# Mock ranked lists
# Must have > 10 overlap with pathways to run
seqs <- paste0("SEQ", 1:15, "-p")
rank_vec <- setNames(rnorm(15), seqs)
ranks_list <- list(contrast1 = sort(rank_vec, decreasing = TRUE))

# Mock pathways
# Create pathways that overlap with the sequences
pathways <- list(
  PathwayA = seqs[1:10],
  PathwayB = seqs[6:15]
)

# Run PTM-SEA (adjusting min_size for small mock data)
results <- run_ptmsea(ranks_list, pathways, min_size = 1, pvalueCutoff = 1.0)
#> Running PTM-SEA for 'contrast1' (15 sites, 15 overlap)
```
