# Run Motif-based GSEA

Run Motif-based GSEA

## Usage

``` r
run_motif_gsea(
  ranks_list,
  motif_term2gene,
  min_size = 3,
  max_size = 5000,
  n_perm = 1000,
  pvalueCutoff = 0.1
)
```

## Arguments

- ranks_list:

  Named list of named numeric vectors (stats per contrast).

- motif_term2gene:

  data.frame from \`scan_motifs\`.

- min_size:

  Minimum set size.

- max_size:

  Maximum set size.

- n_perm:

  Number of permutations.

- pvalueCutoff:

  P-value cutoff.

## Value

List of GSEA results.

## Examples

``` r
# Mock data: 100 sequences, first 10 weighted to include PKA motifs
seqs <- paste0("SEQ", 1:100, "-p")
ranks <- rnorm(100)
# Make top ranked sequences PKA-like (MotifA)
ranks[1:5] <- c(3.0, 2.5, 2.0, 1.9, 1.8)
names(ranks) <- seqs
ranks <- sort(ranks, decreasing = TRUE)
rankList <- list(contrast1 = ranks)

# Define MotifA as PKA-like sequences (matches SEQ1-SEQ5)
motif_term2gene <- data.frame(
  term = rep("MotifA_PKA", 5),
  gene = seqs[1:5],
  stringsAsFactors = FALSE
)

# Run GSEA
res <- run_motif_gsea(rankList, motif_term2gene, min_size = 3, pvalueCutoff = 1.0)
#> 
```
