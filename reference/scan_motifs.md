# Scan sequences for kinase motifs

Scans a vector of sequences against a database of regex motifs.
Generates a TERM2GENE data.frame for GSEA/Enrichment analysis.

## Usage

``` r
scan_motifs(sequences, motif_db = get_kinase_motifs(), center_pos = 8L)
```

## Arguments

- sequences:

  Character vector of sequences (e.g., 15-mers).

- motif_db:

  data.frame of motifs (must have 'motif_name' and 'pattern' columns).
  Defaults to \`get_kinase_motifs()\`.

- center_pos:

  Integer. Index of the central phosphosite in the sequence. If NULL,
  matches anywhere. If provided, the regex match MUST include/overlap
  this position. Default is 8 (for 15-mers).

## Value

data.frame with columns 'term' (motif name) and 'gene' (sequence ID).
Sequence ID is the original sequence with "-p" appended (std PTMSEA
format).

## Examples

``` r
# Create some mock sequences (15-mers)
# PKA_AKT_consensus pattern: "R.R..[ST]....."
# To center S at 8, the pattern R.R..S must match indices 3-13.
# We'll create sequences with R at 3 and 5, and S at 8.

seqs <- c(
  "GGRRRGSEVVVAAAA", # Matches PKA (R at 3,5, S at 8)
  "LARRRASVAQLTTAA", # Matches PKA
  "AARARAASVAAAAAA", # Matches PKA
  "RRRRAASVAAAAAAA", # Matches PKA and others
  "GGGGGGSGGGGGGGG" # Non-match control
)

motif_db <- get_kinase_motifs()
res <- scan_motifs(seqs, motif_db)
head(res)
#>                term              gene
#> 1 PKA_AKT_consensus AARARAASVAAAAAA-p
#> 2             CAMK2 AARARAASVAAAAAA-p
```
