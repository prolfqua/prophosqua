# Trim flanking sequence to specified width

Trims a flanking sequence centered on the phosphosite (position 8 in
15-mer). Removes equal amino acids from N and C termini.

## Usage

``` r
trim_flanking_seq(seq, trim_to = 15L)
```

## Arguments

- seq:

  Character. Flanking sequence (15-mer expected).

- trim_to:

  Integer. Target width: 15 (no trim), 13, or 11. Default 15.

## Value

Character. Trimmed sequence.
