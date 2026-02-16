# Validate sequence window alignment

Filters phosphosite data to keep only rows where the central residue of
the SequenceWindow matches the reported modified amino acid (modAA).

## Usage

``` r
validate_sequence_window(
  data,
  seq_col = "SequenceWindow",
  mod_col = "modAA",
  center_pos = 8L
)
```

## Arguments

- data:

  Data frame with PTM results containing SequenceWindow and modAA
  columns

- seq_col:

  Character. Name of sequence window column. Default "SequenceWindow"

- mod_col:

  Character. Name of modified amino acid column. Default "modAA"

- center_pos:

  Integer. Position of the central residue (1-indexed). Default 8

## Value

Filtered data frame with only valid sequence windows

## Details

Sequence windows should be centered on the modified residue. This
function validates that the character at the center position matches the
reported modification site. Both are converted to uppercase for
comparison.

## Examples

``` r
data <- data.frame(
  SequenceWindow = c("AAASAAAA", "BBBSBBB", "CCCACCC"),
  modAA = c("S", "S", "S")
)
validate_sequence_window(data)
#> [1] SequenceWindow modAA         
#> <0 rows> (or 0-length row.names)
# Only first two rows kept (center matches S)
```
