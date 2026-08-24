# Drop Sequence Windows a Motif Scan Cannot Use

A window is usable only if it is a full, uninterrupted stretch of
residues. Windows padded with underscores because the site sits near a
protein terminus, and windows shorter than seven residues, carry too
little context for a motif match and are dropped rather than scanned.

## Usage

``` r
filter_sequence_windows(data)
```

## Arguments

- data:

  Data frame with a \`SequenceWindow\` column.

## Value

\`data\` reduced to usable windows, upper-cased.
