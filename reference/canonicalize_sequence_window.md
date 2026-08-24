# Give the Flanking Sequence Column its Canonical Name

The flanking sequence arrives as \`SequenceWindow\` from a prolfquapp
DEA and as \`PTM_FlankingRegion\` from some search engines' output. Both
mean the same thing.

## Usage

``` r
canonicalize_sequence_window(data)
```

## Arguments

- data:

  Data frame that should carry a flanking sequence column.

## Value

\`data\` with the column named \`SequenceWindow\`.
