# How Much of the Data PTMsigDB Annotates

The two sides are matched on trimmed flanking sequence, which is the
site identifier PTMsigDB keys on. A low overlap bounds how many
signatures can be tested at all, and is the first thing to check when
few pathways pass the size filter.

## Usage

``` r
ptmsigdb_overlap(data, pathways, trim_to)
```

## Arguments

- data:

  Result table with a \`SequenceWindow\` column.

- pathways:

  PTMsigDB signatures.

- trim_to:

  Flanking width to trim to before matching.

## Value

A list with the summary \`stats\` table and the two counts behind it.
