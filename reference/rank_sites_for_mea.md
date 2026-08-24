# Rank the Sites of One Contrast for Motif Enrichment

Motif enrichment walks a single ranked list, so a sequence window may
appear only once. Where several sites share a window the most extreme
statistic is kept — the one that would move the enrichment score most —
rather than an average, which would let two opposing sites cancel into
an uninformative middle.

## Usage

``` r
rank_sites_for_mea(data, stat_column, contrast_name)
```

## Arguments

- data:

  Data frame of one analysis, already filtered.

- stat_column:

  Column to rank on.

- contrast_name:

  Contrast to select.

## Value

A two-column table, \`SequenceWindow\` and \`statistic.site\`, ordered
from the largest statistic down.
