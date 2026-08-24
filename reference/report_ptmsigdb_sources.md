# Report the Sub-Source Breakdown of a PTMsigDB Filter

Prints how many signatures each sub-source contributed and how many
survived the filter, so a run's log says what was kept without the
caller having to inspect the database.

## Usage

``` r
report_ptmsigdb_sources(merged, filtered)
```

## Arguments

- merged:

  All merged pathways.

- filtered:

  The pathways that survived filtering.

## Value

Invisibly, \`NULL\`.
