# Write a Pathway List as GMT

GMT is one line per set: the set name, a description field, then its
members, all tab separated. PTMsigDB carries no per-set description, so
the field is written as \`NA\` — the placeholder the format expects
rather than an empty column, which some readers would collapse.

## Usage

``` r
write_gmt(pathways, file)
```

## Arguments

- pathways:

  Named list of character vectors.

- file:

  Path to write.

## Value

Invisibly, \`file\`.

## Examples

``` r
gmt <- tempfile(fileext = ".gmt")
write_gmt(list(KINASE_A = c("SITE1", "SITE2")), gmt)
readLines(gmt)
#> [1] "KINASE_A\tNA\tSITE1\tSITE2"
```
