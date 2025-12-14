# Make FASTA database summary

Generates database statistics including number of sequences and amino
acid frequencies. Extracted from prozor package to reduce dependencies.

## Usage

``` r
make_fasta_summary(resDB, old = FALSE, as_string = TRUE)
```

## Arguments

- resDB:

  A list of sequences from
  [`prozor::readPeptideFasta`](https://rdrr.io/pkg/prozor/man/readPeptideFasta.html)

- old:

  Logical; use old (slower) method for AA frequency calculation

- as_string:

  Logical; return formatted string (TRUE) or list (FALSE)

## Value

If `as_string=TRUE`, an array of strings to pass to
[`cat()`](https://rdrr.io/r/base/cat.html). If `as_string=FALSE`, a list
with nrSequences, lengthSummary, and aafreq.

## Examples

``` r
if (FALSE) { # \dontrun{
fasta <- prozor::readPeptideFasta("path/to/file.fasta")
cat(make_fasta_summary(fasta))
} # }
```
