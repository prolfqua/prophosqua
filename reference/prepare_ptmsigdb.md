# Download, Merge and Filter PTMsigDB

Fetches PTMsigDB for human and mouse, merges the two into one signature
set, keeps only the requested sub-sources, and trims the flanking
sequences to the width the site records of an analysis use. The result
is the signature database PTM-SEA scores against.

## Usage

``` r
prepare_ptmsigdb(
  output_dir = "data/ptmsigdb",
  keep_sources = "KINASE-PSP",
  trim_to = 15
)
```

## Arguments

- output_dir:

  Directory for the filtered database and the download cache.

- keep_sources:

  Sub-sources to keep, as a character vector.

- trim_to:

  Flanking width to trim the sequences to; one of 11, 13 or 15.

## Value

Invisibly, a list with the \`rds\` and \`gmt\` paths written and the
\`pathways\` themselves.

## Details

Human and mouse are merged rather than chosen between because the
signatures are keyed on flanking sequence, not on organism: a site
conserved between the two contributes the same sequence from either, and
a signature curated in only one of them would otherwise be lost for the
other.

Sub-sources present in PTMsigDB v2.0.0 are \`KINASE-PSP\`
(PhosphoSitePlus kinase-substrate, experimental), \`KINASE-iKiP\` (iKiP,
computational), \`PATH-NP\` / \`PATH-WP\` / \`PATH-BI\` (NetPath,
WikiPathways, Broad Institute pathways), \`PERT-PSP\` and
\`PERT-P100-\*\` (perturbations), and \`DISEASE-PSP\`.

## Examples

``` r
# Downloads from PhosphoSitePlus; needs network access.
if (FALSE) { # \dontrun{
prepare_ptmsigdb("data/ptmsigdb", keep_sources = "KINASE-PSP", trim_to = 15)
} # }
```
