# Prepare Kinase-Library Inputs from the Combined PTM Results

Writes the two file kinds the kinase-library tool consumes: one list of
unique sequence windows for motif scanning, and one ranked file per
contrast for motif enrichment analysis.

## Usage

``` r
prep_kinaselib_inputs(xlsx_file, output_dir, analysis_type, sheet, stat_column)
```

## Arguments

- xlsx_file:

  The combined \`PTM_results.xlsx\`.

- output_dir:

  Where to write, normally the analysis \`KinaseLib\` directory.

- analysis_type:

  \`DPA\`, \`DPU\` or \`CF\`; names the output files.

- sheet:

  Sheet of \`xlsx_file\` to read.

- stat_column:

  Column to rank the sites on.

## Value

Invisibly, a character vector of the files written.

## Details

The workbook, sheet and ranking statistic are the same ones PTM-SEA and
the KinaseLib GSEA read, so all three enrichment analyses rank the sites
identically and their results can be compared.

## Examples

``` r
# Needs a combined workbook from a pipeline run.
if (FALSE) { # \dontrun{
prep_kinaselib_inputs(
  "PTM_results.xlsx", "PTM_DPA/KinaseLib",
  analysis_type = "DPA", sheet = "DPA", stat_column = "statistic.site"
)
} # }
```
