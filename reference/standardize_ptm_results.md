# Standardize a PTM Result Table

DPA, DPU and CF arrive with different column names for the same three
quantities, because each is produced by a different route. This brings
all three onto the one vocabulary every downstream report reads:
\`diff.site\` for the log2 fold change, \`FDR.site\` for its adjusted
p-value and \`statistic.site\` for its test statistic.

## Usage

``` r
standardize_ptm_results(data, analysis_type)
```

## Arguments

- data:

  Data frame of results for one analysis type.

- analysis_type:

  One of \`"dpa"\`, \`"dpu"\` or \`"cf"\`, case-insensitive.

## Value

\`data\` reduced to the standard columns, renamed.

## Details

Columns absent from \`data\` are dropped rather than raising: a DEA run
that did not record \`SequenceWindow\`, say, still produces a usable
table, only without the motif reports. That silence is the reason a
column can be present in \`Result_DPU.xlsx\` and missing from every
report — check the selection below before suspecting the analysis.

## Examples

``` r
dpu <- data.frame(
  protein_Id = "P1", site = "P1~S1", contrast = "a_vs_b",
  gene_name.site = "GENE", diff_diff = 1.5, FDR_I = 0.01, tstatistic_I = 4
)
standardize_ptm_results(dpu, "dpu")
#>   protein_Id  site contrast gene_name diff.site FDR.site statistic.site
#> 1         P1 P1~S1   a_vs_b      GENE       1.5     0.01              4
```
