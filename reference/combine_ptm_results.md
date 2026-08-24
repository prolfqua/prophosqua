# Combine the Three PTM Analyses into One Workbook

Reads the DPA, DPU and CF result workbooks, standardizes their columns,
and writes them beside the normalized abundances as a single multi-sheet
workbook and the same content as an RDS. Every downstream report of the
pipeline reads this workbook rather than the three separate ones.

## Usage

``` r
combine_ptm_results(
  dpa_xlsx,
  dpu_xlsx,
  cf_xlsx,
  protein_parquet,
  site_parquet,
  output_xlsx,
  output_rds
)
```

## Arguments

- dpa_xlsx, dpu_xlsx, cf_xlsx:

  The three result workbooks.

- protein_parquet, site_parquet:

  Normalized abundances of the two DEA runs.

- output_xlsx, output_rds:

  Where to write the combined result.

## Value

Invisibly, the list that was written.

## Details

Three abundance sheets are produced: protein abundances, site
abundances, and site abundances corrected by their protein — the last
matching what the CorrectFirst analysis models, so a reader can check a
CF result against the values behind it.

## Examples

``` r
# Needs the outputs of a full pipeline run.
if (FALSE) { # \dontrun{
combine_ptm_results(
  "PTM_DPA/Result_DPA.xlsx",
  "PTM_DPU/Result_DPU.xlsx",
  "PTM_CF_DPU/CorrectFirst_PTM_usage_results.xlsx",
  "DEA_total/Results_WU_x/lfqdata_normalized.parquet",
  "DEA_phospho/Results_WU_y/lfqdata_normalized.parquet",
  "PTM_results.xlsx", "PTM_results.rds"
)
} # }
```
