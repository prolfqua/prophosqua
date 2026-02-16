# Prepare data for N-to-C plotting

Prepares PTM data for N-to-C visualization by ensuring required columns
exist and standardizing column names.

## Usage

``` r
prepare_ntoc_data(data, analysis_type = "dpa")
```

## Arguments

- data:

  Data frame with PTM results

- analysis_type:

  Character. Analysis type: "dpa", "dpu", or "cf"

## Value

Data frame prepared for N-to-C plotting functions

## Details

The function handles different analysis types which may have different
column naming conventions:

- DPA: Uses diff.site, FDR.site columns

- DPU: Uses diff_diff, FDR_I columns

- CF: Uses diff, FDR columns
