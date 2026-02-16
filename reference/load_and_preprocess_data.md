# Load and preprocess DEA data from Excel

Reads an Excel file and validates that required columns are present.

## Usage

``` r
load_and_preprocess_data(
  file_path,
  required_cols,
  sheet_name = "diff_exp_analysis"
)
```

## Arguments

- file_path:

  Path to the Excel file

- required_cols:

  Character vector of required column names

- sheet_name:

  Name of the sheet to read (default: "diff_exp_analysis")

## Value

Data frame with the loaded data
