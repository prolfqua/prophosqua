# Filter contaminant proteins from data

Removes rows matching FGCZ contaminants, contam\_ prefix, or rev\_
prefix.

## Usage

``` r
filter_contaminants(data, protein_id_col = "protein_Id")
```

## Arguments

- data:

  Data frame with protein data

- protein_id_col:

  Name of the protein ID column (default: "protein_Id")

## Value

Filtered data frame
