# Read kinase motifs from CSV

Helper function to load a custom motif database from a file.

## Usage

``` r
read_kinase_motifs_csv(file)
```

## Arguments

- file:

  Path to CSV file. Expected columns: "kinase_group", "motif_name",
  "pattern", "description".

## Value

data.frame suitable for \`scan_motifs\` or
\`get_kinase_motifs(extra_motifs=...)\`
