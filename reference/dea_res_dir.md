# Get DEA results directory path

Constructs the path to a DEA results directory based on project
structure.

## Usage

``` r
dea_res_dir(project_dir, WU, date)
```

## Arguments

- project_dir:

  Path to the project directory

- WU:

  Work unit identifier (e.g., "phospho", "prot")

- date:

  Date string used in folder naming (e.g., "20240101")

## Value

Character string with full path to the results directory
