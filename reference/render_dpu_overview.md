# Render the Phospho and Protein Integration Overview

Renders the integration overview from the DPU result object. This report
takes R objects rather than file paths as parameters, so it cannot go
through \[render_ptm_report()\]: the objects are built here from the
saved result and handed over directly.

## Usage

``` r
render_dpu_overview(
  input_rds,
  output_dir,
  project_id = "PTM_analysis",
  work_unit_id = "DPU_Integration"
)
```

## Arguments

- input_rds:

  The \`combined_test_diff.rds\` written by the DPU computation.

- output_dir:

  Directory to write \`Result_DPU.html\` to.

- project_id:

  Project identifier shown in the report header.

- work_unit_id:

  Work unit identifier shown in the report header.

## Value

Invisibly, the path of the rendered file.

## Details

The template and its bibliography are copied into a private render
directory first. Rendering in place would write knitr's intermediates
and both copied files into the project root, and would let two
concurrent renders overwrite each other's intermediates.

## Examples

``` r
# Needs the DPU result of a pipeline run.
if (FALSE) { # \dontrun{
render_dpu_overview("PTM_DPU/combined_test_diff.rds", "PTM_DPU")
} # }
```
