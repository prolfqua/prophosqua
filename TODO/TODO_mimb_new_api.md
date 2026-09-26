# MiMB chapter on the new API

The chapter ([_MiMBIntegratedPTM.Rmd](../manuscript/_MiMBIntegratedPTM.Rmd)) runs on DEAs from prolfquapp 2.10.5 and reads them with `compute_dpa_dpu()` and `compute_cf_dea()`. Left open: the deposited results on Zenodo, which come from the old DEAs.

## Data

- inputs: Zenodo 10.5281/zenodo.15879865, `PTM_experiment_FP_22_Maculins_and_QC.zip`, downloaded by the chapter when `PTM_example/` is absent
- DEAs: made by the chapter itself with `dodea = TRUE`; `dodea = FALSE` reads the folders of 2026-09-26 in [inst/PTM_example_analysis_v2](../inst/PTM_example_analysis_v2) (git-ignored)
- the old DEAs (2026-04-01, 2026-07-29) predate prolfquapp 2.10.5 and have no `AnnData.h5ad`

## What changed in the chapter

- DEA paragraph: the default model is now `lm_impute`, so a feature that cannot be fitted is refitted at the limit of detection; the old text said missing values are not imputed
- DEA outputs: `AnnData.h5ad` added to the list
- `definePaths`: DEA folders from `prolfquapp::zipdir_name(..., date = )` instead of `dea_xlsx_path()`
- `b4loaddata` to `b8computeptmusage`: `compute_dpa_dpu(..., remove_contaminants = TRUE)` instead of the xlsx readers, `filter_contaminants()`, the left join, the match rates and `test_diff()`
- correct-first: `compute_cf_dea(..., remove_contaminants = TRUE)` instead of parquet, yaml and inline modelling; the text describes its model, `lm_impute` with the median protein abundance added back; volcano from `prolfqua::ContrastsPlotter` on `cf$results`
- unchanged: N-to-C plots, sequence logos, the integration report, the Excel export, tables

## Why the old code could not stay

On the fresh DEAs the site DEA names proteins by accession and the total DEA by FASTA id, so the chapter's own join matched 0 of 86,932 sites (84,000 on the April DEAs), DPU was empty, and `n_to_c_usage_multicontrast()` stopped on zero proteins. `compute_dpa_dpu()` maps both to the accession: 96.7% of sites match.

## Open

- the results deposited on Zenodo (10.5281/zenodo.15830988: integration report and xlsx) come from the old DEAs and old code; they need a new version to match the chapter
- `n_to_c_usage_multicontrast()` and `n_to_c_expression_multicontrast()` fail on an input with no significant protein instead of returning no plots
- the protein-id mapping belongs in the preprocess functions of prolfquapp and prolfquappPTMreaders (see [TODO_api_cleanup.md](TODO_api_cleanup.md))
