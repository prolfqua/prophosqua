# Changelog

## prophosqua 0.3.0

- The PTM pipeline’s analysis code now lives here. Every computation the
  workflow performs is a documented, tested package function, and every
  script and report template it runs is installed with the package under
  `inst/application`:
  [`compute_dpa_dpu()`](https://prolfqua.github.io/prophosqua/reference/compute_dpa_dpu.md),
  [`compute_cf_dea()`](https://prolfqua.github.io/prophosqua/reference/compute_cf_dea.md),
  [`combine_ptm_results()`](https://prolfqua.github.io/prophosqua/reference/combine_ptm_results.md),
  [`prepare_ptmsigdb()`](https://prolfqua.github.io/prophosqua/reference/prepare_ptmsigdb.md),
  [`prep_kinaselib_inputs()`](https://prolfqua.github.io/prophosqua/reference/prep_kinaselib_inputs.md),
  [`compute_ptmsea()`](https://prolfqua.github.io/prophosqua/reference/compute_ptmsea.md),
  [`compute_kinaselib_gsea()`](https://prolfqua.github.io/prophosqua/reference/compute_kinaselib_gsea.md)
  and
  [`compute_mea()`](https://prolfqua.github.io/prophosqua/reference/compute_mea.md),
  reached from the workflow through `CMD_*.R` front ends and the one
  `inst/application/bin/ptm.sh` wrapper, which takes the step as its
  first argument: `ptm.sh dpa_dpu`, `ptm.sh render`, `ptm.sh help`. An
  analysis project no longer carries a copy of any of it, so there is
  nothing left to edit in a working directory that a rerun would
  silently discard.
  [`copy_ptm_shell_script()`](https://prolfqua.github.io/prophosqua/reference/copy_ptm_shell_script.md)
  places the wrapper in a working directory for running the same steps
  by hand.

- Computing and reporting are separate steps. `Analysis_DPA_DPU.Rmd`,
  `Analysis_CorrectFirst_DEA.Rmd`, `Analysis_PTMSEA.Rmd`,
  `Analysis_KinaseLibrary.Rmd` and `Analysis_MEA.Rmd` render from a
  saved result object instead of producing it, so correcting a caption
  or a sentence costs a render rather than a full reanalysis: the
  CorrectFirst report went from about 85 seconds and a refit of every
  site model to about 15 seconds, and the enrichment reports no longer
  repeat their permutation tests.

- The three enrichment analyses always write their result workbook and
  RDS, with zero rows when nothing was enriched, instead of writing
  nothing at all. Their export used to sit inside a chunk that only ran
  when something was found, so an empty enrichment left a stale workbook
  from an earlier run looking current, and no workflow could declare the
  files as expected outputs.

- [`standardize_ptm_results()`](https://prolfqua.github.io/prophosqua/reference/standardize_ptm_results.md)
  replaces the `standardize_results()` script function, with the column
  selection of all three analyses and their order covered by tests. A
  column absent from an analysis is still dropped silently, which is why
  a column can be present in `Result_DPU.xlsx` and missing from every
  report; the selection is now documented and visible in one place.

- [`render_ptm_report()`](https://prolfqua.github.io/prophosqua/reference/render_ptm_report.md)
  renders any of the installed reports, and
  [`render_dpu_overview()`](https://prolfqua.github.io/prophosqua/reference/render_dpu_overview.md)
  the integration overview. Both render from the install path into a
  private intermediates directory, so no project directory collects
  knitr leftovers and two concurrent renders of one template cannot
  overwrite each other’s intermediates.

- `prophosqua` now imports arrow, bookdown, optparse, prolfquapp,
  rmarkdown and writexl, which the moved code needs at run time.

- The nine helpers that locate and normalize prolfquapp DEA outputs are
  now package functions
  ([`get_dea_xlsx()`](https://prolfqua.github.io/prophosqua/reference/get_dea_xlsx.md),
  [`get_dea_file()`](https://prolfqua.github.io/prophosqua/reference/get_dea_file.md),
  [`get_dea_parquet()`](https://prolfqua.github.io/prophosqua/reference/get_dea_parquet.md),
  [`get_dea_yaml()`](https://prolfqua.github.io/prophosqua/reference/get_dea_yaml.md),
  [`get_sample_name_column()`](https://prolfqua.github.io/prophosqua/reference/get_sample_name_column.md),
  [`get_dea_sample_name_column()`](https://prolfqua.github.io/prophosqua/reference/get_dea_sample_name_column.md),
  [`canonicalize_dea_sample_column()`](https://prolfqua.github.io/prophosqua/reference/canonicalize_dea_sample_column.md),
  [`canonicalize_uniprot_ids()`](https://prolfqua.github.io/prophosqua/reference/canonicalize_uniprot_ids.md)
  and
  [`get_dea_ptm_site_info()`](https://prolfqua.github.io/prophosqua/reference/get_dea_ptm_site_info.md)),
  documented with runnable examples and covered by tests. They
  previously lived in a script each analysis project carried its own
  copy of, where a stale copy could shadow a package function of the
  same name and change results silently.

- The N-to-C and sequence-logo reports no longer
  [`source()`](https://rdrr.io/r/base/source.html) a script from the
  calling project’s `src/` directory. They use the package’s own
  functions, so a report cannot pick up a different implementation
  depending on which directory it runs in.

- Sites whose estimate rests on limit-of-detection imputation are marked
  as imputed again in the N-to-C figures. The imputation flag moved from
  `modelName` to `estimate_type`, and the plots still read the old
  column, so every site was drawn as observed; the plotting code now
  reads `estimate_type` and tolerates its absence rather than failing.

- The integrated-PTM manuscript vignette now cites the FragPipe
  computational-platform paper for FragPipe. All three FragPipe
  citations previously resolved to an unrelated paper about an
  insecticidal protein, which reached print in the Methods in Molecular
  Biology chapter proof.

- Manuscript bibliography corrections: accented author names are
  restored throughout (previously every non-ASCII character had been
  stripped, printing “Schfer”, “Villn”, “Ylmaz”); the two Zenodo
  software deposits no longer parse the surname “Wolski” as a given
  name; ggseqlogo and the Delom & Chevet reference gained their missing
  author and title; missing article numbers were added; and both cited
  preprints were updated to their published versions.

- Every figure and table in the N-to-C, sequence logo, PTM-SEA,
  kinase-library GSEA and MEA reports now carries a caption naming the
  plotted quantity, what points or tiles represent, the axes and colour
  encoding, the grouping and the filtering that produced it, so a
  caption identifies its figure on its own. Per-contrast figures name
  their contrast.

- The MEA report now opens with a full method Overview: the software
  chain behind the numbers (kinase-library motif scoring, GSEApy
  pre-ranked GSEA, this R report), an explanation of how enrichment at
  the extremes of the ranked site list is computed and read, the
  interpretation limits of motif-predicted substrate sets, and
  literature and tool links. Every figure and table carries a caption
  naming the metric, axes, encoding and preprocessing.

- MEA result tables and the MEA Excel export replace the ambiguous
  `size` column with `set_size` (ranked sites matching the kinase motif)
  and `n_leading` (the leading-edge subset driving the score). `size`
  previously held the leading-edge count while being labelled as the
  number of substrates.

- Multi-contrast N-to-C figures now carry the protein description
  (identifier, length, number of sites and of not localized sites) once,
  as figure title and subtitle, instead of repeating a truncated title
  above every contrast panel, and the per-panel legends are collected
  into a single legend. Panels stay readable with six or more contrasts.

- Package documentation builds now cover the reproducible analysis
  vignettes; the data-intensive methods and quality-control reports
  remain in their dedicated Snakemake workflow.

- The quality-control report now honors supplied input paths without
  downloading example data, uses the current `LFQData` accessors, and
  normalizes channel totals to the first observed channel.

- The methods report now resolves its non-DEA render against the current
  precomputed April 2026 analysis results, uses the current `LFQData`
  mutation API, and no longer overwrites the packaged example dataset
  when rendered.

- Began tracking user-visible changes in `NEWS.md`. For changes before
  this version, see the git history.
