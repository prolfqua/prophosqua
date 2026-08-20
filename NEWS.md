# prophosqua 0.3.0

- The nine helpers that locate and normalize prolfquapp DEA outputs are now package
  functions (`get_dea_xlsx()`, `get_dea_file()`, `get_dea_parquet()`, `get_dea_yaml()`,
  `get_sample_name_column()`, `get_dea_sample_name_column()`,
  `canonicalize_dea_sample_column()`, `canonicalize_uniprot_ids()` and
  `get_dea_ptm_site_info()`), documented with runnable examples and covered by tests. They
  previously lived in a script each analysis project carried its own copy of, where a stale
  copy could shadow a package function of the same name and change results silently.
- The N-to-C and sequence-logo reports no longer `source()` a script from the calling
  project's `src/` directory. They use the package's own functions, so a report cannot pick
  up a different implementation depending on which directory it runs in.
- Sites whose estimate rests on limit-of-detection imputation are marked as imputed again in
  the N-to-C figures. The imputation flag moved from `modelName` to `estimate_type`, and the
  plots still read the old column, so every site was drawn as observed; the plotting code now
  reads `estimate_type` and tolerates its absence rather than failing.

- The integrated-PTM manuscript vignette now cites the FragPipe computational-platform paper for FragPipe. All three FragPipe citations previously resolved to an unrelated paper about an insecticidal protein, which reached print in the Methods in Molecular Biology chapter proof.
- Manuscript bibliography corrections: accented author names are restored throughout (previously every non-ASCII character had been stripped, printing "Schfer", "Villn", "Ylmaz"); the two Zenodo software deposits no longer parse the surname "Wolski" as a given name; ggseqlogo and the Delom & Chevet reference gained their missing author and title; missing article numbers were added; and both cited preprints were updated to their published versions.
- Every figure and table in the N-to-C, sequence logo, PTM-SEA, kinase-library GSEA and MEA reports now carries a caption naming the plotted quantity, what points or tiles represent, the axes and colour encoding, the grouping and the filtering that produced it, so a caption identifies its figure on its own. Per-contrast figures name their contrast.
- The MEA report now opens with a full method Overview: the software chain behind the numbers (kinase-library motif scoring, GSEApy pre-ranked GSEA, this R report), an explanation of how enrichment at the extremes of the ranked site list is computed and read, the interpretation limits of motif-predicted substrate sets, and literature and tool links. Every figure and table carries a caption naming the metric, axes, encoding and preprocessing.
- MEA result tables and the MEA Excel export replace the ambiguous `size` column with `set_size` (ranked sites matching the kinase motif) and `n_leading` (the leading-edge subset driving the score). `size` previously held the leading-edge count while being labelled as the number of substrates.
- Multi-contrast N-to-C figures now carry the protein description (identifier, length, number of sites and of not localized sites) once, as figure title and subtitle, instead of repeating a truncated title above every contrast panel, and the per-panel legends are collected into a single legend. Panels stay readable with six or more contrasts.
- Package documentation builds now cover the reproducible analysis vignettes; the data-intensive methods and quality-control reports remain in their dedicated Snakemake workflow.
- The quality-control report now honors supplied input paths without downloading example data, uses the current `LFQData` accessors, and normalizes channel totals to the first observed channel.
- The methods report now resolves its non-DEA render against the current precomputed April 2026 analysis results, uses the current `LFQData` mutation API, and no longer overwrites the packaged example dataset when rendered.
- Began tracking user-visible changes in `NEWS.md`. For changes before this version, see the git history.
