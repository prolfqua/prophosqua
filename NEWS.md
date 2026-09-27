# prophosqua 0.3.0

- Reports, the workbook and the enrichments use only the rows whose site estimate is observed; every analysis is still computed from all `lm_impute` models, and MuData keeps every row. `get_tables()` of `PTM_statistics` and `PTM_results` returns the observed DPA, DPU and CF rows by default and every row with `estimates = "all"`, so PTM-SEA, the kinase inputs and the Kinase GSEA rank observed sites only. `get_estimate_counts()` counts the rows of each analysis and contrast by estimate type before the filter; the workbook carries it as the `estimate_counts` sheet, and each Summary tab of the statistics report shows it. The report filters with the exported `observed_site_estimates()`. The CorrectFirst the reports, the workbook (`CF`, `abundances_site_cf`, `CF_intensities`) and the enrichments show is the variant corrected with the protein DEA's imputed values, `correct_first_protein_imputed`, returned by the new `PTM_statistics$get_cf_reported()`; CF and the site-imputed variant stay in MuData. `test_diff()` returns only the matched site and protein pairs: the site-only and protein-only rows it appended, which carried no usage difference, and their `measured_In` column are gone, so DPU has no more rows than there are matched sites. The DPU availability table leaves the report, and the unused match-rate section leaves `_Overview_PhosphoAndIntegration_site.Rmd`. The statistics report's CorrectFirst DPU section opens with a Summary tab, as DPA and DPU do, which replaces the Model summary tab. The three Summary tabs share one layout, the estimate counts beside the observed rows by modified residue; the contrasts and the site-protein match rates, which concern all three analyses, move to the Overview. Every report figure that grows with the contrasts is sized by them, so reports with many contrasts no longer squeeze their panels: the sequence and difference logos by their rows of panels; the statistics volcanoes and site-versus-protein plots, and the enrichment volcanoes, by ggplot2's panel grid, a single-contrast figure's size per panel; the enrichment dot plots and heatmaps by their contrast columns and term rows. The displayed width grows with the figure, so panels and text keep their size on the page.

- Code that corrected other packages' output is gone. Decoys are dropped once, when the DEA pair is built, by the decoy pattern each DEA recorded and `prolfqua::is_decoy()`, instead of a hard-coded `^rev_` in CorrectFirst and the workbook. The formula and the contrasts come from `DEAResultReader`'s `formula` and `contrast_definitions`. The kinase assignments name the windows they were given, so `term2gene` is no longer upper-cased. The sequence windows are cut from the FASTA by the PTM readers, centred on the modified residue and padded with X at a protein terminus, so `validate_sequence_window()` and the `toupper()`/`trimws()` of windows are removed, and terminus windows are recognized by the X padding the readers use rather than by underscores, which no reader writes; on o43037 and the MiMB example, 848 and 448 site windows are padded and now stay out of the kinase inputs and the sequence logos. The MiMB chapter follows.

- Enrichment results never enter MuData, and protsea serializes them. Each analysis keeps its completed PTM-SEA, Kinase GSEA and MEA as protsea documents, `result_ptm_sea.json.gz`, `result_kinase_gsea.json.gz` and `result_mea.json.gz`, written by `protsea::write_gsea_result_json()` or, for the MEA, by the kinase-library tool; the two kinase-library preparations stay gzipped CBOR. `PTM_results.h5mu` holds the statistics and only the names and checksums of these files, and reading it refuses a missing or changed file. The documents no longer carry a `prophosqua` member: the tables derived from each result (`all_clean`, `all_results`, `gsea_info`, `mea_clean`, `summary_df`) are rebuilt from its decoded gseaResult objects, and the MEA tables come from the kinase-library document rather than from its separate result table. prophosqua no longer packs gseaResult objects itself, `MotifEnrichment` is gone (`MEA` follows `KinaseAssignments`), `compute_ptm_enrichment_cbor()` and `assemble_ptm_cbor()` are `compute_ptm_enrichment()` and `assemble_ptm_results()`, and the command `enrich_cbor` is `enrich`. The bundled example is a directory, `inst/extdata/ptm_results_example/`, laid out as a pipeline run lays it out; the embedded form that put the documents into MuData is removed.

- The PTM-SEA and kinase-library GSEA stages stop with the name of a missing set-size or permutation parameter (`gsea$min_size`, `gsea$max_size`, `gsea$n_perm`, `kinaselib$gsea_max_size`). The kinase-library maximum no longer defaults to 5000, and `clusterProfiler::GSEA` reads a missing size as no limit, so an unset key used to change the result silently.

- `compute_cf_dea()` on a site DEA fitted without `lm_impute` reports that the DEA must be rerun with the `lm_impute` model, instead of failing inside dplyr.

- `clusterProfiler`, `fgsea` and `fgczQuartoTemplate`, which the enrichment stages and `render_ptm_report()` call, move from Suggests to Imports; `limpa` and `htmltools` are no longer suggested.

- The example DEA pairs use the sample name as the file name and no longer carry a `raw_file` column. This needs prolfqua >= 1.8.0, whose `setup_analysis()` accepts one column as both `sample_name` and `file_name`.

- The MiMB chapter runs on DEAs from prolfquapp 2.10.5, which it reads with `compute_dpa_dpu()` and `compute_cf_dea()`. Its fixed DEA date for `dodea = FALSE` is 20260926. The xlsx route it used can no longer pair sites with proteins, because the site DEA now names proteins by accession and the total DEA by FASTA id; `compute_dpa_dpu()` maps both to the accession. `dea_xlsx_path()`, `dea_res_dir()`, `load_and_preprocess_data()` and `filter_contaminants()`, which only the chapter called, are removed, and prophosqua no longer depends on `readxl`.

- Contaminants are kept, as the DEAs keep them. The `remove_contaminants` parameter of the PTM inputs, and the argument of the same name of `compute_dpa_dpu()` and `compute_cf_dea()`, drops the sites and proteins the DEAs flag as contaminants (`CON`) from both experiments, so DPA, DPU, CorrectFirst and the workbook see the same features. DPA/DPU used to drop contaminant proteins by a hard-coded pattern that could not match the site DEA's accessions, so their contaminant sites stayed in.

- The exported API is what the PTM pipeline, the MiMB chapter and ptmbenchmark call: 34 exports instead of 85, and six command scripts instead of fourteen. Removed with no remaining caller: the DEA-folder readers of the old xlsx, parquet and yaml outputs (`get_dea_xlsx()`, `get_dea_parquet()`, `get_dea_yaml()`, `get_sample_name_column()`, `get_dea_sample_name_column()`, `canonicalize_dea_sample_column()`, `get_dea_ptm_site_info()`); the file-based enrichment entry points and their commands (`compute_ptmsea()`, `compute_kinaselib_gsea()`, `compute_mea()`, `prep_kinaselib_inputs()`, `prepare_ptmsigdb()`, `write_gmt()`, `mea_gsea_result_data()`, `read_mea_ranks()`, `write_gsea_result_json()`, `CMD_PTMSEA.R`, `CMD_KINASELIB_GSEA.R`, `CMD_MEA.R`, `CMD_PREP_KINASELIB.R`, `CMD_PREP_PTMSIGDB.R`); `compute_dpa_dpu_h5ad()`, `compute_cf_dea_h5ad()`, `read_ptm_anndata_pair()`, `CMD_DPA_DPU.R` and `CMD_CF_DEA.R`; `render_dpu_overview()` with `CMD_DPU_OVERVIEW.R`; `ptm_dpa_dpu_report_data()`, `ptm_enrichment_report_data()`, `plot_enrichment_dotplot()`, `export_gsea_plots_pdf()`, `summarize_enrichment_results()`, `summarize_significant_sites()`, `prepare_ntoc_data()`, `n_to_c_expression()`, `n_to_c_usage()`, `explode_multisites()`, `split_ptmsigdb_pathways()`, `run_ptmsea_up_down()`, `ptmsea_ora_prep()`, `prepare_n_to_c_data()`; the `combined_test_diff_example` dataset and `ptmsigdb_kinase.rds.zip`; `create_top_index.Rmd` and `integration_MSStats_multisite_singlesite.R`. The rank, PTMsigDB and identifier helpers the pipeline uses internally are no longer exported. `report_file()` is exported, since the pipeline resolves its templates with it.

- `PTM_statistics$new(inputs)` computes DPA, DPU and CorrectFirst; the `CF` and `DPA_DPU` stage classes, which were never persisted on their own, are gone, and so is `build()`: a stage is constructed with `Type$new(source, ...)`. The six enrichment stages share one base class, `PTM_enrichment`, and each keeps only the fields a report, the workbook or the kinase-library tool reads; the summary tables no reader used (`data_info`, `ptmsigdb_summary`, `overlap_stats`, `prep_info`, `results_info`, `kl_info`, `assignment_stats`, `kinase_stats`, `ranks_info`, `analysis_inputs`, `has_results`, and the copies of `ranks`, `pathways` and `term2gene`) are no longer computed. `import_ptm_h5mu()` takes `parameters` and `ptmsigdb` and imports the reference data itself, so `CMD_IMPORT_H5MU.R` no longer reaches into the namespace. `compute_cf_dea()` fits the contrasts the DEAs recorded and takes no annotation file. `render_ptm_report()` renders the two Quarto reports only. The site key is `site` everywhere; the `protein_Id_site` spelling of older DEAs is no longer accepted.

- The enrichment documents are built by protsea, which owns the string_gsea format: prophosqua's copy of the document builder and its JSON writer are removed, and PTM-SEA ranks sequence windows without PTMsigDB's `-p` suffix, which is dropped from the database sites at matching, so the documents need no identifier rewriting. prophosqua checks a stored document's wrapper, checksum and identity only; the checks of protsea's own structure (running scores, hit positions, term fields) are gone.

- prophosqua reads a DEA artifact only through `prolfquapp::DEAResultReader` and no longer decodes `uns/prolfquapp`, the DEA `varm` blocks or the DEA layers itself. Every table it takes from a DEA is keyed by the artifact's feature keys (`protein_Id`, and `site` for the site DEA) and tables are merged by joins on those keys and the sample key, never by row position: the DEA results are joined to the feature annotation, the site-imputed CF variant takes the imputed value where the site's imputation route is `fitted` by a join, and a stored result block keeps its `site` column so reading it back joins the site annotation on it. `compute_dpa_dpu()` reads each DEA folder's `AnnData.h5ad`, as `compute_cf_dea()` does, instead of the DE workbook. The CF variant matrices are samples x site ids. The example and test DEA pairs are made by prolfquapp's own DEA on synthetic abundances instead of a hand-written h5ad. Needs the prolfquapp 2.10.5 reader.

- CorrectFirst values reach the workbook one way only: the `abundances_site_cf` sheet is the model's own `wide_data`, with the sample's protein median added back, instead of a second `site − protein` computation that had drifted from it. The old xlsx/parquet path that held a third copy, `combine_ptm_results()` with its reader `read_normalized_abundances()` and the `CMD_COMBINE_RESULTS.R` command, is removed; the pipeline builds the workbook from MuData. Every CF and `correct_first_protein_imputed` result row counts the imputed protein values its fit used (`n_protein_imputed`). `correct_first_site_protein_imputed` is no longer fitted and is kept as its corrected-abundance layer only: its filled site values are the site model's own predictions, so a fit on them reused the site's data and understated the variance (on o43037, 2,646 sites at FDR < 0.05 against 292 for the protein-imputed variant on the same sites).

- Remove the R Markdown analysis vignettes (`Analysis_*.Rmd`) and their bibliography: the PTM pipeline renders only the two Quarto reports, `ptm_statistics.qmd` and `ptm_enrichment.qmd`, so knitr is no longer a vignette builder. The example helpers only those vignettes called, `compute_kinaselib_gsea_example()` and `compute_mea_example()`, go with them, and so do their bundled inputs `mea_results.zip` and `term2gene.csv.zip`.

- CorrectFirst fits the model formula the DEAs recorded (`uns$prolfquapp$formula`) instead of a fixed `~ G_`, and computes two variants beside it from the prolfquapp `imputedData` layers: `correct_first_protein_imputed` (observed site minus imputed protein) and `correct_first_site_protein_imputed` (site filled from its own lm fit only, sites refitted at the LOD keep their gaps, minus imputed protein). CF and both variants live in one `enriched_CF` modality, which replaces `cf` and has its own sites, those of `enriched` whose protein the proteome quantified, in the enriched order: `X` holds CF, each variant is a layer named by its result key, and `varm` holds the results of CF and the protein-imputed variant. The DPU enrichment documents move to `enriched`, beside the DPU results. Reports, enrichment and delivery still read the current CF. The per-sample median of the observed protein abundances is added back instead of a constant 20, and the CF refit at the LOD uses the lower quartile of corrected values seen once in a group. `compute_cf_dea()` reads each DEA folder's `AnnData.h5ad`, and both DEAs need prolfquapp >= 2.10.5 with the `lm_impute` model, which writes `imputedData`.

- The current CF keeps the site-contrast pairs its LOD refit estimates, flagged `estimate_type == "lod_imputed"`, instead of dropping them as `MissingInOneCondition`; `model_counts` counts pairs by `estimate_type`, and `n_before` equals the number of result rows. Sites without a corrected value in any sample are dropped before the fit.

- Every stage handoff is a gzipped CBOR envelope again, the completed PTM-SEA, Kinase GSEA and MEA stages included. Their envelope carries the string_gsea document as JSON text, so the pipeline exchanges one format and the results archive ships `result_*.cbor.gz`.

- Keep the enrichment payload out of `PTM_results.h5mu`. The final MuData records where its stage artifacts are, in `uns/prophosqua/enrichment_cbor`, and reads them back from beside itself; it no longer copies the enrichment documents and kinase preparations into the file. On one three-analysis order that is 1834 MB of a 2345 MB file, every byte of which already existed as a stage artifact on disk. Moving the MuData without its analysis folders now fails naming the missing file.
- Store a completed PTM-SEA, Kinase GSEA or MEA stage as gzipped JSON through protsea, which owns the string_gsea format, so `result_*.json.gz` opens with any gzip reader and parses with `protsea::read_gsea_json()`. The kinase preparations have no such format and stay gzipped CBOR. Both compress about 2.5x. A stored document names the statistics it was computed from, the binding the CBOR envelope already carried.

- Follow the fgczQuartoTemplate figure conventions in both Quarto reports: figures inherit the template's compact size and 40% width and rely on the lightbox for the full-resolution view, multi-panel chunks (rank distributions, enrichment map and tree) are laid out in two columns, and the renderer no longer forces `fig_retina = 1`, so the zoomed image carries twice the on-page resolution. Network and enrichment-map term labels wrap at underscores and slashes and the graph legends move below the panel, so no label is cut at the figure edge.

- Export one `PTM_results.xlsx` workbook from the final MuData, with the statistics, CorrectFirst intensity/annotation, and all nine enrichment result tables. Stop writing separate per-analysis Excel and RDS files.

- Show the selected DPA, DPU, or CorrectFirst ranking as the sole input in the enrichment report's visual overview, with a larger opening sentence naming the analysis.

- Keep leading-edge site lists in MuData JSON but omit them from the enrichment report tables, making the standalone HTML smaller and easier to scan.

- Retain every tested PTM-SEA and Kinase GSEA term in the final MuData JSON, including terms above the report FDR threshold. The enrichment report now has a graphical overview and separate explanations of all three methods.

- Render the enrichment QMD for one selected analysis at a time, with no running-score plots; the JSON retains those curves. Kinase GSEA now uses a separate maximum set size (default 5000) so large substrate sets are tested. When enrichplot cannot construct a similarity tree, the report still shows the enrichment map.
- Render the two installed PTM Quarto reports directly from MuData, with the statistics report accepting `PTM_statistics.h5mu` before enrichment finishes.
- Align the report-template dependency with protsea and prolfquapp so CI can resolve the full report stack.

- Enrichment and kinase preparation steps now exchange compact CBOR artifacts instead of full MuData snapshots. Final assembly validates all nine JSON documents and adds them to one `PTM_results.h5mu`.

- N-to-C sticks now carry a head at the site log2 fold change, an open head and dashed stick mark imputed estimates, significant sites are labelled with residue and position by `ggrepel`, and residues use the colour-blind safe Okabe-Ito palette. The legends are titled Residue, Estimate and Protein.

- Load with current DOSE releases, where the registered `gseaResult` class is
  no longer exported directly from the namespace.
- Keep the MiMB manuscript and its source assets outside `vignettes/` so package vignette builds process only recognized vignette sources.

- Add the first FGCZ Quarto multitab report, `ptm_statistics.qmd`, covering DPA, DPU, and CorrectFirst DPU from one final MuData input. Each analysis now includes compact sequence and difference logos generated during rendering, with relative residue positions and site counts. Sequence-logo sections place the plot first and their counts table second in a third-level tabset. The FDR and absolute log2 fold-change thresholds are report parameters. Package vignette builds use a deterministic two-contrast, 72-site final H5MU with mixed sequence backgrounds and visible position-zero differences, plus synchronized `fgczQuartoTemplate` assets.

- Add the companion FGCZ Quarto multitab report, `ptm_enrichment.qmd`. DPA, DPU, and CorrectFirst each expose PTM-SEA, Kinase GSEA, and MEA summaries, plots, complete searchable result tables, and running-score plots decoded exclusively from the nine versioned JSON documents in final MuData. Native MEA JSON now preserves GSEApy's source running scores and hit positions in the same schema used by the R GSEA methods.

- Extend the MuData enrichment report with ranked-site ridge plots, gene-set networks, enrichment maps, and term-similarity trees reconstructed from the same portable JSON documents. Plot-specific dimensions keep dense network views readable and wide volcano panels compact; volcano labels now point inward and remain inside the figure.

- Store every PTMSEA, KinaseGSEA, and MEA result as validated JSON in MuData. The final typed result exposes complete document getters, verifies schema and checksums when reading, and reconstructs temporary clusterProfiler objects through protsea without persisting `gseaResult` objects.

- Declare the `limpa` dependency used by the existing CorrectFirst vignette so an isolated `R CMD check` can rebuild it.

- Preserve missing gene annotations as `NA` in terminal DPA/DPU RDS exports, matching the legacy workbook reader.

- Preserve the legacy model, estimate-type, contrast, and sample-design column positions in terminal delivery exports.

- Describe the MuData report inputs and terminal delivery files consistently in the DPA/DPU and CorrectFirst reports.

- Complete R6 stages now carry paired DEA, independent DPA/DPU and CorrectFirst, preparation, enrichment, and final results through MuData. Reports read final MuData; Excel/RDS delivery exports run last. Preserve existing joins, contrasts, and numerical computations. Replace the single-site PTM H5AD writer with the MuData commands.

- Read current prolfquapp DEA artifacts and compose complete paired-input, CorrectFirst, and DPA/DPU stages in MuData.

- A paired site and total-proteome AnnData analysis can now be written as one
  `PTM_results.h5ad` with `compute_ptm_results_h5ad()`. The new artifact keeps
  the site measurements and upstream DEA results intact, adds feature-aligned
  DPA, moderated and unmoderated DPU, and CorrectFirst matrices, and records a
  versioned result mapping plus hashes and schema versions for both inputs.
- DPA, DPU, and CorrectFirst can now consume an explicit pair of site-level
  and total-proteome `AnnData.h5ad` files written by prolfquapp. The reader
  validates schema, experiment roles, sample identity, and shared design before
  reconstructing the existing statistical inputs; the legacy DEA-directory
  entry points remain available during numerical comparison.
- Fixed the MEA ranked lists in the enrichment JSON: the `.rnk` files carry a
  header row, which was read as data, so every MEA gene pool began with a bogus
  `SEQUENCEWINDOW` entry, all ranking statistics were strings instead of numbers,
  and every rank was off by one. Regenerate `MEA_*_results.json` to pick up the
  fix.
- The three enrichment CMD scripts (PTM-SEA, kinase-library GSEA, MEA) now also
  write their results as `*.json` in the string_gsea `GSEAResult` structure
  (per contrast: the ranked sequence windows as a shared pool, per term the
  full mapped membership plus the leading edge), via the new exported
  `gsea_result_data()`, `mea_gsea_result_data()`, `read_mea_ranks()` and
  `write_gsea_result_json()`. Pool ids are the canonical upper-case sequence
  windows of the differential table; the form submitted to the enrichment
  (PTMsigDB `-p` suffix, kinase-library lower-case phospho residue) is kept in
  `input_label`. Downstream tools such as ptm3d consume these files.
- `test_diff()` now explicitly selects a coherent moderated or unmoderated
  standard-error/degrees-of-freedom pair. `compute_dpa_dpu()` keeps moderated
  DPU as the published result and also returns an unmoderated comparison table,
  together with the count of rows that cannot be tested because their raw
  degrees of freedom are invalid. Older DEA outputs without the new raw
  contrast columns are rejected and must be regenerated with current
  `prolfqua`.
- The analysis reports are the package's vignettes again. All seven -- DPA/DPU,
  CorrectFirst, PTM-SEA, KinaseLibrary, MEA, N-to-C and seqlogo -- live in
  `vignettes/`, the one place an analysis lives, and the vignette machinery
  installs them into `doc/`, from where a pipeline run renders them with its own
  parameters. There is no second copy under `inst/application` to drift from the
  documented one; only the index page and the overview include, which are not
  analyses, still ship there.
- Every report renders without a pipeline run. Each declares parameter defaults
  that fall back to the example data the package bundles, so `make
  build-vignettes` builds all seven, and the DPA/DPU and CorrectFirst reports --
  which previously defaulted to hardcoded paths inside a project directory --
  now compute their example from a synthetic pair of DEA output directories via
  the same functions the pipeline calls.
- The PTM pipeline's analysis code now lives here. Every computation the workflow performs
  is a documented, tested package function, and every script and report template it runs is
  installed with the package under `inst/application`: `compute_dpa_dpu()`,
  `compute_cf_dea()`, `combine_ptm_results()`, `prepare_ptmsigdb()`,
  `prep_kinaselib_inputs()`, `compute_ptmsea()`, `compute_kinaselib_gsea()` and
  `compute_mea()`, reached from the workflow through `CMD_*.R` front ends and the one
  `inst/application/bin/ptm.sh` wrapper, which takes the step as its first argument:
  `ptm.sh dpa_dpu`, `ptm.sh render`, `ptm.sh help`. An analysis project no longer carries a
  copy of any of it, so there is nothing left to edit in a working directory that a rerun
  would silently discard. `copy_ptm_shell_script()` places the wrapper in a working
  directory for running the same steps by hand.
- Computing and reporting are separate steps. `Analysis_DPA_DPU.Rmd`,
  `Analysis_CorrectFirst_DEA.Rmd`, `Analysis_PTMSEA.Rmd`, `Analysis_KinaseLibrary.Rmd` and
  `Analysis_MEA.Rmd` render from a saved result object instead of producing it, so
  correcting a caption or a sentence costs a render rather than a full reanalysis: the
  CorrectFirst report went from about 85 seconds and a refit of every site model to about
  15 seconds, and the enrichment reports no longer repeat their permutation tests.
- The three enrichment analyses always write their result workbook and RDS, with zero rows
  when nothing was enriched, instead of writing nothing at all. Their export used to sit
  inside a chunk that only ran when something was found, so an empty enrichment left a
  stale workbook from an earlier run looking current, and no workflow could declare the
  files as expected outputs.
- `standardize_ptm_results()` replaces the `standardize_results()` script function, with the
  column selection of all three analyses and their order covered by tests. A column absent
  from an analysis is still dropped silently, which is why a column can be present in
  `Result_DPU.xlsx` and missing from every report; the selection is now documented and
  visible in one place.
- `render_ptm_report()` renders any of the installed reports, and `render_dpu_overview()`
  the integration overview. Both render from the install path into a private
  intermediates directory, so no project directory collects knitr leftovers and two
  concurrent renders of one template cannot overwrite each other's intermediates.
- `prophosqua` now imports arrow, bookdown, optparse, prolfquapp, rmarkdown and writexl,
  which the moved code needs at run time.

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
- Keep enrichment round-trip tests compatible with the gseaResult move from DOSE to enrichit in Bioconductor 3.23.
