# prophosqua API cleanup

The cleanup pass is done: prophosqua's R code is half its size, and the exported API is what the pipeline, the MiMB chapter and ptmbenchmark call. What remains are fixes that belong upstream, listed below for a decision.

## Measured against the start of the cleanup

- R/: 8,412 → 4,127 lines (−51%), 5,523 → 2,882 code lines (−48%)
- R files: 41 → 29
- Exports: 85 → 35; CMD scripts: 14 → 6
- tests/testthat: 2,945 → 2,068 lines

## Verified

- `make check`: 0 errors, 0 warnings, 0 notes
- o43037 `ptm-pipeline run`: 29/29 jobs, exit 0
- DPA, DPU, CF and the eight statistics sheets of `PTM_results.xlsx`: identical to the run before the cleanup
- enrichment: same terms and set sizes; NES r ≥ 0.9999, the permutations being unseeded
- the run caught one regression, a CBOR manifest check that refused every assembled file; fixed, and a test now reads an assembled file back

## What the pass found and removed

Code with no caller, measured by a call graph from the six CMD scripts the Snakefile runs, the two Quarto reports, the MiMB manuscript, `_QCReport.qmd` and ptmbenchmark:

- the xlsx/parquet/yaml DEA-folder readers (`get_dea_xlsx()` and six more)
- the file-based enrichment entry points and their eight CMD scripts
- `compute_*_h5ad()`, `read_ptm_anndata_pair()`: a second path next to MuData
- the embedded string_gsea builder (`gsea_result_json.R`), a copy of protsea's
- report helpers for the removed R Markdown vignettes, `render_dpu_overview()`, `create_top_index.Rmd`

Ad hoc fixes and compat shims:

- `site_column()`: the `protein_Id_site` spelling of old DEAs
- `canonicalize_sequence_window()`: a `PTM_FlankingRegion` no reader emits
- `.ptm_delivery_missing()`: "" → NA to match the former Excel delivery
- `derive_contrasts()` and the annotation file of `compute_cf_dea()`: the DEAs record their contrasts
- the `G_` column CF required after it moved to the DEAs' formula
- dual MEA column spellings; a default for a setting the pipeline always sets
- placeholder result fields (`data_info = data.frame(value = 1)`) written only to satisfy validators

Indirections:

- `build(Type, ...)`, which was `Type$new(self, ...)`
- `CF` and `DPA_DPU` stage classes, never persisted alone, with a deep compare of inputs to compose them
- six copy-paste enrichment classes, now one base class `PTM_enrichment`
- per-stage `write_h5mu()`/loaders for stages the pipeline never writes
- string_gsea validators (≈180 lines) re-checking protsea's own format, three times per document

## Decided

- contaminants: kept by default, as the DEAs keep them; the PTM parameter `remove_contaminants` (and the argument of `compute_dpa_dpu()`, `compute_cf_dea()`) drops the features the DEAs flag `CON`, once, when the pair is built. The hard-coded pattern of DPA/DPU is gone; it dropped o43037's 29 contaminant proteins but none of their 42 sites, whose ids are accessions.
- protein ids: `canonicalize_uniprot_ids()` stays for now. The total DEA reports `sp|P12345|NAME`, the site DEA `P12345`; the fix belongs upstream, in the preprocess functions of prolfquapp and prolfquappPTMreaders, to be taken up later.
- enrichment storage: MuData holds quant and statistics only. PTM-SEA, Kinase GSEA and MEA are protsea `.json.gz` documents beside `PTM_results.h5mu`, which keeps their names and checksums; the kinase-library preparations stay CBOR. The embedded form and the bundled `.h5mu` example are gone.

## For a decision: fixes at the source

- `DEAResultReader` could expose the formula as a string and the contrasts as a named vector; prophosqua unwraps the metadata tables today.
- string_gsea validation belongs in protsea; `decode_gsea_json()` checks nothing.
- `prolfquapp::read_h5mu()` prints `stack imbalance` warnings: a PROTECT mismatch in the C code of the HDF5 layer under anndataR, not in prophosqua.

## Kept, used outside the pipeline

- MiMB chapter: `compute_dpa_dpu()`, `compute_cf_dea()`, the multi-contrast N-to-C plots, `plot_seqlogo_with_diff()`, `copy_phospho_integration()`; its xlsx helpers went when it moved to the new DEAs
- QC report: `make_fasta_summary()`
- ptmbenchmark: `compute_dpa_dpu()`, `compute_cf_dea()`; its call now passes no annotation file

## ptm-pipeline

- `template/helpers.py`: the five DEA-path helpers are removed; templates resolve through `prophosqua::report_file()`
- the three skills are parked as `SKILL.md.bak` until they are brought up to date
- new configs carry `remove_contaminants: false`
