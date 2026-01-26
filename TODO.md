# TODO: DRY prophosqua and ptm-pipeline

**Date:** 2026-01-22

## Problem

Significant overlap between prophosqua vignettes and ptm-pipeline templates:

| prophosqua | ptm-pipeline | Overlap |
|------------|--------------|---------|
| `Analysis_PTMSEA.Rmd` | `Analysis_PTMSEA.Rmd` | ~70% |
| `Analysis_Motif_Enrichment.Rmd` | `Analysis_MEA.Rmd` | ~60% |
| `Analysis_KinaseLibrary.Rmd` | `Analysis_KinaseLibrary_GSEA.Rmd` | ~60% |
| `Vis_MEA.Rmd` | (integrated) | visualization |

### Key Differences

- **prophosqua vignettes**: Use package example data, run full analysis, basic clusterProfiler plots
- **ptm-pipeline templates**: Use YAML params, read Excel/CSV, enhanced visualization via `enrichment_viz_utils.R`

## Recommended Solution

### 1. Move visualization functions to prophosqua

Create `R/enrichment_viz.R` with functions from `ptm-pipeline/template/src/enrichment_viz_utils.R`:

```r
# Functions to move:
plot_enrichment_dotplot()
plot_enrichment_volcano()
plot_enrichment_heatmap()
export_gsea_plots_pdf()
extract_gsea_results()
summarize_enrichment_results()
```

### 2. Simplify prophosqua vignettes

Keep vignettes as minimal documentation:
- Load example data
- Run one analysis
- Show one plot using new viz functions

Remove duplicate DPA/DPU sections - just show API usage.

### 3. ptm-pipeline templates become thin wrappers

```r
library(prophosqua)

# 1. Load data (parameterized)
data <- readxl::read_xlsx(params$xlsx_file)

# 2. Run analysis (prophosqua)
results <- run_ptmsea(...)

# 3. Visualize (prophosqua - moved there)
plot_enrichment_heatmap(results, ...)

# 4. Export (template-specific)
writexl::write_xlsx(...)
```

## Implementation Steps

- [ ] Copy `enrichment_viz_utils.R` to `prophosqua/R/enrichment_viz.R`
- [ ] Add roxygen documentation and `@export` tags
- [ ] Update NAMESPACE (devtools::document)
- [ ] Refactor `Analysis_PTMSEA.Rmd` vignette to use new functions
- [ ] Test vignette builds
- [ ] Update ptm-pipeline templates to use `prophosqua::plot_*` functions
- [ ] Remove `enrichment_viz_utils.R` from ptm-pipeline (or keep as fallback)
- [ ] Repeat for other vignette pairs

## Files Reference

### prophosqua vignettes
- `/Users/wolski/projects/prophosqua/vignettes/Analysis_PTMSEA.Rmd`
- `/Users/wolski/projects/prophosqua/vignettes/Analysis_Motif_Enrichment.Rmd`
- `/Users/wolski/projects/prophosqua/vignettes/Analysis_KinaseLibrary.Rmd`
- `/Users/wolski/projects/prophosqua/vignettes/Vis_MEA.Rmd`

### ptm-pipeline templates
- `/Users/wolski/projects/ptm-pipeline/template/src/Analysis_PTMSEA.Rmd`
- `/Users/wolski/projects/ptm-pipeline/template/src/Analysis_MEA.Rmd`
- `/Users/wolski/projects/ptm-pipeline/template/src/Analysis_KinaseLibrary_GSEA.Rmd`
- `/Users/wolski/projects/ptm-pipeline/template/src/enrichment_viz_utils.R`

## Benefits

| Aspect | Before | After |
|--------|--------|-------|
| Viz functions | ptm-pipeline only | Reusable in prophosqua |
| Vignettes | Duplicate analysis code | Simple examples |
| Templates | Copy-paste from vignettes | Call prophosqua directly |
| Maintenance | Fix bugs in 2 places | Fix once |
