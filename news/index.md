# Changelog

## prophosqua 0.3.0

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
