# Prepare data for PTM-SEA ORA analysis

Prepares enriched and background site lists for over-representation
analysis. When fc_col is provided, significant sites are split into up
and down sets. Only sites that intersect with PTMsigDB are included.

## Usage

``` r
ptmsea_ora_prep(
  data,
  ptmsigdb_sites,
  score_col,
  fc_col = NULL,
  seq_window_col = "SequenceWindow",
  threshold = 0.05,
  trim_to = c("11", "13", "15")
)
```

## Arguments

- data:

  data.frame with flanking sequences, scores, and optionally fold
  changes

- ptmsigdb_sites:

  Character vector. Site IDs from PTMsigDB (stripped, without ;u/;d).
  Use `unique(gsub(";[ud]$", "", unlist(pathways)))`.

- score_col:

  Character. Column for significance score (e.g., "FDR_I").

- fc_col:

  Character or NULL. Column for fold change. If provided, returns
  separate enriched_up and enriched_down. Default NULL.

- seq_window_col:

  Character. Column for flanking sequences. Default "SequenceWindow".

- threshold:

  Numeric. Significance threshold for score_col. Default 0.05.

- trim_to:

  Character. Trim width: "11" (default), "13", or "15".

## Value

List with enriched (or enriched_up/enriched_down if fc_col) and
background.
