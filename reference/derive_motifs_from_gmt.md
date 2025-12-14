# Derive motifs from PTMsigDB GMT file

Analyzes a GMT file (containing lists of aligned substrate sequences for
each kinase) and generates a "consensus" regular expression for each
kinase based on residue frequency.

## Usage

``` r
derive_motifs_from_gmt(
  gmt_file,
  freq_threshold = 0.3,
  min_seqs = 10,
  trim_to = 15
)
```

## Arguments

- gmt_file:

  Path to PTMsigDB GMT file (e.g., "ptm.sig.db.all.flanking...gmt").

- freq_threshold:

  Numeric (0-1). Minimum frequency for a residue (or group) to be
  included in the regex. Default 0.3 (30 percent).

- min_seqs:

  Integer. Minimum number of sequences required to derive a motif for a
  kinase. Default 10.

- trim_to:

  Integer. Width of the analyzed window (must be odd). Default 15.

## Value

data.frame with columns: kinase_group, motif_name, pattern, description
