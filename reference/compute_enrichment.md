# Compute the Enrichment Analyses of the PTM Pipeline

The three enrichment analyses — PTM-SEA against PTMsigDB, GSEA against
kinase-library motifs, and motif enrichment analysis of the pre-ranked
site lists — used to be computed inside the reports that display them.
Permutation tests are the slowest part of the pipeline, so a fix to a
caption cost a full recomputation, and their result workbooks were
written from inside a conditional chunk and could therefore not be
declared as rule outputs.

## Details

These functions do the computing. Each returns everything its report
needs, under the names the report uses, and each \*\*always\*\* produces
a result — empty tables and \`has_results = FALSE\` when nothing was
enriched, rather than no result at all. That is what makes the workbooks
declarable: an enrichment that found nothing is a fact to record, not a
missing file.
