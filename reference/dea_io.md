# Locate and Normalize prolfquapp DEA Outputs

Helpers that resolve the files a prolfquapp DEA run writes, and
normalize the two identifiers that differ between runs: the sample
column, whose name the run declares in its YAML, and the protein
identifier, which may arrive as a bare accession or in FASTA
\`sp\|P12345\|PROT_HUMAN\` form.

## Details

A DEA output directory holds a \`Results_WU\_\<workunit\>\` subdirectory
with the results workbook, the normalized abundances, and the analysis
configuration.
