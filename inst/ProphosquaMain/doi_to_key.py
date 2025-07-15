#!/usr/bin/env python3
import bibtexparser

# Load your existing .bib
with open('doi_references.bib') as fh:
    bib = bibtexparser.load(fh)

# For each entry, lowercase the DOI and use it as the ID
for e in bib.entries:
    if 'doi' in e:
        doi_lower = e['doi'].lower()
        e['doi'] = doi_lower      # normalize the doi field
        e['ID']  = doi_lower      # set the entry key to the lowercase DOI

# Write out a new .bib
with open('doi_references_key.bib', 'w') as fh:
    bibtexparser.dump(bib, fh)
