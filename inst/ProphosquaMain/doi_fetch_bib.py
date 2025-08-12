#!/usr/bin/env python3
import subprocess
import requests
import sys

output_entries = []

# Read DOIs (one per line) from dois.txt
with open('dois.txt') as fh:
    dois = [line.strip() for line in fh if line.strip()]

for doi in dois:
    bibtex_entry = ""

    # 1) Try Crossref via doi2bib
    try:
        bibtex_entry = (
            subprocess
            .check_output(['doi2bib', doi], stderr=subprocess.DEVNULL)
            .decode('utf-8', errors='ignore')
            .strip()
        )
    except subprocess.CalledProcessError:
        # doi2bib failed, fall back
        pass

    if bibtex_entry:
        output_entries.append(bibtex_entry)
        continue

    # 2) Fallback to DataCite / Zenodo via content negotiation
    try:
        resp = requests.get(
            f'https://doi.org/{doi}',
            headers={'Accept': 'application/x-bibtex'},
            timeout=10
        )
        if resp.status_code == 200 and resp.text.strip():
            output_entries.append(resp.text.strip())
        else:
            print(f"Warning: No BibTeX entry found for DOI {doi}", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Failed to fetch DOI {doi}: {e}", file=sys.stderr)



# in output entries replace &amp; with &
output_entries = [entry.replace('&amp;', '&') for entry in output_entries]

# Write combined output to references.bib
with open('doi_references.bib', 'w') as out_fh:
    out_fh.write('\n\n'.join(output_entries))

print(f"\nDone — wrote {len(output_entries)} entries to references.bib")


