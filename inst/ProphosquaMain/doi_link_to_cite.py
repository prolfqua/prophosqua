#!/usr/bin/env python3
import re, sys

if len(sys.argv) != 3:
    print("Usage: doi_link_to_cite.py input.md output.md")
    sys.exit(1)

inp, outp = sys.argv[1], sys.argv[2]
text = open(inp).read()

# match [DOI](https://doi.org/DOI) and replace with [@DOI]
pattern = re.compile(r'\[([^]]+)\]\(https?://doi\.org/\1\)')
new_text = pattern.sub(r'[@\1]', text)

# Prepend YAML header
header = """---
bibliography: doi_references_key.bib
---

"""

new_text = header + new_text


with open(outp, 'w') as f:
    f.write(new_text)

print(f"Wrote citations to {outp}")
