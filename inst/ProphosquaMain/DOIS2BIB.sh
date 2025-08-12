grep -E -o '10\.[0-9]{4,}/[A-Za-z0-9._\/-]+' Methods_in_MB_PTM2.md | sort -u > dois.txt

wc dois.txt

python doi_fetch_bib.py
python doi_to_key.py


python doi_link_to_cite.py Methods_in_MB_PTM2.md Methods_in_MB_PTM2_cite.md
pandoc Methods_in_MB_PTM2_cite.md \
  --citeproc \
  -o Methods_in_MB_PTM2_cite.docx

