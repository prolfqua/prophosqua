grep -E -o '10\.[0-9]{4,}/[A-Za-z0-9._\/-]+' Methods_in_MB_PTM2.md | sort -u > dois.txt
grep -E -o '10\.[0-9]{4,}/[A-Za-z0-9._\/-]+' Methods_in_MB_Supplement_Wetlab.md | sort -u >> dois.txt


grep -E -o '10\.[0-9]{4,}/[A-Za-z0-9._\/-]+' ../../vignettes/Supplementary_Material_v2.Rmd | sort -u > dois.txt

wc dois.txt

python doi_fetch_bib_v2.py dois.txt


