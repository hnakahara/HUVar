# CSpec web-page PDFs (Released specifications)

The PDF that each CSpec page offers ("PDF" button on
`https://cspec.genome.network/cspec/ui/svi/doc/GNxxx`), for every specification with
status Released in `../cspec_api/manifest.tsv` (128 files, fetched 2026-11-17).

Unlike the API JSON (`../cspec_api/`), the PDF contains the whole page: the general
"VCEP Specifications" paragraph of each criterion, the comments, and the
per-strength text (e.g. the RASopathy rule "PM1 and PM5 may be used in conjunction
at moderate levels", the IL2RG PM2 FAF cutoff 0.000124).

- `GNxxx.pdf` — the PDF as rendered by the CSpec PDF service
  (`/cspec/ui/svi/pdf?url=<page>?layout=print`).
- `text/GNxxx.txt` — plain text from the PDF (poppler `pdftotext`), for grep and
  `git diff`.
- `manifest.tsv` — GN id, genes, status, the specification version the PDF was made
  from, download date, size, SHA-256, title.

Not read by the engine. Refresh (only PDFs whose specification version changed
are downloaded again; each PDF takes about 10 s to render):

    python scripts/fetch_cspec_api.py --pdf          # API JSON, then the PDFs
    python scripts/fetch_cspec_api.py --pdf-only     # PDFs only, from ../cspec_api/manifest.tsv
    python scripts/fetch_cspec_api.py --pdf-only --pdf-ids GN043,GN129
    python scripts/fetch_cspec_api.py --pdf-only --pdf-all   # every status

Size: about 214 MB for the PDFs, 4 MB for the text copies. The repository's
`.gitignore` excludes `*.pdf`, so only `text/`, `manifest.tsv` and this README are
tracked; the PDFs stay local (re-download them with the commands above).
