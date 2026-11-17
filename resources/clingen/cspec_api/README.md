# CSpec API snapshot (latest specifications)

Every ClinGen criteria specification (GN001 ... ) as returned by the CSpec API:

    https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/GNxxx

- `GNxxx.json` — the API response, saved unchanged. Specifications of every status
  are kept (Released, In Prep, Submitted, Deleted). Ids the API does not know
  (e.g. GN233) are absent.
- `manifest.tsv` — one row per file: GN id, genes, status, version, last update,
  VCEP, number of criteria, comparison with the June 2026 snapshot
  (`../cspec_json`, the source of the engine tables), title, SHA-256.
  `vs_2026-06-05`: `new` / `version a -> b` / `content` (same version, different
  criteria text or applicability) / `same`. Versions are compared without trailing
  zeros (the API now writes "2.0" for "2.0.0").
- `FETCHED_AT` — UTC time of the download and the number of files.

This folder is for review only; the engine does not read it. The tables in
`resources/shared/` are still built from `../cspec_json`.

Refresh (overwrites the files here, so `git diff` shows what changed):

    python scripts/fetch_cspec_api.py

The API JSON carries the criteria text per strength only. The general
"VCEP Specifications" paragraph and the comments shown on the CSpec web page
(`https://cspec.genome.network/cspec/ui/svi/doc/GNxxx`) are not included; they are
in the page PDFs in `../cspec_pdf/` (`python scripts/fetch_cspec_api.py --pdf`).
