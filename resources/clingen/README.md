# ClinGen VCEP allele-frequency thresholds (BA1 / BS1)

`disease_prevalence.tsv` supplies per-gene, disease-specific BA1/BS1 cutoffs
consumed by `acmg_classifier.criteria.allele_frequency.DiseaseThresholds`
(used by the BA1 and BS1 evaluators). Genes absent from this file fall back to
the flat defaults (BA1 = 0.05, BS1 = 0.005).

## Schema (tab-separated, header required)

| column          | required | meaning |
|-----------------|----------|---------|
| `gene_symbol`   | yes      | HGNC symbol, the lookup key (matches VEP `gene_symbol`). |
| `inheritance`   | optional | `AD` / `AR` (substring match; `AR` → recessive √ formula, else dominant). |
| `prevalence`    | optional | Disease prevalence as a proportion (e.g. `0.0005` for 1/2000). |
| `allelic_het`   | optional | Max allelic contribution (single-variant), 0–1. Blank → 1.0. |
| `genetic_het`   | optional | Max genetic contribution (this gene), 0–1. Blank → 1.0. |
| `penetrance`    | optional | Penetrance 0–1. Blank → cannot compute (falls back to default). |
| `bs1_threshold` | optional | **Direct override.** If set, used verbatim for BS1. |
| `ba1_threshold` | optional | **Direct override.** If set, used verbatim for BA1. |
| `source_vcep`   | optional | Provenance (e.g. `RASopathy VCEP v2.1`). Not read by the tool. |
| `cspec_url`     | optional | Link to the criteria specification. Not read by the tool. |
| `notes`         | optional | Free text. Not read by the tool. |

Resolution order **per criterion, per gene**:
1. `*_threshold` override column → use as-is.
2. Else compute from biological parameters (Whiffin/Ware 2017):
   - `G = prevalence × genetic_het × allelic_het / penetrance`
   - dominant: `maxAF = G / 2` ; recessive: `maxAF = √G`
   - `BS1 = max(maxAF, 0.0005)` (0.05% floor)
   - `BA1 = min(0.05, 10 × maxAF)` (5% ceiling)
3. Else flat default (BS1 0.005 / BA1 0.05).

## Recommended: auto-generate from cspec GN*.json exports

The cspec registry exports each spec as machine-readable JSON-LD
(`GNxxx.json`), containing the gene(s), mode(s) of inheritance, and the BA1/BS1
filtering-AF thresholds. To build the whole table automatically:

1. From <https://cspec.genome.network/>, export each spec as JSON into this
   directory (`resources/clingen/GN0001.json`, `GN0115.json`, …). The JSON is
   also available from the API, e.g.
   `https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/GN115`.
2. Generate the table:
   ```bash
   python scripts/build_disease_thresholds.py \
       --json-dir resources/clingen \
       --out resources/clingen/disease_prevalence.tsv --released-only
   ```
   One row per gene, BA1/BS1 taken from the *Applicable* evidence-strength
   descriptions. Rows whose description packs multiple cutoffs (some multi-MOI
   specs) are flagged `verify` in `notes` and keep the most conservative
   (highest) value — spot-check those. PM2 is not emitted (the tool derives PM2
   from raw AF separately).

This regenerates the whole file from JSON, so it overwrites any hand-curated
rows — keep curation in the JSON source or in override columns of dedicated rows.

## How to populate (authoritative, verified)

**These are clinical thresholds — do not invent values.** Each row must trace
to a published source.

1. **Preferred — VCEP-published cutoffs (override columns):**
   ClinGen Criteria Specification Registry → <https://cspec.genome.network/>.
   Open each VCEP's rule specifications, read the BA1 / BS1 gnomAD thresholds,
   and the in-scope gene list. Put the published numbers in
   `bs1_threshold` / `ba1_threshold`, one row per gene, and record
   `source_vcep` + `cspec_url`. Many VCEPs apply one threshold to all genes in
   scope; some are gene-specific — follow the spec.

2. **Computed (when a VCEP gives parameters, not a final cutoff):**
   Fill `prevalence` (Orphanet / GeneReviews), `inheritance` (OMIM / PanelApp),
   `penetrance` (GeneReviews; default 1.0 if unknown), and the heterogeneity
   factors when available. The tool (or the cardiodb Allele Frequency App,
   <https://www.cardiodb.org/allelefrequencyapp/>) derives maxAF.

Leave heterogeneity columns blank to assume 1.0 (most conservative — highest
threshold, hardest to fire BA1/BS1).

## Deploying to a data directory

The runtime path is `<data-dir>/shared/disease_prevalence.tsv`
(`Config.disease_prevalence_tsv`). Copy this curated file there:

```bash
mkdir -p <data-dir>/shared
cp resources/clingen/disease_prevalence.tsv <data-dir>/shared/
```

## Provenance / distribution

Values are transcribed from publicly available ClinGen VCEP criteria
specifications (cspec.genome.network) and standard epidemiology resources
(Orphanet, OMIM, GeneReviews). Keep `source_vcep` / `cspec_url` populated so
every threshold remains auditable and attributable to its VCEP.
