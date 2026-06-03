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
| `af_basis`      | optional | `males` → BA1/BS1 compare against the **male (XY) allele frequency** (gnomAD `AF_XY`) instead of the overall population FAF. For X-linked genes whose VCEP states the cutoff "in males" (RPGR, RS1, ABCD1, SLC6A8, OTC). Blank → overall FAF95. Falls back to overall FAF when the gnomAD DB predates the `af_xy` column. |
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

Independently, `af_basis` selects **which** gnomAD frequency the resolved cutoff
is compared against: `males` → `AF_XY` (X-linked "in males" genes), blank →
overall FAF95.

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
   descriptions. PM2 is not emitted (the tool derives PM2 from raw AF
   separately).

### How the threshold is parsed from each description

The free-text BA1/BS1 descriptions are not uniform, so the extractor applies a
few precedence rules (all verified against the released specs):

- **Sub-population rule wins.** `"present at ≥X in any sub-population"` is the
  VCEP's operative gnomAD cutoff and overrides a generic `"above 0.05%"`
  headline (Rett/Angelman-like panels → BA1 `0.000083`, not `0.0005`).
- **Legacy 5% dropped.** `"above 5% in ESP / 1000 Genomes / ExAC"` is pre-gnomAD
  boilerplate; it is ignored when a gnomAD-specific number is also present
  (KCNQ1 BA1 → `0.004`, not `0.05`), but kept when it is the only value
  (RPGR-style males "5%").
- **Range bands take the lower edge.** `"between X and Y"` → the BS1 cutoff is
  `min(X, Y)`, independent of which bound is written first (RPE65, RUNX1).
- **`af_basis=males`** is set when the description says "in males"/"hemizygous".

### Multi-spec genes and `--override`

When a gene appears in several specs, the **more gene-specific** spec wins (a
single-gene VCEP supersedes a grouped panel — fixes FOXG1/MECP2/etc.). On a
specificity tie across distinct diseases (e.g. RYR1: Malignant Hyperthermia vs
Congenital Myopathies; ACTA1), the conflict cannot be auto-resolved, so the
build defaults to the **most conservative** cutoff — the highest BA1 (then
highest BS1) — which **minimises false-positive benign calls**. Such rows are
flagged in `notes`; pin a disease-appropriate value with `--override`:

```bash
python scripts/build_disease_thresholds.py \
    --json-dir resources/clingen \
    --out resources/clingen/disease_prevalence.tsv --released-only \
    --override "RYR1:ba1=0.0038,bs1=0.0007"
```

`--override GENE:field=val[,field=val]` is repeatable, applied after multi-spec
resolution, accepts `ba1` / `bs1` / `af_basis` / `inheritance`, can add a gene
absent from every spec, records `manual override` in `notes`, and fails loudly
on an unknown field.

> **X-linked "in males" genes need gnomAD `AF_XY`.** Rows with `af_basis=males`
> are compared against the male (XY) allele frequency. That requires the gnomAD
> DuckDB to carry the `af_xy` column (built from VCF `AF_XY`); a DB built before
> that column was added still works — the evaluators fall back to the overall
> FAF. Rebuild the gnomAD DB (`scripts/setup_data.py`) to activate male-AF.

This regenerates the whole file from JSON, so it overwrites any hand-curated
rows — keep curation in the JSON source, the `--override` flags, or override
columns of dedicated rows.

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
