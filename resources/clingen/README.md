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
| `af_basis`      | optional | Which gnomAD frequency BA1/BS1 compare against: `males` → the **male (XY) allele frequency** (gnomAD `AF_XY`) for X-linked "in males" genes (RPGR, RS1, ABCD1, SLC6A8, OTC); `popmax` → the **point grpmax/popmax allele frequency** for VCEPs that define the cutoff on the point AF rather than the 95%-CI FAF (RUNX1, GAA, MYOC, GAMT, BMPR2, MECP2, FOXG1, TCF4, UBE3A, RYR1, USH2A, CDH1, F8, FBN1, SCN2A); blank → overall FAF95. `ACMG_POPMAX_AF_BASIS=false` forces every `popmax` gene back to FAF95. Falls back to overall FAF when the gnomAD DB predates the needed column. |
| `ba1_hom_count` | optional | Homozygote+hemizygote count at/above which BA1 fires regardless of frequency (VCEP OR-clause: SLC6A8/OTC `10`). Blank → frequency-only. |
| `pp2`           | optional | VCEP PP2 decision: `applicable` → PP2 fires for missense in this gene; `not_applicable` → PP2 suppressed (VCEP declined it). Blank → no VCEP covers the gene; PP2 uses its statistical heuristic. |
| `pp2_requires`  | optional | Co-criteria PP2 is conditional on, comma-separated (e.g. `PM2,PP3` for BMPR2). PP2 is suppressed post-hoc unless all listed criteria are also triggered for the variant. Blank → unconditional. |
| `pvs1`          | optional | VCEP PVS1 decision: `not_applicable` → PVS1 suppressed because loss-of-function is not the disease mechanism (gain-of-function / dominant-negative genes: MYOC, RASopathy, cardiomyopathy, the activating PIK3 genes, VWF, …). `applicable` → the VCEP explicitly applies PVS1, which **establishes LoF as the mechanism** (ACTA1/RYR1 Congenital Myopathies): the decision tree then skips the ClinVar/LOEUF heuristic and the undercuration strength caps. Blank → the decision tree runs with the heuristic (APC has a gene-specific tree in code). |
| `bs1_exclude`   | optional | Bare protein change(s) the VCEP bars from BS1 regardless of frequency — a recurrent disease allele (e.g. MYOC `p.Gln368Ter`). Comma-separated; `*` and `Ter` compare equal. |
| `pm2_subpop`    | optional | Highest-subpopulation correction for the deflated low-AC FAF95: `point` (RUNX1 — also require the GrpMax point AF ≤ threshold) or `ci95` (Cardiomyopathy/HCM, LGMD — require the upper 95% CI of the GrpMax AF ≤ threshold, from gnomAD `ac_grpmax`/`an_grpmax`). Blank → FAF/AF as-is. |
| `pm2_zygosity`  | optional | Homozygote/hemizygote ceiling `<scope>:<max>` (`hom` / `hemi` / `homhemi`) that PM2 also requires (SLC6A8 `homhemi:0`, OTC `homhemi:1`, the SCID genes / GATM / GAMT `hom:0`, ABCD1 `hemi:0`). Blank → no zygosity gate. |
| `pm2_subset`    | optional | `non_cancer` → judge PM2 absence on the gnomAD non-cancer subset (ENIGMA BRCA1/2). Blank → overall callset. |
| `pm2_min_depth` | optional | Minimum gnomAD mean read depth for "absent" to be callable (ENIGMA BRCA1/2 `25`); requires the gnomAD coverage DuckDB. Blank → no depth gate. |
| `bp7_phylop`    | optional | BP7 phyloP "highly conserved" cutoff (e.g. `0.1` / `0.2` / `0` / default `2.0`), or `na` when the VCEP declared conservation non-informative (TP53, BRCA1/2, RUNX1, MYOC, the LCA & SCID genes — BP7 then skips the conservation gate). phyloP100way. Blank → global `bp7_phylop_max` default. |
| `bp7_intronic`  | optional | `noncanonical` → BP7 applies to any intronic position except the canonical ±1,2 (RASopathy / PIK3 panels, RUNX1, MYOC, VHL). Blank → the Walker deep-intronic (+7/−21) default. |
| `revel_pp3_supporting` / `revel_pp3_moderate` / `revel_pp3_strong` | optional | VCEP gene-specific REVEL **PP3** cutoffs (`REVEL ≥ value` fires at that strength). Only consulted when `--insilico-tool revel`. A spec that grants PP3 only at Supporting fills just `*_supporting`, capping the gene at Supporting. All blank → the genome-wide Bergquist 2024 default tiers apply. |
| `revel_bp4_supporting` / `revel_bp4_moderate` / `revel_bp4_strong` | optional | VCEP gene-specific REVEL **BP4** cutoffs (`REVEL ≤ value` fires at that strength). Same single-tier capping and default-fallback rules as the PP3 columns. |
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
- **`pp2`** records each VCEP's gene-level PP2 decision from the PP2 criteria
  code: `applicable` when a PP2 strength is Applicable, `not_applicable` when the
  VCEP carries a PP2 code but declined it, blanket-negated the description
  (KCNQ1: "Not applicable due to … z-score 1.83"), or excluded the gene by name
  (GN018: "applicable to MTOR, PIK3CA and AKT3 but **not PIK3R2**"). Resolved to
  the **most gene-specific** spec (a single-gene VCEP supersedes a grouped
  panel); on a specificity **tie** across distinct diseases the conservative
  `not_applicable` wins (it suppresses PP2 → fewer false pathogenic-supporting
  calls). A gene whose disease-appropriate decision is `applicable` but loses a
  tie (e.g. ACTA1: GN147 applicable vs GN169 not-applicable) can be pinned with
  `--override ACTA1:pp2=applicable`. This authoritative list is what the PP2
  evaluator uses to avoid over-assigning the (gene-level) criterion; genes no
  VCEP covers fall back to a ClinVar/gnomAD statistical heuristic.
- **`pp2_requires`** captures co-criteria a VCEP makes PP2 conditional on
  ("PM2_supporting and PP3 must be met" → `PM2,PP3` for BMPR2). The registry
  suppresses PP2 post-hoc unless every listed criterion is also triggered.
- **`revel_pp3_*` / `revel_bp4_*`** are mined from the PP3 / BP4 criteria-code
  descriptions: each *Applicable* strength tier (Supporting / Moderate / Strong;
  Very Strong folds into Strong) contributes its REVEL cutoff. Numbers are read
  only from REVEL-anchored text (other tools' cutoffs — SpliceAI, CADD,
  AlphaMissense — are skipped), a range/band yields the **firing edge** (lower
  bound for PP3, upper for BP4), and a monotonicity guard drops contradictory
  tiers (e.g. a Moderate cutoff below Supporting from a data-entry typo).
  Resolved to the **most gene-specific** spec, like PM2. Specs that cite REVEL
  without a number (e.g. SCN1A "follow ClinGen recommendations") leave the
  columns blank → the genome-wide Bergquist 2024 default applies.

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
