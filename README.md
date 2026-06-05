# HUHVar

**Hiroshima University Hospital Variant classification tool** (package name:
`acmg-classifier`)


[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-%3E=3.11-blue.svg)](pyproject.toml)
[![Status](https://img.shields.io/badge/status-research_use_only-orange.svg)]()

**Fully local ACMG/AMP variant pathogenicity classifier** implementing both the
original 2015 combinatorial rules (Richards et al., *Genet Med* 2015) and the
Bayesian point-based framework (Tavtigian et al., *Genet Med* 2018/2020;
Bergquist et al., 2024).

The tool annotates each variant in an input VCF using only on-disk databases
(no network calls at classification time), evaluates all 28 ACMG criteria, and
writes a TSV containing per-criterion strength, supporting evidence, and the
final classification under both frameworks side-by-side.

> ⚠️ **Research use only.**
> **This software is intended for research use only and not for clinical
> diagnostic purposes.** It is not a medical device and must not be used for
> primary clinical diagnosis without independent expert review.
> See [Disclaimer](#disclaimer).

---

## Table of Contents

- [Quick start](#quick-start)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Data setup](#data-setup)
- [Usage](#usage)
  - [classify — bulk VCF classification](#classify--bulk-vcf-classification)
  - [explain — single-variant detail](#explain--single-variant-detail)
  - [validate / status / setup](#validate--status--setup)
  - [Manual evidence supplement](#manual-evidence-supplement)
- [Output format](#output-format)
- [Classification model](#classification-model)
- [Commercial use](#commercial-use)
- [Project layout](#project-layout)
- [Configuration via environment variables](#configuration-via-environment-variables)
- [Known limitations](#known-limitations)
- [Citing](#citing)
- [Contributing](#contributing)
- [License](#license)
- [Disclaimer](#disclaimer)

---

## Quick start

**Step 1: Create conda environment and install dependencies**

```bash
conda create -n acmg python=3.12
conda activate acmg
conda install -c bioconda -c conda-forge samtools tabix bcftools htslib ensembl-vep=111
```

**Step 2: Clone and install HUHVar**

```bash
git clone https://github.com/hnakahara/HUHVar.git
cd HUHVar
pip install -e .
```

**Step 3: Download and build annotation databases**

> ⏱️ This step downloads ~350 GB and builds local indexes. **Expect 1–2 days
> to complete** depending on network throughput and CPU; run it in `tmux` /
> `screen` or as a background job.

```bash
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh38 --workers 12
```

**Step 4: Install the curated gene-rule tables**

Copy the VCEP-derived gene rule tables into `data/shared/` (the criteria
evaluators read them from there):

```bash
mkdir -p /path/to/download/directory/data/shared
cp resources/gene_inheritance.tsv \
   resources/clingen/disease_prevalence.tsv \
   resources/clingen/pm1_hotspots.tsv \
   /path/to/download/directory/data/shared/
```

**Step 5: Verify installation**

```bash
acmg-classify status --data-dir /path/to/download/directory/data
```

**Step 6: Classify variants**

```bash
acmg-classify classify input.vcf -o results.tsv --assembly GRCh38 --data-dir /path/to/download/directory/data --workers 12
```

---

## Features

- **Two classifiers in parallel**: ACMG 2015 Table 5 combinatorial rules and
  Tavtigian 2020 Bayesian point system reported side-by-side, so reviewers can
  spot disagreements at a glance.
- **All 28 ACMG criteria** auto-evaluated (PVS1, PS1–PS4, PM1–PM6, PP1–PP5,
  BA1, BS1–BS4, BP1–BP7), with ClinGen SVI strength adjustments and the
  Bergquist-2024 3-point extension.
- **PVS1 decision tree** (Abou Tayoun et al. 2018) applied for LoF variants
  including NMD prediction, last-exon rescue, and biological-relevance gating.
- **Inheritance-aware PM2** (BS1/BS2 also) thresholds switch between dominant
  and recessive frequencies using a per-gene inheritance table.
- **Per-gene ClinGen VCEP rules** mined from the cspec specifications
  (`resources/clingen/disease_prevalence.tsv`, `pm1_hotspots.tsv`): disease-
  specific BA1/BS1 cutoffs, PP2 applicability, PM5 Grantham-distance gating,
  PM1 hotspot regions, inheritance-aware BS2 (incl. a dominant heterozygote
  rule), the PS1 splice extension (canonical vs non-canonical), and BP1/BP3
  applicability. X-linked "in males" genes are compared against the gnomAD
  male (XY) allele frequency. See `resources/clingen/README.md`.
- **gnomAD frequencies** use the v4.1 **joint** (combined exome+genome) release
  on GRCh38; on GRCh37 (no joint release) the exome and genome callsets are both
  loaded and merged per-field, and BA1/BS1/PM2 compare against the GrpMax
  filtering allele frequency (FAF95).
- **In silico prediction**: **ESM1b** (default; Brandes 2023, MIT-licensed,
  commercial-use ready) or **AlphaMissense** (non-commercial) for missense.
  Splice prediction defaults to **OpenSpliceAI** (Chao 2025, GPL-3.0; runtime
  inference via the `openspliceai` CLI); **SpliceAI** (Illumina-licensed) is an
  opt-in alternative. The splice predictor overrides the missense call —
  including on missense variants — when its score crosses the high-impact
  threshold (≥ 0.20). The missense predictors use Bergquist 2024 Table 2
  strengths.
  > _OpenSpliceAI reuses SpliceAI's 0–1 delta scale; lacking an OddsPath
  > calibration of its own, its PP3 is awarded at the conservative Supporting
  > tier (vs SpliceAI's Moderate). SQUIRLS and MMSplice are retained in the
  > code but disabled (SQUIRLS' precomputed DB is no longer downloadable;
  > MMSplice has a dependency conflict)._
- **Fully offline** at classification time (OpenSpliceAI runs locally too).
  Local databases include: Ensembl reference genome, VEP cache, gnomAD
  (DuckDB), ClinVar (VCF + derived SQLite for PS1/PM5), ESM1b / AlphaMissense,
  RepeatMasker, and optional SpliceAI (opt-in, Illumina-licensed).
- **Manual evidence supplement**: per-variant TSV that can **override any
  criterion** (functional studies, family segregation, expert strength
  adjustments, etc.). Two combination modes (`--supplement-mode`): `merge`
  (default — curator entries override/add on top of the tool's calls) or
  `manual-only` (listed variants are classified purely from the supplement;
  variants not listed fall back to the tool's automated calls).
- **GRCh37 and GRCh38** both supported with separate database trees.
- **Deterministic / reproducible**: no remote API calls, no randomness in the
  scoring logic; output is byte-stable for a fixed (code, data) pair.

---

## Requirements

### System

| Component | Version | Notes |
|-----------|---------|-------|
| Linux / WSL2 | — | Tested on Rocky Linux 9.6; Windows native is **not** supported (uses bash & POSIX tools) |
| Python | ≥ 3.11 | 3.12 recommended |
| Disk | ~ 350 GB free | gnomAD dominates the footprint (v4.1 joint on GRCh38; exomes + genomes on GRCh37) |
| RAM | ≥ 16 GB | DuckDB build of gnomAD spikes briefly to ~10 GB |
| CPU | ≥ 4 cores | `--workers` defaults to 4 |

### External binaries (must be in `PATH`)

Install via `conda` / `mamba` (recommended):

```bash
conda create -n acmg -c bioconda -c conda-forge \
    python=3.12 \
    samtools tabix bcftools htslib \
    ensembl-vep=111
conda activate acmg
```

Required:

- **samtools** (>= 1.17) — FASTA indexing
- **tabix / bgzip** (htslib) — VCF indexing
- **bcftools** (>= 1.17) — VCF normalisation
- **vep** (Ensembl Variant Effect Predictor 111) — transcript annotation

Optional:

- **wget** or **curl** — data download (one of them is required)

### Python dependencies

Declared in `pyproject.toml`:

```
click pydantic pydantic-settings structlog rich tenacity
cyvcf2 pysam duckdb
```

`pip install -e .` installs everything in one step.

---

## Installation

```bash
# 1. Clone
git clone https://github.com/hnakahara/HUHVar.git
cd HUHVar

# 2. Create & activate environment (conda example)
conda create -n acmg -c bioconda -c conda-forge \
    python=3.12 samtools tabix bcftools htslib ensembl-vep=111
conda activate acmg

# 3. Install the package (editable install for development)
pip install -e .

# Optional: install the default splice predictor (OpenSpliceAI, GPL-3.0).
# Required for --splice-tool openspliceai (the default); skip to run without
# splice scoring or to use --splice-tool spliceai instead.
pip install -e ".[openspliceai]"
```

---

## Data setup

A one-shot setup script downloads and builds everything required. The gnomAD
per-chromosome VCFs are fetched from the **Google Cloud and AWS mirrors in
parallel** (one concurrent download from each, byte-identical), and interrupted
downloads **resume** on re-run (size-verified against the remote, so a partial
file left by Ctrl+C is completed rather than skipped):

```bash
# Default: GRCh38 only, downloads everything (~ 350 GB, takes hours)
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh38 --workers 12

# Both assemblies
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh38 --workers 12
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh37 --workers 12

# If you already have a reference genome / gnomAD VCFs locally, point to them
python scripts/setup_data.py --data-dir /path/to/download/directory/data \
    --genome-fasta /db/reference/GRCh38/hg38.fa \
    --gnomad-vcf-dir /db/gnomad/v4.1/joint/vcf \
    --workers 12

# Skip the gigantic gnomAD download (~300 GB)
python scripts/setup_data.py --data-dir /path/to/download/directory/data --skip-gnomad --workers 12

# Pick specific chromosomes for gnomAD (testing/partial setup)
python scripts/setup_data.py --data-dir /path/to/download/directory/data \
    --gnomad-chromosomes chr1 chr2 chrX
```

### `setup_data.py` options

| Flag | Description |
|------|-------------|
| `--data-dir PATH` | Root data directory (default `./data`) |
| `--assembly {GRCh38,GRCh37}` | Genome assembly |
| `--genome-fasta PATH` | Use an existing FASTA instead of downloading |
| `--gnomad-vcf-dir PATH` | Use existing gnomAD `*.vcf.bgz` files |
| `--gnomad-chromosomes CHR ...` | Subset of chromosomes (default all 24) |
| `--workers N` | Build parallelism for the ClinVar (XML parse, max 24) and gnomAD (DuckDB) steps (default = CPU - 1) |
| `--skip-gnomad` | Skip gnomAD download (~ 300 GB) |
| `--skip-genome` | Skip reference FASTA download (~ 880 MB) |
| `--skip-vep-cache` | Skip VEP cache download (~ 14 GB) |
| `--skip-esm1b` | Skip ESM1b LLR archive download / SQLite build (~ 1.34 GB) |

After setup, the expected layout is:

```
data/
├── shared/                          # curated gene-rule tables (copied from resources/)
│   ├── gene_inheritance.tsv         #   gene → AD/AR/XL
│   ├── disease_prevalence.tsv       #   per-gene VCEP rules (BA1/BS1, PP2, PM5, BS2, PS1, BP1/BP3, …)
│   └── pm1_hotspots.tsv             #   per-gene PM1 hotspot residues/regions
├── vep_cache/                       # VEP indexed cache, both assemblies
├── esm1b/                           # (optional) protein-coordinate, shared across assemblies
│   └── esm1b_llr.sqlite             #   built from Brandes 2023 archive
└── GRCh38/                          # (mirror at GRCh37/)
    ├── genome/GRCh38.p14.fa(+.fai)
    ├── clinvar/clinvar_GRCh38.vcf.gz(+.tbi)
    ├── clinvar/clinvar_ps1_pm5_GRCh38.sqlite
    ├── gnomad/gnomad_v4.1_joint.duckdb   # GRCh37: gnomad_v2.1.1_exome_genome.duckdb
    ├── gnomad/gnomad_v4.1_constraint.tsv
    ├── alphamissense/AlphaMissense_hg38.tsv.gz
    ├── repeats/repeatmasker_dfam_hg38.bed.gz
    ├── (default splice) openspliceai/2000nt/   # OSAI_MANE model dir (see note below)
    └── (optional)        spliceai/spliceai_scores.raw.{snv,indel}.hg38.vcf.gz
```

> **OpenSpliceAI models are not downloaded by `setup_data.py`.** Install the
> CLI (`pip install openspliceai`) and place the OSAI_MANE model files under
> `data/<asm>/openspliceai/<flanking-size>nt/` (default `2000nt/`), or point
> `--openspliceai-model-dir` at them. Download models from
> <https://github.com/Kuanhao-Chao/openspliceai>. If neither the CLI nor a
> model directory is present, splice scoring is skipped (no splice evidence).

Validate the install at any time:

```bash
acmg-classify validate --data-dir /path/to/download/directory/data --assembly GRCh38
acmg-classify status   --data-dir /path/to/download/directory/data
```

---

## Usage

### `classify` — bulk VCF classification

```bash
acmg-classify classify input.vcf -o results.tsv --assembly GRCh38
```

Full option list:

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `VCF` (positional) | path | — | Input VCF (`.vcf`, `.vcf.gz`); must exist |
| `-o, --output PATH` | path | stdout | Output TSV path |
| `--data-dir PATH` | path | `./data` | Data root |
| `--assembly {GRCh37,GRCh38}` | str | auto-detect from VCF header → fallback `GRCh38` | Force a specific assembly |
| `--insilico-tool {esm1b,alphamissense}` | str | `esm1b` | Missense predictor used for PP3/BP4 |
| `--splice-tool {openspliceai,spliceai}` | str | `openspliceai` | Splice predictor used for PP3/BP4/PVS1. `openspliceai` (default, GPL-3.0) runs OSAI_MANE at inference time; `spliceai` (Illumina-licensed) uses precomputed VCFs. The splice call takes precedence over the missense call — including on missense variants — when its score ≥ 0.20. |
| `--openspliceai-model-dir PATH` | path | `<data-dir>/<asm>/openspliceai/<flank>nt/` | OSAI_MANE model directory |
| `--openspliceai-flanking-size N` | int | `2000` | Model context length (must match the downloaded model: 80/400/2000/10000) |
| `--spliceai-dir PATH` | path | `<data-dir>/<asm>/spliceai/` | Override SpliceAI VCF directory (only when `--splice-tool spliceai`) |
| `--supplement PATH` | path | — | Manual evidence TSV (see below) |
| `--supplement-mode {merge,manual-only}` | str | `merge` | How `--supplement` combines with the tool's calls (see [Manual evidence supplement](#manual-evidence-supplement)) |
| `--workers N` | int | `4` | Parallel workers for annotation |

Example for a panel sequencing run with manual segregation evidence:

```bash
acmg-classify classify panel.vcf.gz \
    -o panel_results.tsv \
    --assembly GRCh38 \
    --insilico-tool esm1b \
    --splice-tool openspliceai --openspliceai-flanking-size 2000 \
    --supplement manual_evidence.tsv --supplement-mode merge \
    --workers 12
```

### `explain` — single-variant detail

Print a human-readable breakdown of one variant without writing to disk:

```bash
acmg-classify explain chr17 7674221 G A --assembly GRCh38 --data-dir /path/to/download/directory/data
```

Useful for debugging a specific classification or for ad-hoc review.

### `validate` / `status` / `setup`

```bash
# Check that all required data files exist
acmg-classify validate --data-dir /path/to/download/directory/data --assembly GRCh38
acmg-classify validate --data-dir /path/to/download/directory/data --assembly both

# Show DB versions / build dates
acmg-classify status --data-dir /path/to/download/directory/data

# (Re)run downloads — same effect as scripts/setup_data.py without
# the advanced flags
acmg-classify setup --data-dir /path/to/download/directory/data --assembly GRCh38
```

### Manual evidence supplement

A tab-separated file with these columns:

```
variant_id	criterion	strength	evidence
chr17:43044295:G:A	PS3	Strong	PMID:12345678 functional study showing loss of BRCA1 HR activity
chr17:43044295:G:A	PP1	Supporting	3 affected family members segregating variant
```

- `variant_id` must match the canonical `chrom:pos:ref:alt` used in the output.
- `criterion` is any ACMG code (e.g. `PVS1`, `PS1`, `PS3`, `PM2`, `PP1`, `BP1`).
- `strength` is one of `VeryStrong`, `Strong`, `ThreePoint`, `Moderate`,
  `Supporting`.
- `evidence` is a free-text rationale shown in the `*_evidence` column.

Manual entries can override **any** criterion — including ones the tool
evaluates automatically (PVS1, PS1, PM2, BP1, …), not just the curation-only
ones. How they combine with the tool's calls is controlled by
`--supplement-mode`:

| Mode | Behaviour |
|------|-----------|
| `merge` (default) | Keep the tool's automated calls; for every criterion a curator names, **override its strength** (e.g. PVS1 Strong → Moderate) or **add it** if the tool left it not-met. Other criteria are untouched. |
| `manual-only` | For variants **listed** in the supplement, discard all automated evidence and classify **purely from the supplement**. Variants **not listed** fall back to the tool's automated calls. |

The override is applied **before** the ACMG combination rules (PVS1↔PP3
mutual exclusion, BA1/BS1/PM2 exclusivity, PP2 / PM5 gene gating), so those
rules operate on the curator-adjusted evidence. The audit trail records the
change (e.g. `[manual override Strong→Moderate] …`).

This is the recommended path for incorporating unpublished functional data,
segregation analysis, expert-panel de novo assessments (PS2/PM6), or expert
strength adjustments to any automated criterion.

A reference example lives at `tests/fixtures/sample_supplement.tsv`.

---

## Output format

The result TSV has the following columns (one row per variant):

### Identity & call

`variant_id`, `chrom`, `pos`, `ref`, `alt`, `filter`,
`transcript`, `gene`, `hgvs_c`, `hgvs_p`

### Classifications (both frameworks in parallel)

| Column | Meaning |
|--------|---------|
| `classification_2015` | One of `Pathogenic`, `LikelyPathogenic`, `VUS`, `LikelyBenign`, `Benign` — per ACMG 2015 Table 5 |
| `rules_2015` | The triggered criteria that drove the 2015 verdict (e.g. `PVS1 + PM2 + PP3`) |
| `bayesian_score` | Integer point total (Tavtigian 2020) |
| `classification_bayesian` | Bayesian verdict using the thresholds described below |

Bayesian thresholds:

```
≥ +10  → Pathogenic
+6..+9 → Likely Pathogenic
 0..+5 → VUS
−1..−6 → Likely Benign
≤ −7   → Benign
```

The benign side is asymmetric (prior P = 0.10): less evidence is needed to reach
benign, per Tavtigian 2020. `BA1` triggers stand-alone Benign regardless of the
sum.

### Annotation snapshot

`gnomad_af`, `gnomad_ac`, `gnomad_an`,
`gnomad_faf95_popmax`, `gnomad_popmax_af`, `gnomad_popmax_pop`,
`gnomad_pli`, `gnomad_loeuf`,
`clinvar_variation_id`, `clinvar_significance`, `clinvar_stars`,
`alphamissense_score`, `alphamissense_classification`,
`esm1b_llr`,
`splice_tool`, `splice_score`,
`in_repeat`, `repeat_class`

`esm1b_llr` is populated when `--insilico-tool esm1b` (default);
`alphamissense_*` is populated when `--insilico-tool alphamissense`. The other
column is left empty for the non-active tool. `splice_tool` / `splice_score`
reflect the active `--splice-tool` (`openspliceai` by default).

### Per-criterion columns

For each of the 28 ACMG codes (`PVS1`…`BP7`) three columns are emitted:

| Suffix | Content |
|--------|---------|
| (none) | `1` if triggered (and not suppressed), `0` otherwise |
| `_strength` | One of `VeryStrong`, `Strong`, `ThreePoint`, `Moderate`, `Supporting`, `NotMet`, `Indeterminate` |
| `_evidence` | Human-readable rationale (e.g. `gnomAD FAF95_popmax=0.0421 >= 0.05`, `ClinVar 2 stars: Pathogenic`) |

### Trailing

`warnings` — pipeline-level non-fatal issues, semicolon separated.

### Skipped records (sidecar)

Sites with `ALT='.'` (no variant) are filtered out of the main TSV and written
to a `*_skipped.tsv` alongside the main output, so input/output counts can be
reconciled.

---

## Classification model

### ACMG 2015 (rule-based)

Implements Table 5 of Richards et al. *Genet Med* 2015 with the following
ClinGen SVI updates:

- **PM2 downgraded to Supporting** (SVI recommendation, 2020)
- **PP5 / BP6 NOT applied** (SVI 2018: ClinVar reputable-source assertions
  decoupled from primary data risk double-counting with PS1/PS3/PM5)

See `src/acmg_classifier/classification/classifier_2015.py`.

### Bayesian point system (Tavtigian 2020 + Bergquist 2024)

Strength → points:

| Strength | Pathogenic | Benign |
|----------|------------|--------|
| Very Strong | +8 | −8 |
| Strong | +4 | −4 |
| Three Point (2024) | +3 | −3 |
| Moderate | +2 | −2 |
| Supporting | +1 | −1 |

Posterior probability of pathogenicity is implicit in the sum; the integer
thresholds match the natural-scale fit reported in the 2020 paper.

See `src/acmg_classifier/classification/classifier_bayesian.py`.

### PVS1 decision tree

The Abou Tayoun 2018 LoF decision tree is implemented at
`src/acmg_classifier/pvs1/` and decides whether PVS1 should fire as Very
Strong, Strong, Moderate, Supporting, or be entirely suppressed based on:

- NMD predicted vs escape
- Last-exon / final 50 nt of penultimate exon rules
- Rescue transcript existence
- Region of biological relevance
- Single-exon gene logic

### In-silico aggregation (PP3 / BP4)

Strengths are calibrated to Bergquist et al. *Genet Med* 2024 Table 2.

**Missense predictor** — pick one with `--insilico-tool`:

- **ESM1b** (default, `--insilico-tool esm1b`): Strong / ThreePoint / Moderate /
  Supporting for PP3 (LLR `≤−24.0 / ≤−14.0 / ≤−12.2 / ≤−10.7`); ThreePoint /
  Moderate / Supporting for BP4 (LLR `≥8.8 / ≥−3.2 / ≥−6.3`). Lower LLR ⇒
  more pathogenic. MIT-licensed → commercial-use ready (see
  [Commercial use](#commercial-use)).
- **AlphaMissense** (`--insilico-tool alphamissense`): Strong / ThreePoint /
  Moderate / Supporting for PP3 (`≥0.990 / ≥0.972 / ≥0.906 / ≥0.792`);
  ThreePoint / Moderate / Supporting for BP4 (`≤0.070 / ≤0.099 / ≤0.169`). No
  Strong BP4 category. Scores are CC BY-NC-SA 4.0 (non-commercial).

**Splice predictor** — default **OpenSpliceAI**. The splice call takes
precedence over the missense predictor when its max Δscore ≥ 0.20; below that
threshold the missense predictor's call is retained. The same Δscore also feeds
the PVS1 splice branch (≥ 0.20 supports splice-LoF).

- **OpenSpliceAI** (default): PP3 max_delta `≥0.20` → **Supporting**; BP4 /
  BP7 max_delta `≤0.10` → no-impact (Supporting). Shares SpliceAI's 0–1 delta
  scale, but lacking its own OddsPath calibration, PP3 is capped at Supporting
  (vs SpliceAI's Moderate). GPL-3.0; runtime inference via `openspliceai`.
- **SpliceAI** (`--splice-tool spliceai`): PP3 max_delta `≥0.20` → Moderate;
  BP4 max_delta `≤0.10` → Supporting (Walker *Am J Hum Genet* 2023).
- _SQUIRLS / MMSplice: retained in code but disabled (see overview note)._

### Commercial use

The tool itself is Apache-2.0, and the **defaults are commercial-use ready**:
ESM1b (MIT, missense) and OpenSpliceAI (GPL-3.0, splice) replace the
non-commercial AlphaMissense and the Illumina-licensed SpliceAI as defaults.

```bash
# This is the default configuration — shown explicitly for clarity.
acmg-classify classify input.vcf -o results.tsv \
    --assembly GRCh38 --data-dir /path/to/download/directory/data \
    --insilico-tool esm1b --splice-tool openspliceai
```

Notes:
- **AlphaMissense** (`--insilico-tool alphamissense`) is CC BY-NC-SA 4.0 —
  non-commercial only.
- **OpenSpliceAI** is GPL-3.0. Running it locally to produce classifications
  (a service/report) does not trigger GPL source-disclosure; redistributing
  the tool would. Do **not** use the Illumina SpliceAI model weights bundled in
  the openspliceai repo (CC BY-NC 4.0) for commercial work — use OSAI_MANE.
- **SpliceAI** (`--splice-tool spliceai`) remains a separate
  Illumina-licensed option.
- gnomAD, ClinVar, VEP are commercially permissive.

---

## Project layout

```
HUHVar/
├── src/acmg_classifier/
│   ├── cli.py                       # Click entry point (acmg-classify)
│   ├── config.py                    # Runtime configuration & derived paths
│   ├── pipeline/                    # Top-level orchestration
│   ├── annotation/                  # VEP / DB → AnnotationData aggregator
│   ├── local_db/                    # gnomAD / ClinVar / RepeatMasker / VEP runner
│   ├── pvs1/                        # PVS1 decision tree
│   ├── criteria/
│   │   ├── pathogenic/              # PVS1 / PS1-4 / PM1-6 / PP1-5
│   │   └── benign/                  # BA1 / BS1-4 / BP1-7
│   ├── classification/
│   │   ├── classifier_2015.py       # ACMG 2015 Table 5 combinatorial rules
│   │   └── classifier_bayesian.py   # Tavtigian 2020 point system
│   ├── models/                      # Pydantic data classes (VariantRecord, etc.)
│   ├── io/                          # VCF reader, TSV/report writer
│   ├── setup/                       # Programmatic setup / validate / status
│   └── utils/
├── scripts/
│   ├── setup_data.py                # Database download & build script
│   ├── build_disease_thresholds.py  # cspec GN*.json → disease_prevalence.tsv
│   └── build_pm1_hotspots.py        # cspec_summary.json → pm1_hotspots.tsv
├── resources/                       # tracked curated rule tables (copied to data/shared/)
│   ├── gene_inheritance.tsv         #   gene → AD/AR/XL
│   └── clingen/                     #   VCEP cspec exports + derived tables
│       ├── disease_prevalence.tsv   #     per-gene VCEP rules
│       └── pm1_hotspots.tsv         #     per-gene PM1 hotspots
├── tests/
│   ├── unit/                        # pytest unit tests
│   ├── integration/                 # end-to-end pipeline tests (require data/)
│   └── fixtures/                    # sample.vcf, sample_supplement.tsv
├── data/                            # (NOT versioned) annotation databases
├── pyproject.toml
├── LICENSE
├── NOTICE
└── README.md
```

---

## Configuration via environment variables

`Config` (Pydantic `BaseSettings`) reads any `ACMG_*` environment variable
or a `.env` file in the working directory. Examples:

```bash
export ACMG_DATA_DIR=/db/acmg
export ACMG_ASSEMBLY=GRCh38
export ACMG_WORKERS=8
export ACMG_INSILICO_TOOL=esm1b
export ACMG_SPLICE_TOOL=openspliceai
export ACMG_SUPPLEMENT_MODE=merge   # or manual-only

# Criterion tunables (sensible defaults; override only to trade precision/recall)
export ACMG_PM5_MIN_STARS=1        # min ClinVar review stars for a PM5 comparator
export ACMG_BS2_MIN_HOMALT=5       # BS2 healthy-homozygote count (recessive)
export ACMG_BS2_MIN_HEMI=5         # BS2 healthy-hemizygote count (X-linked)
export ACMG_BS2_MIN_HET=5          # BS2 healthy-carrier count (dominant)
```

CLI flags take precedence over environment variables.

---

## Known limitations

- **Windows native is not supported.** `setup_data.py` shells out to bash and
  htslib binaries. Use WSL2 or Linux.
- **gnomAD FILTER handling for BA1/BS1/BS2.** Records with a non-`PASS`
  FILTER currently cause the frequency-based benign criteria to be skipped.
  This is overly conservative for common variants and can lead to false LP
  calls. A fix is in the backlog; in the meantime cross-check any
  Bayesian-LP call where `gnomad_af` is high.
- **AlphaMissense license.** The non-default `--insilico-tool alphamissense`
  uses scores that are CC BY-NC-SA 4.0 — commercial use requires direct
  arrangement with DeepMind/Google. The default `esm1b` is MIT-licensed (see
  [Commercial use](#commercial-use)).
- **OpenSpliceAI models are a separate install.** The default `--splice-tool
  openspliceai` needs the `openspliceai` CLI (`pip install openspliceai`) and
  OSAI_MANE model files placed under `data/<asm>/openspliceai/<flank>nt/` (not
  fetched by `setup_data.py`). If absent, splice scoring is silently skipped
  (no splice evidence) and a warning is logged.
- **SpliceAI.** The non-default `--splice-tool spliceai` uses pre-computed
  score VCFs that are not redistributed. Users with an Illumina license can
  place them under `data/<asm>/spliceai/`.
- **No DUP/CNV support.** Only SNV / small INDEL / MNV are classified. SV
  callers should be paired with a dedicated CNV interpreter.
- **Single-sample input.** Multi-sample joint VCFs are accepted but only
  variant-level (not genotype-level) information is consumed.

---

## Citing

If you use this tool in published work, please cite the underlying ACMG
guidelines as well as this software:

```
Richards S, et al. Standards and guidelines for the interpretation of
sequence variants: a joint consensus recommendation of the American
College of Medical Genetics and Genomics and the Association for
Molecular Pathology. Genet Med. 2015;17(5):405-424.

Tavtigian SV, et al. Modeling the ACMG/AMP variant classification
guidelines as a Bayesian classification framework.
Genet Med. 2018;20(9):1054-1066.

Tavtigian SV, et al. Fitting a naturally scaled point system to the
ACMG/AMP-based variant classification framework. Hum Mutat. 2020.

HUHVar (acmg-classifier). Department of Clinical and Molecular Genetics,
Hiroshima University Hospital, 2026.
https://github.com/hnakahara/HUHVar
```

A machine-readable `CITATION.cff` is shipped with the repository for tools
that consume it (Zenodo, GitHub citation widget, etc.).

---

## Contributing

Issues and pull requests are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md)
for the development workflow. By contributing you agree that your
contributions will be licensed under Apache-2.0, per the [NOTICE](NOTICE) file.

---

## License

[Apache License 2.0](LICENSE) © 2026 Department of Clinical and Molecular
Genetics, Hiroshima University Hospital.

Third-party data sources accessed by `scripts/setup_data.py` are governed by
their own terms — see [NOTICE](NOTICE) for attribution and links.

---

## Disclaimer

**This software is intended for research use only and not for clinical
diagnostic purposes.**

HUHVar (`acmg-classifier`) is **not** an in vitro diagnostic device and is **not
certified for clinical decision making**. Any variant interpretation produced
by this tool must be reviewed by qualified clinical personnel before being
used to inform patient care.

The authors and the Department of Clinical and Molecular Genetics, Hiroshima
University Hospital make no warranty as to the accuracy or fitness for purpose
of the output, and accept no liability for clinical decisions made on the
basis of this software, in accordance with the "AS IS" provisions of the
Apache License 2.0.
