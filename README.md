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
- [Testing](#testing)
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
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh38 --gnomad-workers 6
```

**Step 4: Verify installation**

```bash
acmg-classify status --data-dir /path/to/download/directory/data
```

**Step 5: Classify variants**

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
- **In silico prediction**: AlphaMissense (default, non-commercial) or
  **ESM1b** (MIT-licensed, commercial-use ready) for missense; SQUIRLS
  (default) or SpliceAI for splice. SpliceAI overrides the missense call
  when its score crosses the high-impact threshold. Both missense
  predictors use Bergquist 2024 Table 2 strengths.
- **Fully offline** at classification time. Local databases include:
  Ensembl reference genome, VEP cache, gnomAD (DuckDB), ClinVar (VCF +
  derived SQLite for PS1/PM5), AlphaMissense, RepeatMasker, optional SpliceAI.
- **Manual evidence supplement**: per-variant TSV overrides for PS3/PP1/PM3/
  etc. that the auto-pipeline cannot derive (functional studies, family
  segregation, etc.).
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
| Disk | ~ 350 GB free | gnomAD v4.1 exomes dominate the footprint |
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
- **java** (>= 11) — only needed if SQUIRLS is used as the splice predictor

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

# Optional: dev extras (pytest, ruff)
pip install -e ".[dev]"

# 4. Verify the CLI is wired up
acmg-classify --help
```

The console script `acmg-classify` is registered by `pyproject.toml` →
`[project.scripts]` and resolves to `acmg_classifier.cli:cli`.

---

## Data setup

All annotation databases live under a single `data/` directory (NOT committed
to git). A one-shot setup script downloads and builds everything required:

```bash
# Default: GRCh38 only, downloads everything (~ 350 GB, takes hours)
python scripts/setup_data.py --data-dir ./data --assembly GRCh38

# Both assemblies
python scripts/setup_data.py --data-dir ./data --assembly GRCh38
python scripts/setup_data.py --data-dir ./data --assembly GRCh37

# If you already have a reference genome / gnomAD VCFs locally, point to them
python scripts/setup_data.py --data-dir ./data \
    --genome-fasta /db/reference/GRCh38/hg38.fa \
    --gnomad-vcf-dir /db/gnomad/v4.1/exomes/vcf

# Skip the gigantic gnomAD download (~300 GB)
python scripts/setup_data.py --data-dir ./data --skip-gnomad

# Pick specific chromosomes for gnomAD (testing/partial setup)
python scripts/setup_data.py --data-dir ./data \
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
| `--gnomad-workers N` | DuckDB build parallelism (default = CPU - 1) |
| `--skip-gnomad` | Skip gnomAD download (~ 300 GB) |
| `--skip-genome` | Skip reference FASTA download (~ 880 MB) |
| `--skip-vep-cache` | Skip VEP cache download (~ 14 GB) |
| `--skip-esm1b` | Skip ESM1b LLR archive download / SQLite build (~ 1.34 GB) |

After setup, the expected layout is:

```
data/
├── shared/
│   └── gene_inheritance.tsv         # gene → AD/AR/XL (ships in repo)
├── vep_cache/                       # VEP indexed cache, both assemblies
├── esm1b/                           # (optional) protein-coordinate, shared across assemblies
│   └── esm1b_llr.sqlite             #   built from Brandes 2023 archive
└── GRCh38/                          # (mirror at GRCh37/)
    ├── genome/GRCh38.p14.fa(+.fai)
    ├── clinvar/clinvar_GRCh38.vcf.gz(+.tbi)
    ├── clinvar/clinvar_ps1_pm5_GRCh38.sqlite
    ├── gnomad/gnomad_v4.1_exomes.duckdb
    ├── gnomad/gnomad_v4.1_constraint.tsv
    ├── alphamissense/AlphaMissense_hg38.tsv.gz
    ├── repeats/repeatmasker_dfam_hg38.bed.gz
    └── (optional) spliceai/spliceai_scores.raw.{snv,indel}.hg38.vcf.gz
```

Validate the install at any time:

```bash
acmg-classify validate --data-dir ./data --assembly GRCh38
acmg-classify status   --data-dir ./data
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
| `--insilico-tool {alphamissense,esm1b}` | str | `alphamissense` | Missense predictor used for PP3/BP4 |
| `--splice-tool {squirls,spliceai}` | str | `squirls` | Splice predictor used for PP3/BP4. SpliceAI takes precedence over AlphaMissense when its score ≥ 0.20. |
| `--spliceai-dir PATH` | path | `<data-dir>/<asm>/spliceai/` | Override SpliceAI VCF directory |
| `--supplement PATH` | path | — | Manual evidence TSV (see below) |
| `--inheritance {AD,AR,XL,Unknown}` | str | `Unknown` | Pedigree-level inheritance default (per-gene values in `data/shared/gene_inheritance.tsv` take precedence) |
| `--workers N` | int | `4` | Parallel workers for annotation |

Example for a panel sequencing run with manual segregation evidence:

```bash
acmg-classify classify panel.vcf.gz \
    -o panel_results.tsv \
    --assembly GRCh38 \
    --insilico-tool alphamissense \
    --splice-tool squirls \
    --supplement manual_evidence.tsv \
    --inheritance AR \
    --workers 8
```

### `explain` — single-variant detail

Print a human-readable breakdown of one variant without writing to disk:

```bash
acmg-classify explain chr17 43044295 G A --assembly GRCh38
```

Useful for debugging a specific classification or for ad-hoc review.

### `validate` / `status` / `setup`

```bash
# Check that all required data files exist
acmg-classify validate --data-dir ./data --assembly GRCh38
acmg-classify validate --data-dir ./data --assembly both

# Show DB versions / build dates
acmg-classify status   --data-dir ./data

# (Re)run downloads — same effect as scripts/setup_data.py without
# the advanced flags
acmg-classify setup    --data-dir ./data --assembly GRCh38
```

### Manual evidence supplement

A tab-separated file with these columns:

```
variant_id	criterion	strength	evidence
chr17:43044295:G:A	PS3	Strong	PMID:12345678 functional study showing loss of BRCA1 HR activity
chr17:43044295:G:A	PP1	Supporting	3 affected family members segregating variant
```

- `variant_id` must match the canonical `chrom:pos:ref:alt` used in the output.
- `criterion` is any ACMG code (e.g. `PS3`, `PP1`, `BS3`, `BP5`).
- `strength` is one of `VeryStrong`, `Strong`, `ThreePoint`, `Moderate`,
  `Supporting`.
- `evidence` is a free-text rationale shown in the `*_evidence` column.

Manual entries **override** automated evaluations for the same variant /
criterion pair (auto-evaluated PS3 from ClinVar SCV is replaced by the manual
PS3 line, for example). This is the recommended path for incorporating
unpublished functional data, segregation analysis, or expert-panel
de novo assessments (PS2/PM6).

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
−5..+5 → VUS
−9..−6 → Likely Benign
≤ −10  → Benign
```

`BA1` triggers stand-alone Benign regardless of the sum (per Tavtigian 2020).

### Annotation snapshot

`gnomad_af`, `gnomad_ac`, `gnomad_an`,
`gnomad_faf95_popmax`, `gnomad_popmax_af`, `gnomad_popmax_pop`,
`gnomad_pli`, `gnomad_loeuf`,
`clinvar_variation_id`, `clinvar_significance`, `clinvar_stars`,
`alphamissense_score`, `alphamissense_classification`,
`esm1b_llr`,
`splice_tool`, `splice_score`,
`in_repeat`, `repeat_class`

`alphamissense_*` is populated when `--insilico-tool alphamissense` (default);
`esm1b_llr` is populated when `--insilico-tool esm1b`. The other column is
left empty for the non-active tool.

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

- **AlphaMissense** (default): Strong / ThreePoint / Moderate / Supporting
  for PP3 (`≥0.990 / ≥0.972 / ≥0.906 / ≥0.792`); ThreePoint / Moderate /
  Supporting for BP4 (`≤0.070 / ≤0.099 / ≤0.169`). No Strong BP4 category.
- **ESM1b** (`--insilico-tool esm1b`): Strong / ThreePoint / Moderate /
  Supporting for PP3 (LLR `≤−24.0 / ≤−14.0 / ≤−12.2 / ≤−10.7`); ThreePoint /
  Moderate / Supporting for BP4 (LLR `≥8.8 / ≥−3.2 / ≥−6.3`). Lower LLR ⇒
  more pathogenic. Use this path for **commercial deployments** — see
  [Commercial use](#commercial-use).

**Splice predictor** — default **SQUIRLS**. When `--splice-tool spliceai` or
pre-computed SpliceAI VCFs are present at `data/<asm>/spliceai/`, SpliceAI
takes precedence over the missense predictor when its max Δscore ≥ 0.20.
Below that threshold the missense predictor's call is retained.

### Commercial use

The tool itself is Apache-2.0, but **AlphaMissense scores are CC BY-NC-SA 4.0
(non-commercial)**. For commercial deployments, switch the missense
predictor to **ESM1b** (Brandes et al. 2023, MIT-licensed), which has full
Bergquist 2024 strength calibration in this implementation:

```bash
acmg-classify classify input.vcf -o results.tsv \
    --assembly GRCh38 --data-dir ./data \
    --insilico-tool esm1b --splice-tool squirls
```

All other defaults (gnomAD, ClinVar, VEP, SQUIRLS) are commercially
permissive. SpliceAI remains a separate Illumina-licensed option.

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
│   └── setup_data.py                # Database download & build script
├── tests/
│   ├── unit/                        # pytest unit tests
│   ├── integration/                 # end-to-end pipeline tests (require data/)
│   └── fixtures/                    # sample.vcf, sample_supplement.tsv
├── data/                            # (NOT versioned) annotation databases
│   └── shared/gene_inheritance.tsv  # tracked: gene → AD/AR/XL map
├── pyproject.toml
├── LICENSE
├── NOTICE
└── README.md
```

---

## Testing

```bash
# Install dev extras (pytest, pytest-cov)
pip install -e ".[dev]"

# Unit tests only (no data/ required)
pytest tests/unit -v

# With coverage
pytest tests/unit --cov=src/acmg_classifier --cov-report=term-missing

# Full integration suite (requires populated data/)
pytest tests/integration -v
```

Lint:

```bash
ruff check src tests
```

---

## Configuration via environment variables

`Config` (Pydantic `BaseSettings`) reads any `ACMG_*` environment variable
or a `.env` file in the working directory. Examples:

```bash
export ACMG_DATA_DIR=/db/acmg
export ACMG_ASSEMBLY=GRCh38
export ACMG_WORKERS=8
export ACMG_INSILICO_TOOL=alphamissense
export ACMG_SPLICE_TOOL=squirls
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
- **AlphaMissense license.** Scores are CC BY-NC-SA 4.0 — commercial use
  requires direct arrangement with DeepMind/Google. The tool itself is
  Apache-2.0 but the bundled annotation source is not. Switch to
  `--insilico-tool esm1b` for commercial settings (see
  [Commercial use](#commercial-use)).
- **SpliceAI.** Pre-computed score VCFs are not redistributed. Users with an
  Illumina license can place them under `data/<asm>/spliceai/` and pass
  `--splice-tool spliceai`.
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
