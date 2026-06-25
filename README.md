# HUVar

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

**Step 2: Clone and install HUVar**

```bash
git clone https://github.com/hnakahara/HUVar.git
cd HUVar
pip install -e .
```

**Step 3: Download and build annotation databases**

> ⏱️ This step downloads **~1.5 TB** (the gnomAD per-chromosome VCFs dominate)
> and builds local indexes. **Expect 1–2 days to complete** depending on network
> throughput and CPU; run it in `tmux` / `screen` or as a background job.
>
> 💾 The raw gnomAD VCFs are only needed to **build** the local DuckDB — once
> setup finishes you can reclaim most of that space by deleting them. See
> [Reclaiming disk space after setup](#reclaiming-disk-space-after-setup).

```bash
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh38 --workers 12
```

**Step 4: Install the curated gene-rule tables**

Copy the curated tables into `data/shared/` (the criteria evaluators read them
from there). This is two groups: the assembly-independent gene-rule tables in
`resources/shared/`, plus the **assembly-specific** eRepo manual-evidence
supplement from `resources/<assembly>/`:

```bash
mkdir -p /path/to/download/directory/data/shared

# 1) Assembly-independent gene-rule tables (auto-loaded by the evaluators)
cp resources/shared/*.tsv /path/to/download/directory/data/shared/

# 2) Manual-evidence supplement for your assembly (used via --supplement; see below)
#    GRCh38:
cp resources/GRCh38/erepo_manual_criteria_hg38.tsv /path/to/download/directory/data/shared/
#    GRCh37:
# cp resources/GRCh37/erepo_manual_criteria_hg19.tsv /path/to/download/directory/data/shared/
```

The `resources/shared/` tables are read **automatically** from `data/shared/`
on every run:

| File | Used by |
|------|---------|
| `disease_prevalence.tsv` | per-gene VCEP rules for BA1/BS1/BS2, PM2/PM4/PM5, PP2, PVS1 applicability, PS1 splice, BP1/BP3, BP7, REVEL cutoffs |
| `gene_inheritance.tsv` | inheritance-aware PM2 / BS1 / BS2 thresholds |
| `pm1_hotspots.tsv` | PM1 hotspot residues / regions |
| `pm4_regions.tsv` | per-gene PM4 region/strength rules (Strong residues, allow/deny regions, stop-loss strength, conservation / deletion-content gates, ABCA4 nucleotide-phyloP) |
| `ps1_paralog_map.tsv` | PS1 cross-gene paralogue residue map (SCN1A/2A/3A/8A, KCNQ1↔KCNQ2) |
| `vcep_pvs1_splice_exons.tsv` | optional per-skipped-exon PVS1 splice strengths (absent → flat per-gene splice defaults) |

The eRepo file is **not** auto-loaded — it is a ready-made manual-evidence
supplement (ClinGen Evidence Repository curated calls for criteria the tool
cannot derive automatically, e.g. PS4 / PP4). Pass it with `--supplement` at
classify/explain time — see [Manual evidence supplement](#manual-evidence-supplement).
It is assembly-specific, so use the `hg38` file for GRCh38 and the `hg19` file
for GRCh37.

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
- **All 28 ACMG criteria** supported, with ClinGen SVI strength adjustments and
  the Bergquist-2024 3-point extension. The ~18 derivable from local data
  (PVS1, PS1/PS3/PS4, PM1/PM2/PM4/PM5, PP1–PP3, BA1, BS1/BS2, BP1/BP3/BP4/BP7)
  are evaluated **automatically**; the evidence-dependent rest (PS2, PM3, PM6,
  PP4, BS3/BS4, BP2/BP5) are added via the [manual supplement](#manual-evidence-supplement),
  and PP5/BP6 are disabled per ClinGen SVI. See
  [Per-criterion decision basis](#per-criterion-decision-basis).
- **PVS1 decision tree** (Abou Tayoun et al. 2018) applied for LoF variants
  including NMD prediction, last-exon rescue, and biological-relevance gating.
- **Inheritance-aware PM2** (BS1/BS2 also) thresholds switch between dominant
  and recessive frequencies using a per-gene inheritance table.
- **Per-gene ClinGen VCEP rules** mined from the cspec specifications
  (`resources/shared/disease_prevalence.tsv`, `pm1_hotspots.tsv`): disease-
  specific BA1/BS1 cutoffs, PVS1 applicability (gain-of-function genes decline
  it) and the APC-specific PVS1 tree, PP2 applicability, PM5 Grantham-distance
  gating, PM1 hotspot regions, inheritance-aware BS2 (incl. a dominant
  heterozygote rule and a ≥3★ ClinVar fallback), the PS1 splice extension
  (canonical vs non-canonical), BP1/BP3 applicability, per-gene BP7 phyloP /
  intronic-range policy, and PM2 subpopulation / homozygote-count rules.
  X-linked "in males" genes are compared against the gnomAD male (XY) allele
  frequency; a homozygote/hemizygote-count BA1 rule (`ba1_hom_count`, e.g.
  SLC6A8/OTC ≥10) fires independently of frequency. See `resources/clingen/README.md`.
- **gnomAD frequencies** use the v4.1 **joint** (combined exome+genome) release
  on GRCh38; on GRCh37 (no joint release) the exome and genome callsets are both
  loaded and merged per-field. BA1/BS1/PM2 compare against the GrpMax filtering
  allele frequency (FAF95) by default; genes whose VCEP defines the cutoff on the
  point grpmax/popmax allele frequency instead use that (`af_basis=popmax`,
  globally revertible via `ACMG_POPMAX_AF_BASIS=false`). PM2 can additionally
  require a minimum gnomAD read depth (`pm2_min_depth`, ENIGMA BRCA1/2 ≥25) and
  judge absence on the non-cancer subset (`pm2_subset=non_cancer`). On GRCh38 the
  v4.1 release dropped the non-cancer subset, so it is backfilled from a small
  v3.1.2 genomes companion DB built **only for the BRCA1/2 contigs (chr13/chr17)**
  by default — the only genes that use it; GRCh37's v2.1.1 build carries the
  subset inline. The read-depth gate comes from the optional gnomAD coverage
  build (see [Data setup](#data-setup)).
- **In silico prediction**: **ESM1b** (default; Brandes 2023, MIT-licensed,
  commercial-use ready), **AlphaMissense** (non-commercial), or **REVEL**
  (Ioannidis 2016; free for non-commercial use, commercial use needs a separate
  licence) for missense. Splice prediction defaults to **OpenSpliceAI** (Chao
  2025, GPL-3.0; runtime inference via the `openspliceai` CLI); **SpliceAI**
  (Illumina-licensed) is an opt-in alternative. The splice predictor overrides
  the missense call — including on missense variants — when its score crosses
  the high-impact threshold (≥ 0.20). The missense predictors use Bergquist 2024
  Table 2 strengths; for REVEL a VCEP-specified gene cutoff (mined from cspec
  into `disease_prevalence.tsv`) overrides the genome-wide default when present.
  > _**Opt-in auxiliary predictors (BayesDel, CADD).** A few VCEPs define PP3/BP4
  > on BayesDel (ENIGMA BRCA1/2, TP53) or on an agreement of REVEL with CADD /
  > AlphaMissense (CTLA4, PIK3CD, PIK3R1, BMPR2, ABCA4). Both are academic /
  > non-commercial-licensed, so they are **licence-gated**: consulted only when
  > `--insilico-tool` is `revel` or `alphamissense` **and** the data is staged
  > (`--with-bayesdel` / `--with-cadd`); they are **never read under the
  > commercial-safe `esm1b` default**. When active, the per-gene VCEP rule is
  > authoritative (e.g. CTLA4 cannot meet PP3 on REVEL alone — REVEL∧CADD must
  > agree). TP53 uses the VCEP's published per-missense code table
  > (`resources/shared/tp53_pp3_bp4_codes.tsv`, which bakes in the Align-GVGD
  > class this tool does not compute); the assigned code plus its Align-GVGD and
  > BayesDel values are carried into the PP3/BP4 evidence._
  > _OpenSpliceAI reuses SpliceAI's 0–1 delta scale; lacking an OddsPath
  > calibration of its own, its PP3 is awarded at the conservative Supporting
  > tier (vs SpliceAI's Moderate). SQUIRLS and MMSplice are retained in the
  > code but disabled (SQUIRLS' precomputed DB is no longer downloadable;
  > MMSplice has a dependency conflict)._
- **Fully offline** at classification time (OpenSpliceAI runs locally too).
  Local databases include: Ensembl reference genome, VEP cache, gnomAD
  (DuckDB), ClinVar (VCF + derived SQLite for PS1/PM5), ESM1b / AlphaMissense,
  RepeatMasker, optional SpliceAI (opt-in, Illumina-licensed), and optional
  CADD / BayesDel (opt-in auxiliary PP3/BP4 predictors).
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
| Disk | ~ 1.5 TB free during setup | gnomAD VCFs dominate (v4.1 joint on GRCh38; exomes + genomes on GRCh37). The raw VCFs are deletable after the DuckDB build — see [Reclaiming disk space](#reclaiming-disk-space-after-setup) — leaving a much smaller steady-state footprint |
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
git clone https://github.com/hnakahara/HUVar.git
cd HUVar

# 2. Create & activate environment (conda example)
conda create -n acmg -c bioconda -c conda-forge \
    python=3.12 samtools tabix bcftools htslib ensembl-vep=111
conda activate acmg

# 3. Install the package (editable install for development).
# This pulls in OpenSpliceAI (the default splice predictor, GPL-3.0) and its
# bundled grch37/grch38 gene annotations as a regular dependency.
pip install -e .
```

---

## Data setup

A one-shot setup script downloads and builds everything required. The gnomAD
per-chromosome VCFs are fetched from the **Google Cloud and AWS mirrors in
parallel** (one concurrent download from each, byte-identical), and interrupted
downloads **resume** on re-run (size-verified against the remote, so a partial
file left by Ctrl+C is completed rather than skipped):

```bash
# Default: GRCh38 only, downloads everything (~ 1.5 TB, takes hours–days)
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh38 --workers 12

# Both assemblies
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh38 --workers 12
python scripts/setup_data.py --data-dir /path/to/download/directory/data --assembly GRCh37 --workers 12

# If you already have a reference genome / gnomAD VCFs locally, point to them
python scripts/setup_data.py --data-dir /path/to/download/directory/data \
    --genome-fasta /db/reference/GRCh38/hg38.fa \
    --gnomad-vcf-dir /db/gnomad/v4.1/joint/vcf \
    --workers 12

# Skip the gigantic gnomAD download (~1.5 TB)
python scripts/setup_data.py --data-dir /path/to/download/directory/data --skip-gnomad --workers 12

# Pick specific chromosomes for gnomAD (testing/partial setup)
python scripts/setup_data.py --data-dir /path/to/download/directory/data \
    --gnomad-chromosomes chr1 chr2 chrX

# The GRCh38 non-cancer companion (PM2 BRCA1/2) builds ONLY chr13+chr17 by
# default. Rebuild it genome-wide if a future VCEP adds a non-cancer gene:
python scripts/setup_data.py --data-dir /path/to/download/directory/data \
    --gnomad-noncancer-full --only gnomad-noncancer --workers 12

# Refresh ClinVar to the latest weekly release (re-download VCF + XML, rebuild SQLite)
python scripts/setup_data.py --data-dir /path/to/download/directory/data \
    --force-clinvar --only clinvar-vcf clinvar-sqlite --workers 12

# Opt-in auxiliary predictors — CADD (downloads the ~80 GB GRCh38 SNV file and
# normalises it) and BayesDel (staged from a manually-downloaded raw file):
python scripts/setup_data.py --data-dir /path/to/download/directory/data --only cadd --with-cadd
python scripts/setup_data.py --data-dir /path/to/download/directory/data --only bayesdel \
    --with-bayesdel --bayesdel-raw /path/to/BayesDel_noAF.tsv
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
| `--force-clinvar` | Force a fresh ClinVar download/rebuild even when local files exist (ClinVar is a rolling weekly release at a fixed URL). Re-acquires the VCF, the source RCV XML, and the PS1/PM5 SQLite. Combine with `--only clinvar-vcf clinvar-sqlite` to refresh ClinVar alone. |
| `--skip-gnomad` | Skip gnomAD download (~ 1.5 TB) |
| `--skip-gnomad-coverage` | Skip the gnomAD exomes coverage summary download + DuckDB build (per-locus mean read depth; used by the PM2 read-depth gate for ENIGMA BRCA1/2). Downloaded by default. |
| `--skip-gnomad-noncancer` | Skip the gnomAD v3.1.2 non-cancer companion DB (GRCh38). PM2 for ENIGMA BRCA1/2 then falls back to the overall AF. |
| `--gnomad-noncancer-chromosomes CHR ...` | Contigs to build the non-cancer companion DB from. **Default: only the BRCA1/2 contigs (`chr13 chr17`)** — the sole genes whose VCEP judges PM2 absence on the non-cancer subset, so the full v3.1.2 genomes download is avoided. GRCh37 ignores this (its v2.1.1 build carries the subset inline). |
| `--gnomad-noncancer-full` | Build the non-cancer companion DB genome-wide (all 24 contigs) instead of just `chr13`/`chr17`. Use when a future VCEP adds a non-cancer-subset gene outside the BRCA loci. |
| `--skip-genome` | Skip reference FASTA download (~ 880 MB) |
| `--skip-vep-cache` | Skip VEP cache download (~ 14 GB) |
| `--skip-esm1b` | Skip ESM1b LLR archive download / SQLite build (~ 1.34 GB) |
| `--skip-openspliceai` | Skip the OSAI_MANE model download (all 4 flanking sizes). Default tool, so skip only when supplying models another way. The `openspliceai` CLI itself ships as a package dependency. |
| `--with-revel` | Download REVEL (~ 600 MB zip) and build the per-assembly TSV for `--insilico-tool revel` (opt-in; ESM1b is the default tool) |
| `--with-cadd` | Build the per-assembly CADD TSV (opt-in auxiliary PP3/BP4 predictor for CTLA4/PIK3CD/PIK3R1/BMPR2/ABCA4). Downloads the CADD whole-genome SNV file (**GRCh38 v1.7 ~80 GB**) and normalises it. Consulted only under `--insilico-tool revel/alphamissense`. |
| `--cadd-raw PATH` | Use an already-downloaded CADD `whole_genome_SNVs.tsv.gz` (`Chrom Pos Ref Alt RawScore PHRED`) instead of downloading. |
| `--with-bayesdel` | Build the per-assembly BayesDel TSV (opt-in auxiliary PP3/BP4 predictor for ENIGMA BRCA1/2). Requires `--bayesdel-raw`. Consulted only under `--insilico-tool revel/alphamissense`. _(TP53 uses the shipped VCEP code table and only needs this flag toggled on at classification time — see below.)_ |
| `--bayesdel-raw PATH` | Path to a raw BayesDel score file (download the **no-AF** scores from <https://fenglab.chpc.utah.edu/BayesDel/BayesDel.html>). Layout: `chrom pos ref alt … score(last column)`. |

> **TP53 PP3/BP4** uses the ClinGen TP53 VCEP's published per-missense code table
> (`resources/shared/tp53_pp3_bp4_codes.tsv`, shipped with the repo and copied to
> `data/shared/` like the other curated tables — no download needed). It encodes
> the Align-GVGD class this tool does not compute. It is gated with the BayesDel
> family, so enable it at classification time with `--with-bayesdel` under
> `--insilico-tool revel/alphamissense`. Rebuild the table from a newer VCEP
> spreadsheet with `python scripts/build_tp53_codes.py --xlsx <file>`.

After setup, the expected layout is:

```
data/
├── shared/                          # curated gene-rule tables (copied from resources/)
│   ├── gene_inheritance.tsv         #   gene → AD/AR/XL
│   ├── disease_prevalence.tsv       #   per-gene VCEP rules (BA1/BS1, PP2, PM5, BS2, PS1, BP1/BP3, …)
│   ├── pm1_hotspots.tsv             #   per-gene PM1 hotspot residues/regions
│   ├── pm4_regions.tsv              #   per-gene PM4 region/strength rules
│   ├── ps1_paralog_map.tsv          #   PS1 cross-gene paralogue residue map
│   ├── vcep_pvs1_splice_exons.tsv   #   per-skipped-exon PVS1 splice strengths
│   └── tp53_pp3_bp4_codes.tsv       #   ClinGen TP53 VCEP per-missense PP3/BP4 codes (aGVGD+BayesDel)
├── vep_cache/                       # VEP indexed cache, both assemblies
├── esm1b/                           # (optional) protein-coordinate, shared across assemblies
│   └── esm1b_llr.sqlite             #   built from Brandes 2023 archive
└── GRCh38/                          # (mirror at GRCh37/)
    ├── genome/GRCh38.p14.fa(+.fai)
    ├── clinvar/clinvar_GRCh38.vcf.gz(+.tbi)
    ├── clinvar/clinvar_ps1_pm5_GRCh38.sqlite
    ├── gnomad/gnomad_v4.1_joint.duckdb   # GRCh37: gnomad_v2.1.1_exome_genome.duckdb
    ├── gnomad/vcf/                       # raw downloaded VCFs (~1.5 TB) — DELETABLE after build
    ├── gnomad/gnomad_v4.1_constraint.tsv
    ├── gnomad/gnomad_v4.0_exomes_coverage.duckdb   # per-locus mean DP (PM2 read-depth gate); GRCh37: v2.1
    ├── gnomad/gnomad_v3.1.2_non_cancer.duckdb      # GRCh38 only: PM2 BRCA1/2 non-cancer AF overlay (chr13/chr17 by default)
    ├── alphamissense/AlphaMissense_hg38.tsv.gz
    ├── (optional) revel/revel_grch38.tsv.gz(+.tbi)       # built with --with-revel
    ├── (optional) cadd/cadd_grch38.tsv.gz(+.tbi)         # built with --with-cadd (~80 GB source)
    ├── (optional) bayesdel/bayesdel_grch38.tsv.gz(+.tbi) # built with --with-bayesdel (--bayesdel-raw)
    ├── repeats/repeatmasker_dfam_hg38.bed.gz
    ├── (default splice) openspliceai/{80,400,2000,10000}nt/   # OSAI_MANE 5-model ensembles
    └── (optional)        spliceai/spliceai_scores.raw.{snv,indel}.hg38.vcf.gz
```

> **OpenSpliceAI is set up automatically.** The `openspliceai` CLI is a package
> dependency (installed by `pip install -e .`), and it bundles the grch37/grch38
> gene annotations the `-A` flag resolves — so no annotation file is downloaded.
> `setup_data.py` then downloads the OSAI_MANE 5-model ensemble for **all four
> flanking sizes** (80 / 400 / 2000 / 10000 nt) from the JHU CCB FTP into
> `data/<asm>/openspliceai/<flanking-size>nt/`. The classifier uses the
> `--openspliceai-flanking-size` model (default `2000`); the OpenSpliceAI authors
> recommend `10000` for best accuracy. Opt out of the model download with
> `--skip-openspliceai` (e.g. when supplying models via `--openspliceai-model-dir`).

### Reclaiming disk space after setup

Most of the ~1.5 TB is the **raw gnomAD VCFs**, which are only consumed while
building the local DuckDB. Once `setup_data.py` reports the gnomAD DuckDB as
built (`gnomad/gnomad_v4.1_joint.duckdb` exists), the raw VCFs are no longer read
at runtime and can be deleted:

```bash
# Safe to remove once the DuckDB build is complete — frees ~1.5 TB
rm -rf /path/to/download/directory/data/GRCh38/gnomad/vcf
# (and likewise data/GRCh37/gnomad/vcf if you built GRCh37)
```

If you later need to rebuild or add chromosomes, re-run `setup_data.py` (it will
re-download only what is missing; the existing DuckDB is reused as-is).

**ClinVar VCF** (`clinvar/clinvar_<assembly>.vcf.gz`) is *optional* at runtime:
PS1/PM5 read the separate `clinvar_ps1_pm5_<assembly>.sqlite`, and the VCF is
used only to add the public ClinVar classification to the output for human
review. You may delete it if you don't need that annotation column — classification
results are unaffected:

```bash
# Optional: removes the ClinVar classification column from output, not the calls
rm -f /path/to/download/directory/data/GRCh38/clinvar/clinvar_GRCh38.vcf.gz*
```

Keep everything else (DuckDB, ClinVar SQLite, VEP cache, ESM1b, OpenSpliceAI
models, constraint/AlphaMissense tables) — those are read on every run.
> If the model directory is absent, splice scoring is skipped.

> **Keeping ClinVar current.** ClinVar is published as a rolling weekly release at
> a fixed URL, so a normal re-run *skips* an existing ClinVar VCF/SQLite and they
> grow stale. Use `--force-clinvar` to re-download the VCF + source RCV XML and
> rebuild the PS1/PM5 SQLite from the latest release; pair it with
> `--only clinvar-vcf clinvar-sqlite` to refresh ClinVar without touching the
> large gnomAD/VEP/ESM1b data.

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
| `--insilico-tool {esm1b,alphamissense,revel}` | str | `esm1b` | Missense predictor used for PP3/BP4 (`revel` needs `setup_data.py --with-revel`) |
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

**Adding manual evidence.** `explain` accepts curator evidence the same way
`classify` does — inline via `--evidence` (repeatable, most convenient for a
single variant) and/or from a `--supplement` TSV — combined with the automated
calls per `--supplement-mode` (`merge` default, or `manual-only`):

```bash
# Inline: add a functional-study PS3 and a hotspot PM1
acmg-classify explain chr17 7674221 G A \
  --evidence PS3:strong \
  --evidence PM1:moderate:"PMID 12345 hotspot"

# From a supplement TSV (only rows matching this variant are applied)
acmg-classify explain chr17 7674221 G A --supplement manual_evidence.tsv
```

- `--evidence` format is `CRITERION:STRENGTH[:NOTE]`. Criterion (`PVS1`, `PS3`,
  `PM1`, `BA1`, …) and strength are **case-insensitive**; strength accepts the
  canonical values (`VeryStrong`, `Strong`, `ThreePoint`, `Moderate`,
  `Supporting`) or friendly aliases (`very_strong`, `strong`, `three_point`,
  `moderate`, `supporting`, plus `stand_alone` for BA1). The optional third
  field is a free-text note (it may itself contain `:`).
- See [Manual evidence supplement](#manual-evidence-supplement) for the TSV
  format and how `merge` / `manual-only` behave.

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

> 📦 **Ready-made eRepo supplement.** A curated supplement of ClinGen Evidence
> Repository manual criteria ships per assembly at
> `resources/<assembly>/erepo_manual_criteria_{hg38,hg19}.tsv` (copied to
> `data/shared/` in [Step 4](#quick-start)). It carries calls the tool cannot
> derive automatically (e.g. PS4, PP4) for variants in the eRepo. Use it as-is —
> match the file to your assembly:
> ```bash
> acmg-classify classify input.vcf -o results.tsv --assembly GRCh38 \
>   --data-dir /path/to/download/directory/data \
>   --supplement /path/to/download/directory/data/shared/erepo_manual_criteria_hg38.tsv
> ```
> Combine it with your own curation by concatenating rows into one TSV (or keep
> separate files and run in stages).

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

The same supplement TSV works with `acmg-classify explain` (single variant) via
`--supplement`, and for a single variant you can skip the file entirely and pass
evidence inline with one or more `--evidence CRITERION:STRENGTH[:NOTE]` flags —
see [`explain`](#explain--single-variant-detail). Both honour `--supplement-mode`.

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
`revel_score`,
`splice_tool`, `splice_score`,
`in_repeat`, `repeat_class`

`esm1b_llr` is populated when `--insilico-tool esm1b` (default);
`alphamissense_*` when `--insilico-tool alphamissense`; `revel_score` when REVEL
data is present. All three score columns are emitted for review whenever the
corresponding data file is available, but only the tool selected by
`--insilico-tool` feeds PP3/BP4. `splice_tool` / `splice_score` reflect the
active `--splice-tool` (`openspliceai` by default).

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

### Per-criterion decision basis

Each automatically-evaluated criterion, the data and thresholds it decides on,
and whether a ClinGen **VCEP gene-specific rule** can override the default.
Gene-specific parameters are stored per gene in
`resources/shared/disease_prevalence.tsv` (column names in `code` font) and the
`*_genes.py` / `pvs1/` modules; when a VCEP rule applies it is **authoritative**
and replaces the genome-wide default. The audit-trail `*_evidence` column records
which path fired.

Three cross-criterion rules are enforced after the per-criterion calls:
the allele-frequency criteria are **mutually exclusive** (BA1 > BS1 > PM2 — a
variant gets at most one); **PVS1 and PP3 are not co-applied**; and PVS1 strength
caps interact with ClinVar P/LP-null counts (see [PVS1](#pvs1-decision-tree)).

**Pathogenic**

| Criterion | Decision basis (default) | VCEP gene-specific rule |
|-----------|--------------------------|-------------------------|
| **PVS1** | LoF consequence (nonsense/frameshift/canonical-splice/start-loss/whole-gene-deletion) run through the Abou Tayoun 2018 tree (NMD, last-exon, rescue transcript, single-exon). See [PVS1 decision tree](#pvs1-decision-tree). | **Extensive.** `pvs1` applicability gate (LoF-not-the-mechanism genes withhold PVS1) + gene-specific trees for ~60 genes (`pvs1/vcep_pvs1.py`, APC in `pvs1/apc.py`): codon-range / critical-domain / NMD-boundary gates, per-gene initiation-codon and canonical-splice strengths, and an optional per-skipped-exon splice table. |
| **PS1** | Same amino-acid change as an established P/LP ClinVar variant (≥`pm5_min_stars`). Strong if any comparator is Pathogenic; Moderate if all are only Likely pathogenic. | `ps1_splice` extends PS1 to splice-equivalent variants for genes whose VCEP allows it. |
| **PS3** | Functional evidence: manual supplement (may reach Strong with OddsPath) or ClinVar SCV text-mining (1–2 SCV → Supporting, ≥3 → Moderate cap). | — (strength via supplement). |
| **PS4** | Unrelated affected-proband count from ClinVar P/LP `AffectedStatus=yes` SCVs: ≥10 → Strong, 6–9 → Moderate, 2–5 → Supporting. Gated on rarity (FAF95_popmax < 0.0001 or AC=0) **and** ≥2 probands. | — |
| **PM1** | Mutational hotspot / critical domain from the per-gene cspec table (`pm1_hotspots.tsv`); statistical fallback where no VCEP table exists. | **Yes** — VCEP hotspot residue ranges/residues + strength are authoritative. |
| **PM2** | Rarity on the **raw** gnomAD grpmax AF: dominant < 0.0001, recessive/X-linked < 0.005 → Supporting (SVI default). | **Yes** — `pm2_threshold` / `pm2_strength` / `pm2_basis`; `pm2_subpop` (`point` = also cap GrpMax point AF, `ci95` = upper-95%-CI rule); `pm2_zygosity` homo/hemizygote ceiling (e.g. SLC6A8 0, OTC ≤1, ABCD1 0 hemi); `pm2_subset=non_cancer` and `pm2_min_depth` (ENIGMA BRCA1/2). Gene-specific cSpec wording is hard-coded for a few genes (F8/F9 "absent in males", RYR1 "1 allele allowed", ATM "n=1 in a single subpopulation", PTEN single-vs-multi-allele subpop, RUNX1 GrpMax-FAF-then-all-subpop). |
| **PM4** | Protein-length change from an in-frame indel or stop-loss (outside repeat regions). | `pm4 = not_applicable` withholds PM4 for genes whose VCEP declined it. |
| **PM5** | Different missense at a codon with an established P/LP missense in ClinVar — Moderate (Supporting if comparators are LP-only). | **Yes** — `pm5_grantham` (require ≥ comparator Grantham), `pm5_excludes` (not co-applied with PM1/PS1 for some genes), `pm5_max`, `pm5_lp`, `pm5_min_count` (ACVRL1/ENG require ≥2 distinct same-codon LP/P → Strong). |
| **PP1** | Cosegregation: manual supplement or ClinVar text-mining; capped at Supporting (no meiosis counting from free text). | — |
| **PP2** | Missense in a gene where missense is a common mechanism & benign missense is rare. | **Yes (dominant lever)** — `pp2` applicability is authoritative; `pp2_requires` adds co-requirements (e.g. BMPR2 needs PM2+PP3). ClinVar-stats fallback only when no VCEP covers the gene. |
| **PP3** | Computational deleterious prediction (missense: ESM1b/AlphaMissense/REVEL; splice: OpenSpliceAI/SpliceAI), Bergquist 2024 tiers. See [In-silico aggregation](#in-silico-aggregation-pp3--bp4). | **Yes (REVEL)** — per-gene `revel_pp3_*` cutoffs from cspec override the genome-wide default and cap the gene's strength. **Opt-in auxiliary rules (`--with-bayesdel`/`--with-cadd`, licence-gated to REVEL/AlphaMissense):** BayesDel for ENIGMA BRCA1/2 (domain-gated) and TP53 (VCEP code table); REVEL∧CADD agreement for CTLA4/PIK3CD/PIK3R1; 2-of-3 REVEL/AM/CADD for BMPR2; CADD for ABCA4 synonymous/indel. When active the gene rule is authoritative. |

**Benign**

| Criterion | Decision basis (default) | VCEP gene-specific rule |
|-----------|--------------------------|-------------------------|
| **BA1** | gnomAD FAF95_popmax ≥ cutoff (stand-alone). Default 5%; disease-specific `min(0.05, 10×maxAF)` where parameters exist. | **Yes** — `ba1_threshold` / `af_basis` (`males` = XY AF, `popmax` = point grpmax AF) / `ba1_hom_count` (homo/hemizygote-count rule) per gene. |
| **BS1** | gnomAD FAF95_popmax above the disorder-specific expectation. | **Yes** — `bs1_threshold` / `bs1_strength`; `bs1_exclude` bars a specific recurrent disease allele from BS1 (e.g. MYOC p.Gln368Ter). |
| **BS2** | Observed in healthy adults in gnomAD, **inheritance-aware** (AR→homozygotes, XL→hemizygotes, AD→het carriers). | **Yes** — `bs2` applicability (VCEPs barring population data → withheld), `bs2_count` threshold, `bs2_female_only`, `bs2_hom_only`; a ≥3-star ClinVar BS2 assertion can substitute where the VCEP bars gnomAD. |
| **BP1** | Variant-type-vs-mechanism: applied only for genes whose VCEP names a target consequence. | **Yes (gate)** — `bp1` / `bp1_target` (`missense` for PALB2/APC/BRCA1/2; `truncating` for GoF RASopathy genes), `bp1_strength`, `bp1_exclude`, `bp1_no_splice`. No VCEP decision → not applied. |
| **BP3** | In-frame indel in a repetitive region of unknown function. | **Yes** — VCEP-gated (`bp3`) + `bp3_regions`. |
| **BP4** | Computational no-impact prediction (same tools as PP3), Bergquist 2024 tiers. | **Yes (REVEL)** — per-gene `revel_bp4_*` cutoffs. **Opt-in auxiliary rules** mirror PP3 (BayesDel for BRCA1/2 & TP53; REVEL∧CADD for CTLA4/PIK3CD/PIK3R1; 2-of-3 for BMPR2; CADD for ABCA4 synonymous/indel), licence-gated and authoritative when active. |
| **BP7** | Synonymous / deep-intronic with no predicted splice impact and (default) low conservation. Safe-distance ≥ +7 (donor) / ≤ −21 (acceptor). | **Yes** — `bp7_phylop` conservation cutoff (or `na` where the VCEP deemed conservation non-informative), `bp7_intronic = noncanonical` extends BP7 to any non-±1,2 intronic position. |

**Not auto-applied.** PS2, PM3, PM6, PP4, BS3, BS4, BP2, BP5 require evidence
(de novo confirmation, trans/cis phase, segregation meioses, phenotype
specificity) that cannot be derived from the local databases — add them via
[manual evidence](#manual-evidence-supplement). PP5 / BP6 (reputable-source) are
**deliberately disabled** per ClinGen SVI 2018 to avoid double-counting ClinVar.

### PVS1 decision tree

The Abou Tayoun 2018 LoF decision tree is implemented at
`src/acmg_classifier/pvs1/` and decides whether PVS1 should fire as Very
Strong, Strong, Moderate, Supporting, or be entirely suppressed based on:

- NMD predicted vs escape
- Last-exon / final 50 nt of penultimate exon rules — a truncation there with
  no functional-domain evidence is **N/A** (a critical region must be shown
  removed); a domain in the truncated tail downgrades to Strong
- Rescue transcript existence
- Region of biological relevance
- Single-exon gene logic

Three VCEP overlays sit in front of the generic tree (each authoritative where
it applies, evaluated before the generic tree and its strength caps):

- **Per-gene applicability gate** (`pvs1` column of `disease_prevalence.tsv`).
  Genes whose VCEP declares PVS1 not applicable because loss-of-function is not
  the disease mechanism — gain-of-function / dominant-negative disorders (MYOC,
  the RASopathy panel, the cardiomyopathy genes, the activating PIK3 genes,
  VWF, …) — withhold PVS1 entirely. Conversely, a VCEP that **explicitly applies**
  PVS1 (`pvs1=applicable`) has established LoF as the mechanism, so the decision
  tree skips the ClinVar/LOEUF LoF-mechanism heuristic and the undercuration
  strength caps for those genes (e.g. the Congenital Myopathies VCEP applies
  PVS1 Very Strong to ACTA1/RYR1 null variants — note RYR1's other VCEPs are
  gain-of-function, a multi-disease gene resolved here to the LoF context).
- **Gene-specific decision trees** for ~60 genes (`src/acmg_classifier/pvs1/vcep_pvs1.py`,
  plus APC in `pvs1/apc.py`). Each encodes its VCEP's deviations from the generic
  tree, in one of three forms:
  - **Codon-range bands** — e.g. VHL no-PVS1 before codon 54; CYP1B1 haem-binding
    domain through aa493; TP53 p.Lys351 boundary; GAA codon 916; the InSiGHT
    Lynch genes (MLH1 ≤753, MSH2 ≤891, MSH6 ≤1341, PMS2 ≤798); KCNQ1 1-581 /
    582-620 / 621-676; HBB's β-globin NMD window (codons 24-87). Some genes score
    nonsense vs frameshift differently (HNF1A p.601/p.618, AIPL1 p.328/p.337).
  - **Exon-based NMD** — NMD-predicted → Very Strong, NMD-escape → a fixed
    strength (FBN1, CDH1) or the generic 10%-of-protein rule (ACADVL, GAMT,
    ABCD1, RPGR).
  - **Initiation-codon / canonical-splice / whole-gene-deletion strengths** per
    VCEP (e.g. RPE65 start-loss Strong, GCK Supporting, VHL/DICER1 N/A; CDH1
    canonical-splice Strong default; ACADVL intron-8 GC-donor exclusion).
  - **APC** specifically: truncating variants are PVS1 only within NM_000038.6
    codons 49–2645, and canonical ±1,2 splice / "G→non-G last nucleotide" changes
    use the VCEP's allele-specific strength table (Lists A–E).
- **Optional per-skipped-exon splice override** (`resources/shared/vcep_pvs1_splice_exons.tsv`,
  consumed via `pvs1/vcep_pvs1_exons.py`). For canonical ±1,2 splice variants it
  refines the flat per-gene splice strength by the exon predicted to be skipped
  (in-frame / non-critical exon skips scored Strong or Moderate, e.g. DICER1,
  CDKL5, HNF4A, GAA). Ships with verified entries; exon numbering is generated
  from the MANE GFF by `scripts/build_vcep_pvs1_exons.py`.

The MANE-Select transcript and protein length each gene's bands are defined
against are recorded in `vcep_pvs1.py`.

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
- **REVEL** (`--insilico-tool revel`): Strong / ThreePoint / Moderate /
  Supporting for PP3 (`≥0.932 / ≥0.879 / ≥0.773 / ≥0.644`); Strong / ThreePoint
  / Moderate / Supporting for BP4 (`≤0.016 / ≤0.052 / ≤0.183 / ≤0.290`). No Very
  Strong category. When a ClinGen VCEP specifies its own REVEL cutoff for the
  gene (mined from cspec into the `revel_pp3_*` / `revel_bp4_*` columns of
  `disease_prevalence.tsv`), that gene-specific cutoff is used instead of the
  genome-wide default — and a VCEP that grants only a single strength caps the
  gene at that strength. REVEL scores are free for non-commercial use; see
  [Commercial use](#commercial-use). Requires `setup_data.py --with-revel`.

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
- **REVEL** (`--insilico-tool revel`) scores are **free for non-commercial use
  only**; commercial use requires a **separate licence** (contact the REVEL
  authors / Weiva Sieh — see <https://sites.google.com/site/revelgenomics/>).
  Not downloaded by default; opt in with `setup_data.py --with-revel`.
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
HUVar/
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
export ACMG_BS2_MIN_HOMALT=2       # BS2 healthy-homozygote count (recessive)
export ACMG_BS2_MIN_HEMI=2         # BS2 healthy-hemizygote count (X-linked)
export ACMG_BS2_MIN_HET=3          # BS2 healthy-carrier count (dominant)
export ACMG_POPMAX_AF_BASIS=true   # false → force every gene's BA1/BS1 onto FAF95
                                   # (ignore per-gene af_basis=popmax point estimate)
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
- **REVEL license.** The non-default `--insilico-tool revel` uses REVEL scores
  that are **free for non-commercial use only**; commercial use requires a
  **separate licence** from the REVEL authors (Weiva Sieh —
  <https://sites.google.com/site/revelgenomics/>). The default `esm1b` is
  MIT-licensed and commercial-use ready.
- **OpenSpliceAI is provisioned automatically.** The default `--splice-tool
  openspliceai` needs the `openspliceai` CLI (a package dependency, installed
  with `pip install -e .`, bundling the grch37/grch38 annotations) and OSAI_MANE
  model files (downloaded by `setup_data.py` for all four flanking sizes into
  `data/<asm>/openspliceai/<flank>nt/`; opt out with `--skip-openspliceai`). If
  the model dir is nonetheless absent at classification time, splice scoring is
  silently skipped and a warning logged.
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

HUVar (acmg-classifier). Department of Clinical and Molecular Genetics,
Hiroshima University Hospital, 2026.
https://github.com/hnakahara/HUVar
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

HUVar (`acmg-classifier`) is **not** an in vitro diagnostic device and is **not
certified for clinical decision making**. Any variant interpretation produced
by this tool must be reviewed by qualified clinical personnel before being
used to inform patient care.

The authors and the Department of Clinical and Molecular Genetics, Hiroshima
University Hospital make no warranty as to the accuracy or fitness for purpose
of the output, and accept no liability for clinical decisions made on the
basis of this software, in accordance with the "AS IS" provisions of the
Apache License 2.0.
