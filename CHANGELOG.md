# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **REVEL missense predictor** (Ioannidis et al. 2016) for PP3 / BP4 with
  Bergquist 2024 Table 2 strengths (PP3 `≥0.644 / 0.773 / 0.879 / 0.932`;
  BP4 `≤0.290 / 0.183 / 0.052 / 0.016`; no Very Strong tier). Select with
  `--insilico-tool revel`. When a ClinGen VCEP states a gene-specific REVEL
  cutoff, it overrides the genome-wide default — these are mined from cspec into
  new `revel_pp3_*` / `revel_bp4_*` columns of `disease_prevalence.tsv` (104 of
  131 released genes), and a VCEP granting only one strength caps the gene
  there. New `revel_score` column in the output TSV. `--with-revel` flag in
  `scripts/setup_data.py` downloads REVEL (~600 MB) and builds the per-assembly
  TSV (opt-in; ESM1b remains the default tool).
  **REVEL scores are free for non-commercial use only — commercial use requires
  a separate licence from the REVEL authors.**
- BS2 inheritance-aware homozygote thresholds.
- `--clinvar-workers N` flag in `scripts/setup_data.py` (default 4, max 24) to
  control ClinVar XML parse parallelism.
- **ESM1b missense predictor** (Brandes et al. 2023, MIT-licensed) for
  PP3 / BP4 with Bergquist 2024 Table 2 strengths. Enables commercial
  deployments where AlphaMissense's CC BY-NC-SA 4.0 licence is a blocker.
  Use `--insilico-tool esm1b`. New `esm1b_llr` column in the output TSV.
- `--skip-esm1b` flag in `scripts/setup_data.py` to opt out of the
  ~1.34 GB Brandes archive download.
- VEP runner now passes `--uniprot`, so `ConsequenceInfo.uniprot_id` is
  populated and the ESM1b lookup can key on the SwissProt accession.
- **X-linked male-AF for BA1/BS1.** New `af_basis` column in
  `disease_prevalence.tsv` (`males`) and gnomAD `af_xy` column / `GnomADData.af_xy`
  field. For genes whose VCEP states the cutoff "in males" (RPGR, RS1, ABCD1,
  SLC6A8, OTC) BA1/BS1 now compare against the gnomAD male (XY) allele
  frequency. Requires a gnomAD DB rebuild to populate `af_xy`; older DBs fall
  back to the overall FAF gracefully.
- `--override GENE:field=val[,...]` flag in
  `scripts/build_disease_thresholds.py` to manually pin a gene's BA1/BS1/
  `af_basis`/inheritance after multi-spec resolution (e.g. selecting the
  disease-appropriate RYR1 cutoff).
- **OpenSpliceAI splice predictor** (Chao et al. 2025, GPL-3.0) as the new
  **default** `--splice-tool`. Runs the OSAI_MANE model at inference time via
  the `openspliceai` CLI (`pip install openspliceai`); models live under
  `data/<asm>/openspliceai/<flank>nt/` (default `2000nt`) and are configurable
  via `--openspliceai-model-dir` / `--openspliceai-flanking-size`. Shares
  SpliceAI's 0–1 Δscale; PP3 is awarded at the conservative Supporting tier
  (no OddsPath calibration), BP4/BP7 at `≤0.10`, PVS1 splice-LoF at `≥0.20`.
  Temp-VCF headers include `##contig` lines so OpenSpliceAI's pysam writer
  accepts them.
- **`--supplement-mode {merge,manual-only}`** and supplement override extended
  to **all** criteria (previously only PS3/PP1/PP2 and the curation-only
  criteria honoured the supplement). `merge` (default): curator entries
  override the strength of any named criterion and add criteria the tool left
  not-met. `manual-only`: listed variants are classified purely from the
  supplement; variants not listed fall back to the tool's automated calls.
  Applied before the ACMG combination rules; the audit trail records overrides.
- **Dual-mirror parallel gnomAD download** in `scripts/setup_data.py`: per-
  chromosome VCFs are fetched from the Google Cloud and AWS mirrors in parallel
  (one concurrent download each), with automatic fallback to the other mirror
  on failure. Interrupted downloads now **resume** — the local file size is
  verified against the remote `Content-Length`, so a partial file from a
  Ctrl+C is completed (`wget -c`) instead of being mistaken for complete.

### Changed

- **Default missense predictor is now ESM1b** (was AlphaMissense). ESM1b is
  MIT-licensed, making the out-of-the-box configuration commercial-use ready.
  `--insilico-tool alphamissense` (CC BY-NC-SA 4.0) remains selectable.

- **PP2 now honours ClinGen VCEP applicability first.** New `pp2` column in
  `disease_prevalence.tsv` extracted from each VCEP's PP2 criteria code
  (`applicable` / `not_applicable`, with description-level blanket negations and
  per-gene exclusions like GN018 "but not PIK3R2"). The PP2 evaluator fires for
  `applicable` genes, suppresses `not_applicable` ones, and only falls back to
  the statistical heuristic for genes no VCEP covers — the dominant fix for PP2
  over-assignment (23 VCEP-applicable vs 108 not-applicable genes).
- **PP2 gene-specific co-requirements** (`pp2_requires` column). When a VCEP
  makes PP2 conditional on other criteria (BMPR2 / GN125: "PM2_supporting and
  PP3 must be met" → `PM2,PP3`), the registry suppresses PP2 post-hoc unless
  every required criterion is also triggered for the variant.
- **PP2 statistical fallback tightened** (used only for non-VCEP genes) to curb
  over-assignment: minimum P/LP missense `5 → 10`, maximum benign-missense
  fraction `0.10 → 0.05`, and the gnomAD missense-Z rescue branch now also
  requires a benign-missense rate ≤ 15% (a gene swamped with benign missense is
  no longer PP2-eligible on constraint alone). Thresholds are named constants in
  `clinvar_sqlite.py` for re-tuning.
- **Disease-specific BA1/BS1 thresholds rebuilt from ClinGen VCEP specs.**
  Multi-spec genes now prefer the more gene-specific VCEP (a single-gene panel
  supersedes a grouped panel); genuine ties across distinct diseases (RYR1,
  ACTA1) default to the most conservative — highest — cutoff to minimise
  false-positive benign calls.

- **ClinVar SQLite build is now parallelized** across a worker process pool
  (previously a single-threaded `iterparse`), substantially reducing build
  time on multi-core machines. Worker count defaults to 4 (max 24); override
  with `--clinvar-workers`.
- **Default splice predictor is now OpenSpliceAI** (was: splice evaluation
  disabled / `none`). SpliceAI remains opt-in (`--splice-tool spliceai`,
  Illumina-licensed). The active splice call overrides the missense call —
  including on missense variants — once its score ≥ 0.20.
- GRCh37 transcript IDs in TSV output are now RefSeq (NM_) by default,
  matching GRCh38 behaviour.

### Removed

- **MMSplice integration disabled** due to a dependency conflict. The code is
  retained (commented out) for future re-enablement.
- **SQUIRLS** is no longer CLI-selectable; its precomputed database is no
  longer downloadable. The predictor code is retained for when it returns.

### Fixed

- **AlphaMissense BP4 strength** previously returned `Strong` for
  `score ≤ 0.070`; per Bergquist 2024 Table 2 there is no Strong (-4)
  category for AlphaMissense BP4. Capped at `ThreePoint`. Existing
  Bayesian sums for very-benign AlphaMissense calls drop by 1 point.
- **PM5 fired on truncating comparators.** A pathogenic nonsense/frameshift
  variant that merely shared the residue number (`codon_position`) with the
  query satisfied PM5, which by definition requires a different *missense* at
  the same residue. The same-codon comparator is now restricted to genuine
  missense — removing a large class of false-positive PM5 calls.
- **BA1/BS1 threshold extraction** from VCEP free-text descriptions, which
  previously adopted the wrong cutoff in several specs:
  - KCNQ1 BA1 took the legacy "above 5%" (`0.05`) instead of the gnomAD-specific
    `≥0.004` — now `0.004`.
  - RPE65 / RUNX1 "between X and Y" BS1 bands took the upper bound — now the
    lower edge (RPE65 BS1 `0.008` → `0.0008`).
  - Rett/Angelman-like panels (CDKL5, FOXG1, MECP2, SLC9A6, TCF4, UBE3A) took a
    generic "above 0.05%" headline instead of the operative "≥0.000083 in any
    sub-population" cutoff — now `0.000083`.

### Known issues

- BA1 / BS1 / BS2 currently skip evaluation when the gnomAD FILTER is not
  PASS. This can leave clearly common variants un-classified as benign and
  inflate the Bayesian score. Tracked for a follow-up release.

## [0.1.0] — 2026-05-27

Initial public release.

### Added

- ACMG 2015 Table 5 combinatorial classifier (`classifier_2015.py`).
- Tavtigian 2020 / Bergquist 2024 Bayesian point-based classifier
  (`classifier_bayesian.py`).
- All 28 ACMG criteria evaluators (PVS1, PS1–PS4, PM1–PM6, PP1–PP5,
  BA1, BS1–BS4, BP1–BP7).
- PVS1 decision tree (Abou Tayoun 2018) with NMD prediction, last-exon
  rescue, and biological-relevance gating.
- Fully local annotation pipeline (Ensembl VEP, gnomAD DuckDB, ClinVar VCF
  + SQLite, AlphaMissense, RepeatMasker, optional SpliceAI).
- CLI: `classify`, `explain`, `validate`, `status`, `setup`.
- Manual evidence supplement (TSV override for PS3/PP1/PM3/etc.).
- GRCh37 and GRCh38 dual-assembly support.
- `scripts/setup_data.py` for one-shot database download and build.
