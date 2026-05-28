# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- BS2 inheritance-aware homozygote thresholds.
- SpliceAI fallback when AlphaMissense disagrees with splice prediction
  (SpliceAI ≥ 0.20 takes precedence).
- **ESM1b missense predictor** (Brandes et al. 2023, MIT-licensed) for
  PP3 / BP4 with Bergquist 2024 Table 2 strengths. Enables commercial
  deployments where AlphaMissense's CC BY-NC-SA 4.0 licence is a blocker.
  Use `--insilico-tool esm1b`. New `esm1b_llr` column in the output TSV.
- `--skip-esm1b` flag in `scripts/setup_data.py` to opt out of the
  ~1.34 GB Brandes archive download.
- VEP runner now passes `--uniprot`, so `ConsequenceInfo.uniprot_id` is
  populated and the ESM1b lookup can key on the SwissProt accession.

### Changed

- GRCh37 transcript IDs in TSV output are now RefSeq (NM_) by default,
  matching GRCh38 behaviour.

### Fixed

- **AlphaMissense BP4 strength** previously returned `Strong` for
  `score ≤ 0.070`; per Bergquist 2024 Table 2 there is no Strong (-4)
  category for AlphaMissense BP4. Capped at `ThreePoint`. Existing
  Bayesian sums for very-benign AlphaMissense calls drop by 1 point.

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
