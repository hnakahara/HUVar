# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- BS2 inheritance-aware homozygote thresholds.
- SpliceAI fallback when AlphaMissense disagrees with splice prediction
  (SpliceAI ≥ 0.20 takes precedence).

### Changed

- GRCh37 transcript IDs in TSV output are now RefSeq (NM_) by default,
  matching GRCh38 behaviour.

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
