# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **`setup_data.py --force-clinvar`** — force a fresh ClinVar download/rebuild
  even when local files exist. ClinVar is a rolling weekly release at a fixed
  URL, so a normal re-run skips the existing VCF/SQLite and they go stale; the
  flag re-acquires the VCF, the source RCV XML, and the PS1/PM5 SQLite. Pair
  with `--only clinvar-vcf clinvar-sqlite` to refresh ClinVar alone.

- **PM4 region/strength engine** (`pm4_regions.tsv` + `PM4Regions`): per-gene
  PM4_Strong residues, allow/deny regions, region default, stop-loss strength,
  conservation (`conserved_phylop`) and deletion-content gates, mutual exclusion
  (PVS1/PP3), the `pm4_supporting_max_aa` size downgrade, and ABCA4
  nucleotide-conservation PM4 (`nt_phylop` ≥7.367, >1 nt → Moderate). VHL/RPGR
  region rules; MECP2/CDKL5/FOXG1/TCF4/UBE3A restricted to PM1 functional
  domains.
- **PM1 curated hotspots** for SCN1A/2A/3A/8A (VCEP "PM1 Table" Pathogenic
  Enriched Regions), RYR1 (pore Moderate + MH Supporting), FBN1 (cbEGF Cys /
  Gly-motif / Ca-binding consensus / Cys-creating), VHL, LDLR; ITGA2B/ITGB3
  marked not-applicable; PM1 restricted to missense except GAA/IDUA/CYP1B1.
- **PS1 cross-gene paralogue map** (`ps1_paralog_map.tsv`): SCN1A/2A/3A/8A and
  KCNQ1↔KCNQ2 analogous-residue PS1.
- **PM5 `pm5_min_count`** column — ACVRL1/ENG require ≥2 distinct same-codon
  LP/P comparators → PM5_Strong.
- **BA1 `ba1_hom_count`** column — homozygote/hemizygote-count BA1 (SLC6A8/OTC
  ≥10), independent of frequency.
- **Per-gene BA1/BS1 point-AF basis** (`af_basis=popmax`) for VCEPs that define
  the cutoff on the grpmax/popmax point allele frequency, with a global
  `ACMG_POPMAX_AF_BASIS` (default on) to revert all such genes to FAF95.
- **gnomAD exomes coverage** download + DuckDB build (default-on;
  `--skip-gnomad-coverage`) powering the PM2 read-depth gate (`pm2_min_depth`,
  ENIGMA BRCA1/2 ≥25) and non-cancer-subset PM2 (`pm2_subset=non_cancer`).
- **PM2 gene-specific cSpec wording** for F8/F9 ("absent in males"), RYR1
  ("1 allele allowed"), ATM ("n=1 in a single subpopulation"), PTEN
  (single-vs-multi-allele subpopulation), RUNX1 (GrpMax FAF then all-subpop).
- **PS3 per-gene control** — suppression for genes whose VCEP has no PS3 / only
  animal-model PS3 (PALB2/PDHA1/POLG/CAPN3/ANO5) and a Supporting cap for 15
  Supporting-only VCEPs; manual-supplement PS3 still takes precedence.
- **PVS1 per-gene VCEP applicability gate** (`pvs1` column in
  `disease_prevalence.tsv`). 33 genes whose VCEP declares PVS1 *not applicable*
  because loss-of-function is not the disease mechanism — gain-of-function /
  dominant-negative disorders (MYOC, the RASopathy panel, the cardiomyopathy
  genes, the activating PIK3 genes, RYR1, VWF, …) — now withhold PVS1 even on a
  bona-fide null variant.
- **APC-specific PVS1 decision tree** (Abou Tayoun 2018 / 2023 update,
  `src/acmg_classifier/pvs1/apc.py`): a codon-range gate (truncating variants
  are PVS1 only within NM_000038.6 codons 49–2645) and an allele-specific
  strength table (Lists A–E) for canonical ±1,2 splice and "G→non-G last
  nucleotide" exonic changes.
- **Gene-specific VCEP PVS1 decision trees** (`src/acmg_classifier/pvs1/vcep_pvs1.py`)
  for RPE65, CYP1B1, VHL, GCK, RAG1, ATM, GP9, IDUA, ACVRL1, PAH, HNF1A, GJB2,
  FOXG1, DICER1, PALB2, FBN1, GP1BA, CDH1, AIPL1, ACADVL, TP53, GAA, GAMT,
  HNF4A, RUNX1, CDKL5, RPGR, IL2RG, MECP2, F9 and ABCD1, fixing PVS1
  false-negatives where the generic ClinGen SVI tree under-calls.
  Each spec encodes the VCEP's critical-region / codon-range truncation gate
  (e.g. RAG1 / GP9 / GP1BA / FOXG1 / GJB2 single-exon or single-coding-exon genes
  where NMD is never predicted; CYP1B1 haem-binding domain through aa493; VHL
  no-PVS1 before codon 54; ACVRL1 NMD ≤442 / critical ≤490; IDUA NMD before
  c.1778; ATM through the most-3′ pathogenic p.R3047; PAH c.1285 boundary; DICER1
  NMD cutoff p.1850; PALB2 WD40 reaching the C-terminus; TP53 p.Lys351 boundary;
  GAA codon 916). Three handling modes are supported: fixed codon bands
  (with separate nonsense vs frameshift bands for HNF1A p.601/p.618 and AIPL1
  p.328/p.337), exon-based NMD with a fixed escape strength (FBN1, CDH1) and
  exon-based NMD with the generic 10%-of-protein escape rule (ACADVL, GAMT).
  Initiation-codon strength overrides (RPE65/PAH/ACADVL/GAA Strong, HNF1A/AIPL1/
  TP53 Very Strong, GCK/FOXG1 Supporting, GP1BA/GAMT Moderate, VHL/DICER1 N/A),
  the CDH1 canonical-splice Strong default, the ACADVL intron-8 GC-donor exclusion
  and whole-gene-deletion calls are also encoded. Like APC these run before the
  generic tree and its ClinVar-count strength caps.
- **Optional exon-aware VCEP PVS1 splice refinement**
  (`src/acmg_classifier/pvs1/vcep_pvs1_exons.py`,
  `scripts/build_vcep_pvs1_exons.py`). A reviewer-supplied TSV
  (`resources/shared/vcep_pvs1_splice_exons.tsv`) can override the flat per-gene
  canonical-splice strength per *skipped exon* (donor in intron *n* skips exon
  *n*; acceptor skips exon *n+1*), matching VCEP trees that score in-frame /
  non-critical exon skips at Strong or Moderate (e.g. DICER1, CDKL5). Opt-in: the
  file ships inert (only commented examples), so behaviour is unchanged until a
  reviewer activates rows. `build_vcep_pvs1_exons.py` generates a coordinate-
  accurate coding-exon table (in-frame flag, % protein, NMD geometry) from the
  MANE RefSeq GFF3 so exon numbering matches VEP and avoids the non-coding-exon-1
  off-by-one trap. Generated `resources/GRCh38/vcep_pvs1_exons.tsv` (MANE v1.5) and
  `resources/GRCh37/vcep_pvs1_exons.tsv` (RefSeq GCF_000001405.13, restricted to each
  gene's clinical transcript). Cross-checking the two assemblies validated the
  per-gene CDS lengths (28/30 exact matches) and surfaced four MANE-Select
  transcript/length corrections to `vcep_pvs1.py`:
  - **HNF4A** NM_000457.6 → **NM_175914.5** (452 aa; the MDEP c./p. numbering,
    incl. the p.419 PVS1 cutoff, matches this transcript — confirmed via the
    exon table).
  - **MECP2** NM_004992.4 (e2) → **NM_001110792.2** (e1, 498 aa); the e2 p.E472
    cutoff remapped to e1 p.484 (+12 shared-C-terminus offset).
  - **RPGR** NM_000328.3 → **NM_001034853.2** (ORF15 retinal isoform, 1152 aa);
    its terminal ORF15 exon is ~49% of the protein and glutamylation-critical, so
    every nonsense/frameshift now resolves to PVS1.
  - **GP1BA** 626 → **652 aa**, **RUNX1** 453 → **480 aa** (isoform c), plus
    ATM/GCK/PALB2 transcript-version label refreshes.
- **18 further gene-specific VCEP PVS1 trees** added after a systematic re-scan of
  the ClinGen cspec summary (`resources/clingen/cspec_json/cspec_summary.json`)
  for PVS1 rules with explicit codon / domain / NMD-escape detail: **ADA,
  DCLRE1C, JAK3, IL7R, FOXN1, RAG2** (SCID VCEP — NMD-escape 10% rule; IL7R TM
  domain; FOXN1 TAD bands; RAG2 single-exon core+PHD critical), **CTLA4** (exon-3
  codon 172/173/202 bands), **KCNQ1** (codons 1-581/582-620/621-676, SAD domain),
  **MLH1/MSH2/MSH6/PMS2** (InSiGHT PTC cutoffs 753/891/1341/798 with Moderate
  windows), **OTC** (C-terminal c.1033 Strong), **SLC9A6, TCF4, UBE3A** (codon
  bands; SLC9A6 +10 and UBE3A +20 MANE-isoform offsets), **GUCY2D** (LCA, p.1069
  Strong), **RS1** (X-linked retinoschisis, p.Met1-Cys223). The exon-table
  generator and both assembly tables were extended to cover these (GRCh38: 49
  genes; GRCh37 omits CDKL5/FOXN1/SLC9A6 whose transcripts post-date the 2013
  annotation). Remaining cspec PVS1 genes that only reference an external decision
  tree without codon detail (ABCA4, BMPR2, ENG, F8, GP1BB, ITGA2B, ITGB3,
  SERPINC1, PTEN, SCN1A/1B/2A/3A/8A, BRCA1/2, MYBPC3, PIK3R1, OTOF, MYO15A, GALT,
  HBB, HBA2, NEB) are left to the generic tree pending their attached flowcharts.
- **12 further gene-specific VCEP PVS1 trees** added from supplied decision-tree
  files: **ENG** (≤601 VS), **GP1BB** (TM domain, single-coding-exon like GP9/
  GP1BA), **SCN1B/SCN2A/SCN3A/SCN8A** (Epilepsy VCEP NMD-boundary codons
  204/1591/1586/1582 with 10% escape split), **NEB** (PTC ≤p.8452), **F8**
  (Hemophilia A, C-terminal C2 critical), **PTEN** (≤p.D375), **MYBPC3**
  (Cardiomyopathy, LoF-established, prior to p.1254), **HBB** (β-globin, NMD only
  in codons 24-87 per Peixeiro 2011 — early PTCs escape NMD but are Strong) and
  **HBA2** (α-globin, ≤p.Leu84). Exon-table generator extended to 61 genes
  (GRCh38). Still deferred (image-only files or no codon detail): ABCA4, GALT,
  BRCA1/2, BMPR2, ITGA2B, ITGB3, SERPINC1, PIK3R1, SCN1A, OTOF, MYO15A.
- **BS2 ClinVar expert-panel fallback.** For genes whose VCEP bars gnomAD
  population data for BS2 (CDH1, TP53, SERPINC1, …), a ≥3-star ClinVar review
  that explicitly applied BS2 is harvested (`bs2_evidence` / `bs2_strength`
  columns in the ClinVar DB) and used at the cited strength.
- **Per-gene BP7 conservation & intronic policy.** `bp7_phylop` — the phyloP
  "highly conserved" cutoff (0 / 0.1 / 0.2 / 1.5 / default 2.0) or `na` when the
  VCEP declared conservation non-informative (TP53, BRCA1/2, RUNX1, MYOC, the
  LCA genes, the SCID genes). `bp7_intronic = noncanonical` extends BP7 to any
  intronic position except the canonical ±1,2 (RASopathy / PIK3 panels, plus
  RUNX1 / MYOC / VHL).
- **PM2 highest-subpopulation correction** (`pm2_subpop`): `point` (RUNX1 —
  also require the GrpMax point AF ≤ threshold) and `ci95` (Cardiomyopathy/HCM
  and LGMD — require the upper 95% CI of the GrpMax AF ≤ threshold, reconstructed
  from new gnomAD `ac_grpmax` / `an_grpmax` columns). **PM2 homozygote/
  hemizygote ceiling** (`pm2_zygosity`): PM2 is withheld when gnomAD shows more
  homo-/hemizygotes than the VCEP tolerates (SLC6A8 0, OTC ≤1, the SCID genes /
  GATM / GAMT 0 homozygotes, ABCD1 0 hemizygotes).
- **Variant-level BS1 exclusion** (`bs1_exclude`): a recurrent disease allele
  the VCEP bars from BS1 regardless of frequency (MYOC p.Gln368Ter).
- **Automated OpenSpliceAI setup.** `openspliceai` is now a core dependency
  (installed with the package; it bundles the grch37/grch38 gene annotations the
  `-A` flag resolves, so no annotation file download is needed), and
  `scripts/setup_data.py` downloads the OSAI_MANE 5-model ensemble for **all four
  flanking sizes** (80 / 400 / 2000 / 10000 nt) from the JHU CCB FTP into
  `data/<asm>/openspliceai/<flank>nt/`. The model download runs by default
  (OpenSpliceAI is the default splice tool); opt out with `--skip-openspliceai`.
  Previously both the CLI install and model placement were manual steps.
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

- **PM2 no longer blanket-blocks on a failed gnomAD QC filter.** A filter-failed
  record (most often AC0 — zero high-quality observations) is effectively
  absent / extremely rare, so PM2 is now judged on rarity like any other record
  (a genuinely common filter-failed call still fails the threshold). Recovers a
  large class of false-negative PM2 calls on absent/very-rare null variants.
- **PVS1 last-exon (NMD-escape) refinement.** A truncating variant in the last
  (or penultimate) exon with no functional-domain evidence is now **N/A** (was
  Moderate): without evidence that a critical region is removed, the SVI tree
  does not apply PVS1. A domain in the truncated tail still yields Strong.
- **ITGA2B / ITGB3 PM2 threshold** now reads the operative gnomAD numeric tier
  (`<0.0001`) that sits behind the legacy "Absent from ESP/1000G/ExAC"
  boilerplate, instead of resolving to "must be absent" (threshold 0).
- **Empty Pilot/In-Prep cspec specs no longer shadow a populated spec.** A
  single-gene spec with zero criteriaCodes previously won the "most
  gene-specific" rule and blanked out a populated grouped spec's BA1/BS1; the
  populated spec now wins (recovers thresholds for 19 genes — AKT3/MTOR/PIK3CA/
  PIK3R2, the hearing-loss and mitochondrial panels, …).
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

- **eRepo benchmark false positives/negatives** addressed across criteria:
  - **PM1** — RYR1 pore-only over-firing (Malignant Hyperthermia broad regions
    were merged in), ITGA2B/ITGB3 ("highly polymorphic" → not applicable), and
    RUNX1 in-frame indels (PM1 is missense-only except GAA/IDUA/CYP1B1).
  - **PM4** — ABCA4 single-nucleotide false positive (now >1 nt only), RPGR
    Moderate region narrowed to exons 1–14 (ORF15 repeat excluded), and
    MECP2/CDKL5/FOXG1/TCF4/UBE3A in-frame indels restricted to PM1 domains.
  - **PS3** — text-mining no longer fires on "loss of function" prose, on
    "functional assays have not been reported / results not available", or on
    uncited quantitative assays (a PMID is required for quantitative claims);
    per-gene suppression / Supporting cap added.
  - **PM2** — gene-specific cSpec rules (F8/F9 "in males", RYR1 "1 allele
    allowed", ATM single-subpop allele, PTEN subpop tiers, RUNX1 FAF-priority)
    fix point-AF/FAF inflation false negatives.
  - **BA1→BS1 downgrades** — RYR1 BA1 threshold corrected (0.00697 → 0.0038,
    Malignant Hyperthermia VCEP), per-gene point-AF basis, and the
    homozygote-count rule.
  - **PVS1** — ACTA1/RYR1 null variants no longer blocked: the Congenital
    Myopathies VCEP applies PVS1 Very Strong (multi-spec resolution had taken a
    gain-of-function sibling spec). A VCEP that explicitly applies PVS1 now
    establishes LoF, bypassing the ClinVar/LOEUF heuristic and undercuration cap.
- **PS1/PM5 transcript-numbering collisions and frameshift mis-parsing.** Two
  ClinVar-builder bugs let PS1/PM5 match a spurious comparator: (1) a variant
  sharing only the `amino_acid_change` string on a *different* transcript's
  numbering (e.g. PTEN p.Pro38Leu on MANE vs the long isoform) — now rejected by
  a same-codon genomic-proximity guard (±2 bp); (2) a frameshift HGVS
  `p.Pro38LeufsTer*` parsed as the missense `p.Pro38Leu` (stored as `P38L`) —
  the residue regex now rejects a frameshift/extension tail.
- **RPGR BS1/BA1 cspec typo.** The cspec listed BS1 `≥8.3×10⁻⁵` and a legacy
  "5%" BA1 boilerplate; the VCEP's published classifications use BS1 `>5×10⁻⁶`
  and BA1 `>5×10⁻⁵` (BA1 = 10×BS1, male AF basis). Corrected via a curated
  threshold override.
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
