# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **Somatic oncogenicity (experimental)**: `acmg-classify somatic` / `somatic-explain`
  implement the ClinGen/CGC/VICC oncogenicity SOP (Horak et al. 2022) with
  licence-free data only — NCG 7.1 canonical gene roles plus curated overrides
  (`somatic_gene_roles.tsv`), cancerhotspots.org v2 counts (`cancer_hotspots_v2.tsv`,
  ODbL), CIViC functional evidence / OS2-coded oncogenic assertions for OS2/SBS2
  (CC0, new `setup_data.py` step `civic`), ClinVar somatic-oncogenicity comparators with a different nucleotide change (OS1: ≥ 1 star; OM4: ≥ 2 stars; germline P/LP not used; new ClinVar SQLite columns `oncogenicity` / `onc_review_status` / `onc_stars` — rebuild with `--only clinvar-sqlite --force-clinvar`), gnomAD (OP4/SBS1/SBVS1), HUVar
  PP3/BP4/BP7/PM1 and calibrated functional data (OS2/SBS2). OP2 and any override via
  manual evidence. Builder `scripts/build_somatic_resources.py`. AMP/ASCO/CAP tiers
  are not assigned.
- **PS3/BS3 from calibrated functional assays** (`criteria/functional_assays.py`):
  ENIGMA BRCA1/2 Specifications Table 9 V1.2 (variant-level codes) and MaveDB score
  sets with clinical calibrations (dump 2026-06-24, CC0; primary calibration by
  default, `ACMG_MAVEDB_CALIBRATION=conservative`), strength from the range OddsPath
  per Brnich et al. 2019, per-gene VCEP rules in `functional_evidence_rules.tsv`
  (VCEPs with assay lists / decision trees excluded, e.g. GCK). TP53 uses MaveDB when
  the NCI TP53 Database flowchart is off or lacks the variant. BS3 is now automated
  for these sources. New `setup_data.py` step `mavedb`; builders
  `scripts/build_mavedb_functional.py`, `scripts/build_enigma_table9.py`.
  `ACMG_FUNCTIONAL_ASSAYS=false` restores the previous behaviour.
- **MANE Plus Clinical (reporting)**: VEP now runs with `--mane`, and each
  consequence carries `is_mane_plus_clinical` (recovered by accession from
  `mane_exons.tsv` on GRCh37). The TSV gains `mane_plus_clinical_transcript`,
  `_hgvs_c`, `_hgvs_p`, `_consequence` and `_consequence_differs` after `hgvs_p`,
  the JSON a `mane_plus_clinical` list, and the console report a Plus Clinical
  line. Criteria are still evaluated on MANE Select only (phase 1); evaluation
  on the Plus Clinical transcript is planned separately because several VCEP
  tables (PM1/PM4 hotspots, DYSF comparator rules, REVEL row choice) are keyed
  to MANE Select numbering.
- **MANE Plus Clinical (evaluation)**: each variant in a gene with a Plus Clinical
  transcript is also classified with the criteria evaluated on that transcript
  (`AnnotationData.for_transcript`); the TSV gains
  `mane_plus_clinical_{classification_2015,rules_2015,bayesian_score,classification_bayesian,criteria_changes}`
  and `explain` a `--transcript` option. Consequence, the generic PVS1 tree and the
  isoform ESM1b/REVEL scores (looked up by the orchestrator) follow the evaluated
  transcript; VCEP position rules and ClinVar codon comparisons keep MANE Select
  numbering through `AnnotationData.rule_consequence` (PM1/PS1/PM5 via
  `rule_view()`; position parts of BA1/BS1/BP1/BP3/BP4/PM4/PP3 and the
  gene-specific PVS1 trees). Off with `ACMG_EVALUATE_PLUS_CLINICAL=false`. The
  MANE Select classification is unchanged.
- `setup_data.py` step `mane-summary` (`--force-mane-summary`): downloads the NCBI
  MANE v1.5 summary into `data/shared/` (assembly-independent), validated for the
  `symbol` / `RefSeq_nuc` / `MANE_status` columns and recorded in the data manifest.
  HUVar_app's TransVar service reads it from there.

### Fixed

- **ClinVar (and other tabix) lookups dropped records next to non-ASCII lines.**
  Recent ClinVar releases contain UTF-8 disease names (e.g. "Muir-Torré_syndrome");
  pysam decoded lines as ASCII and the `UnicodeDecodeError` was swallowed by the
  chr-prefix fallback, so every record at a position overlapped by such a line was
  lost (e.g. MLH1 c.2059C>T, ClinVar 90014). Tabix files are now read as UTF-8, and
  the contig is resolved from the index instead of by catching `ValueError`.

### Changed

- **PM2 for the InSiGHT MMR genes (MLH1, MSH2, MSH6, PMS2) uses the overall gnomAD
  allele frequency** (`pm2_basis=global`), as worded in GN115/GN137/GN138/GN139
  ("allele frequency <0.00002 ... in gnomAD v4 dataset"), instead of the GrpMax
  point AF. Agreement with the eRepo PM2 calls for these genes rises from 34/37 to
  36/37 with no additional PM2 assignments. `build_disease_thresholds.py` emits
  `global` for this wording.

### Added (Genome Medicine revision)

- `benchmark/`: scripts used for the analyses in the manuscript (eRepo benchmark,
  agreement with the eRepo classifications, rule calibrations, resource measurements,
  figures), and the benchmark input VCFs `inputs/erepo_hg38_may_2026.vcf` / `inputs/erepo_hg19_may_2026.vcf`.

### Changed (Genome Medicine revision)

- **Expert-panel criteria only with curated evidence**: PS3, PS4, PP1 and BS2 applied
  by a ClinGen expert panel to the variant itself in ClinVar (>=3 stars) are curated
  evidence and are now imported only when curated evidence is in use for the run
  (`classify --supplement`, `explain --evidence`/`--supplement`; new
  `Config.curated_evidence`, set automatically and recorded in the provenance).
  Without curated evidence the classification is purely automated (previously the
  expert-panel record was read unless `--exclude-self-expert-panel` was given).
  The benchmark and CGP results are unchanged: their automated-only runs already
  used `--exclude-self-expert-panel`, and all runs with curated evidence used
  `--supplement`.
- **PVS1/PP3, splice annotations that are not canonical**: VEP calls a variant
  splice_donor/acceptor_variant from the left-aligned VCF allele, so a deletion in
  a repeat can be annotated as canonical although its HGVS (3') change lies beyond
  +/-2 (e.g. LZTR1 c.651+10_651+46del, APC c.531+5_531+8del). Such purely intronic
  variants no longer receive PVS1; their splicing effect is evaluated under PP3.
- **PVS1, canonical splice sites**: for canonical +/-1,2 splice donor/acceptor
  variants, splicing disruption is now presumed, as in the ClinGen SVI PVS1
  decision tree, and the splice-predictor score is reported for reference only.
  Previously a score below 0.20 down-graded PVS1 to Strong/Moderate; with the
  default 80-nt OpenSpliceAI model this affected about a third of canonical sites
  in the eRepo benchmark (e.g. RB1 c.381-1G>A, OpenSpliceAI 0.06 -> PVS1_Moderate;
  SpliceAI 1.00). VCEP-specific PVS1 splice rules still take precedence.
- **PVS1, PALB2 C-terminal truncations**: the PALB2 VCEP boundary is now judged by
  the position of the (new) stop codon rather than the first residue a frameshift
  alters. A PTC upstream of His1184 remains PVS1; a PTC after Tyr1183, including
  frameshifts reading through the native stop (C-terminal extension), is
  PVS1_Strong, as in the VCEP pilot c.3543del (p.Phe1181LeufsTer10). In the eRepo
  benchmark this changes c.3543del and c.3510del (both VCEP PVS1_Strong, VUS).
- **PVS1, alternative-transcript rescue**: only MANE/canonical protein-coding
  transcripts of the gene under evaluation are considered. VEP also reports
  transcripts of neighbouring genes within 5 kb (upstream/downstream_gene_variant,
  e.g. RB1-DT and LPAR6 for RB1), which had been counted as rescuing transcripts
  and down-graded PVS1 from Very Strong to Strong for NMD-predicted null variants
  (85 eRepo variants, e.g. ITGB3, CAPN3, LDLR; RB1, STK11, MEN1 in the CGP cohort).
- **BA1/BS1 without FAF95**: when gnomAD provides no FAF95 (a single observation in
  the grpmax group), BA1 and BS1 are no longer established from the point popmax AF
  (unless the VCEP defines the cutoff on the point AF or the male AF). The VCEP
  frequency rules are stated on FAF95 in populations with more than one
  observation; the fallback had produced 11 BA1/BS1 calls in the eRepo benchmark,
  none applied by the VCEP.
- **BA1/BS1 exclusion variants** for the Congenital Myopathies specifications:
  RYR1 AD (c.14524G>A, c.325C>T), RYR1 AR (c.6721C>T, c.10348-6C>G) and ACTA1 AR
  BA1 (c.541del, c.121C>T; also in the merged ACTA1 row).
- **PM2, Malignant Hyperthermia Susceptibility specification (RYR1, `--cspec
  RYR1=malignant_hyperthermia`)**: PM2 is no longer applied ("Not Applicable: PM2 is
  not used alone", GN012; the VCEP applied PM2 to none of its eRepo variants).
  `pm2_threshold = not_applicable` is now supported in `disease_prevalence*.tsv`.
- **REVEL, MANE row selection**: REVEL v1.3 lists one row per transcript-specific
  amino-acid change, and the same SNV can score very differently across
  transcripts (IDUA c.530T>G: F>C 0.961 on MANE vs S>A 0.400; FBN1 c.7016G>A:
  C>Y 0.973 vs A>T 0.209). `scripts/setup_data.py` now keeps aaref, aaalt and the
  Ensembl transcript IDs (8-column REVEL file), and `revel_db` picks the row of the
  MANE Select ENST, else the row matching the MANE amino-acid change, else the
  highest score. Previously the transcript columns were dropped and the first
  matching row was used. **The REVEL file must be rebuilt** (`setup_data.py
  --with-revel`; an existing 5-column file is renamed `*.5col.tsv.gz` and rebuilt);
  a legacy file still works with the old first-row behaviour.
- **PP3 / BP4, VCEP scope** (`pp3_scope`, new column in `disease_prevalence.tsv`;
  `bp4_scope` extended): the SCID VCEP genes (ADA, DCLRE1C, IL2RG, IL7R, JAK3, RAG1,
  RAG2) no longer receive PP3 for missense variants ("Do not apply to missense
  variants"; PP3 only for synonymous/intronic splice impact) and BP4 is not applied
  for ADA and RAG2 (as already for the other SCID genes); APC, PALB2 and CDH1 no
  longer use a protein-level predictor for missense PP3 (splice prediction only);
  the PIK3/mTOR GOF panel (AKT3, MTOR, PIK3CA, PIK3R2) no longer receives PP3, and
  APC and that panel no longer receive missense BP4. Previously these genes fell
  back to the genome-wide REVEL/AlphaMissense calibration.
- **PVS1, generic SVI tree — NMD escape**: when the truncated/altered region has no
  annotated functional domain ("role of region unknown"), the SVI protein-length
  branch now applies: >10% of the protein altered (counted from the first changed
  residue) → PVS1_Strong, ≤10% → PVS1_Moderate. Previously PVS1 was withheld.
- **PVS1, VCEP-specific rules**: ENIGMA BRCA1/BRCA2 Specifications Table 4 for PTC
  variants (by the PTC codon: BRCA1 exons 8(9)–9(10) and PTC > p.L1854 N/A; BRCA2
  PTC > p.E3309 N/A; otherwise PVS1); ITGA2B transmembrane domain (p.994–1019)
  as the critical region for NMD-escaping truncations (Strong; beyond it the 10%
  rule); GCK stop-loss variants → PVS1 (exon 10 critical, stop-loss reported in
  MODY); frameshifts reading through the stop codon → PVS1 for MECP2, TCF4, CDKL5
  and UBE3A, PVS1_Moderate for SLC9A6.
- **PP3, splice predictions** (Walker 2023): SpliceAI / OpenSpliceAI >= cutoff now
  gives PP3_Supporting (was Moderate for SpliceAI), with VCEP gene-specific cutoffs
  (>= 0.5: ABCD1, ACADVL, ACTA1, ANO5, CAPN3, DNM2, DYSF, MTM1, NEB, RYR1, SGCA/B/D/G,
  SERPINC1, VWF, GP1BA/GP1BB/GP9; > 0.5: GAA, PAH; >= 0.38: RUNX1; > 0.3: HBB, HBA2).
  For missense variants PP3 is the stronger of the protein-level predictor and the
  splice prediction (the splice route no longer pre-empts a stronger REVEL /
  ESM1b / AlphaMissense call).
- **PM2, Monogenic Diabetes VCEP (GCK / HNF1A / HNF4A)**: when gnomAD cannot compute
  the GrpMax FAF (too few copies), PM2_Supporting follows the MDEP copy rule (<= 2
  copies in NFE, <= 1 copy in any other group) instead of the point AF.
- **PM2, Hearing Loss VCEP**: autosomal recessive genes (USH2A, GJB2, SLC26A4, CDH23,
  MYO7A) use the recessive cutoff 0.00007 (the multi-gene spec has no per-gene MOI,
  so the dominant 0.00002 was used).
- **PM5, VCEP conditions** (`criteria/pm5_conditions.py`): splicing must be ruled out
  for the candidate and/or the comparator (gene-specific SpliceAI cutoffs; GUCY2D
  same score band; RPE65 within 0.1; comparator scores from the precomputed SpliceAI
  files), PP3 required (MLH1, MSH2, MSH6, PMS2, RPGR, RS1), REVEL > 0.7 (DYSF), two
  likely-pathogenic comparators = Moderate (F8, F9), MDEP Grantham downgrade to
  PM5_Supporting (GCK, HNF1A, HNF4A), and the SCN1A/2A/3A/8A paralogue route
  (PM5_Supporting for > 1 non-identical P/LP change at the analogous residue).
- **PM5, CDH1 canonical splice variants**: PM5_Supporting when PVS1 is met at Strong
  or higher (HDGC VCEP v3.1).
- **PM4, MYOC** (Glaucoma VCEP v2.1): truncating variants in the olfactomedin domain
  (p.246-502) get PM4 by the share of the protein involved (> 10% Moderate,
  <= 10% Supporting); small in-frame changes there are PM4_Supporting.
- **PM1, HNF1A promoter**: PM1_Supporting for the transcription-factor binding sites
  c.-187 to -195, c.-209 to -227, c.-238 to -259 and c.-276 to -288 (MDEP v3.1).
- **PVS1, VCEP "G to non-G at the last exon nucleotide" rules** (`pvs1/last_nt.py`),
  applied to exonic changes VEP calls missense / synonymous: ATM (HBOP v1.6; intron
  +1..+6 not `gtrrgt`) → the canonical donor weight one level down (PVS1_Strong;
  PVS1_Moderate for in-frame HEAT-repeat exons); MLH1 (InSiGHT; intron not GTRRGT)
  → PVS1_Strong; CDH1 (v3.1) → PVS1_Moderate. The intronic sequence is read from
  the reference FASTA (`genome_fasta`); without it the rule applies and the
  evidence notes that the intron was not checked.
- **PVS1, RUNX1 canonical splice sites**: strength per intron donor/acceptor from the
  Myeloid Malignancy VCEP "Summary of splicing effects" table (c.58/c.59/c.97/c.98
  isoform-c-only sites → N/A; c.351/c.352/c.508 → PVS1; c.509, c.613, c.614, c.805,
  c.806, c.967, c.968 → PVS1_Strong). Previously every canonical site was PVS1.
- **PS1 for splicing = ClinGen SVI Walker et al. 2023 Table 2**: the strength now
  depends on the variant position and its baseline code — outside ±1,2 (needs PP3):
  same-nucleotide comparator P → PS1, LP → PS1_Moderate; at ±1,2 with PVS1: comparator
  in the same ±1,2 dinucleotide P → PS1_Supporting, LP → not applicable; with
  PVS1_Strong/Moderate/Supporting: P → PS1, LP → not applicable; canonical without PVS1
  → not applicable (`registry.apply_ps1_splice_table2`). Previously every same-position
  comparator gave PS1 (P) / PS1_Moderate (LP) regardless of PVS1/PP3. Comparators at the
  other base of the same ±1,2 dinucleotide are now used
  (`clinvar_sqlite.query_splice_dinucleotide_partner`); comparators elsewhere in the
  splice motif are not. Synonymous variants are evaluated for splice PS1 with
  comparators that have no protein change.
- **PS1, BMPR2 splice-motif rule** (GN125 v2.0.0; `criteria/ps1_splice_motif.py`): a
  variant outside ±1,2 in the same donor motif (last 3 exonic bases, +3..+6) or acceptor
  motif (first exonic base, −20..−3) as a ClinVar P / LP variant → PS1_Moderate /
  PS1_Supporting (no PP3 requirement; exonic variants and comparators only when
  synonymous / without protein change; motif positions from the MANE exon structure;
  `clinvar_sqlite.query_plp_window`).
- **PS1 for splicing, VCEP-specific rules take precedence over Walker Table 2**: 29
  genes whose VCEP specifies its own splice PS1 rule (ABCD1; SCID: ADA, DCLRE1C, IL7R,
  JAK3, RAG1, RAG2; Hearing Loss: CDH23, COCH, GJB2, KCNQ4, MYO15A, MYO6, MYO7A, OTOF,
  SLC26A4, TECTA, USH2A; ENIGMA BRCA1/2; HBOP ATM/PALB2; InSiGHT MMR; APC, PTEN, DICER1)
  keep the same-nucleotide comparator at the comparator's strength (P → PS1, LP →
  PS1_Moderate); Table 2 and the ±1,2-dinucleotide comparator apply to the 36 genes
  whose specification cites Walker et al. 2023 (`ps1.PS1_SPLICE_VCEP_RULE`).
- **PS1 splice extension removed for missense-only specifications**: FBN1, LDLR, PIK3CD
  and KCNQ1 mention splicing only as a caveat (`ps1_splice` blank;
  `build_disease_thresholds._PS1_CAVEAT_SENT`).
- **VCEP-classified PS1/PM5 comparators** (`vcep_comparator_rules.tsv`, builder patterns
  extended): LDLR, DICER1 (Pathogenic comparators only — new `pathogenic_only` column),
  RPE65, GUCY2D, AIPL1, RPGR, RS1 (PM5: comparator assessed independently of PM5), KCNQ1
  PS1 (comparator classified without PS1), GP1BA PS1, GP1BB/GP9 PM5. GATM is not
  restricted ("determine its pathogenicity using these specifications").

- **eRepo manual-evidence supplement**: 8 more RNA-based PVS1 calls whose eRepo
  summaries cite RT-PCR / exon-trapping data without a "PVS1 (RNA)" label (GAA ×6,
  ITGA2B c.2828C>T, RPGR c.1245+3A>G); 157 PVS1 (RNA) rows per assembly.
- **PM1 regions re-checked against every VCEP CSpec text** (`pm1_hotspots.tsv`,
  `scripts/build_pm1_hotspots.py` `_CURATED`). Fixed: MECP2 (VCEP numbering is
  NM_004992 e2; MANE NM_001110792 e1 is +12 → MBD 102-174, TRD 314-318), UBE3A
  (catalytic Cys820 of the 852-aa isoform = 840 on MANE NM_130839), F9 Moderate
  (exons 3-5, p.85-174, was missing), BRAF (exon 6 p.238-287 and exon 11
  p.439-478 added), RAF1 (exon 14 p.473-512 and exon 17 p.602-648 added),
  DICER1 Supporting (RNase IIIb range 1682-1846, was two residues), GUCY2D
  Supporting (p.873-951, was missing), ACADVL (VCEP table ranges 1-40, 214-223,
  249-251, 381-382, 481-516 and p.562 added), RPE65 (Ala107-Gly125 range, was two
  residues), ENG (spurious 350-382 range removed; C350/C363/C382 as residues),
  HNF4A Supporting (37-113; 33-36 removed), TNNT2 (89-189 on MANE only), GP1BB
  (Cys147 added), PTPN11 (AA 4 added), F8 Moderate (Ser2283 excluded), RYR1 (MHS
  VCEP N-terminal 1-552 / central 2101-2458 Moderate, Supporting with PS1/PM5;
  C-terminal 4631-4991 Supporting; pore 4800-4950 Moderate).
- **PM1, VCEP external tables**: TP53 PM1 from cancerhotspots.org per amino-acid
  change (≥10 samples → Moderate, 2-9 → Supporting; export in
  `resources/clingen/cancerhotspots_TP53.tsv`, one `alt_aa`-restricted row per
  residue) in addition to the six VCEP codons. VHL from the VCEP "Germline and
  Somatic Hotspots" table: germline hotspots (Stebbins; Chiorean 2022; the
  struck-through somatic residues S65/N78/S80/P86/W88/Y112/W117/R161/C162 used as
  germline, W88 added) Moderate for any missense; the remaining somatic hotspots
  (S68, V74, L89, S111, G114, H115, D121, L135, I151, L158, L169) by
  cancerhotspots.org count of the same amino-acid change
  (`resources/clingen/cancerhotspots_VHL.tsv`): ≥10 → Moderate (none reaches
  10), observed <10 → Supporting, changes not observed → no PM1. SCN1A/2A/3A/8A (PM1 Table v2.20241205) and
  AKT3/MTOR/PIK3CA/PIK3R2 (Brain Malformations Table 4) verified unchanged.
- **PM1 VCEP caveats**: LDLR and KCNQ1 PM1 require PM2; RAG1/RAG2/IL7R PM1 not
  with BA1/BS1/BS2; PAH PM1 not with PP3_Strong; RYR1 Moderate→Supporting with
  PS1/PM5 (`PM1_DOWNGRADE_WITH`); FBN1 PM1_Strong suppresses PS1/PM5
  (`PM1_STRONG_EXCLUDES`); FBN1 N>S and G>A substitutions do not receive PM1
  (`PM1_TOLERATED_SUBSTITUTIONS`).
- **PM5 for protein-truncating variants (PM5_PTC)** (`criteria/pm5_ptc.py`):
  nonsense/frameshift variants receive PM5 where the VCEP specifies it —
  RUNX1 (variant downstream of c.98 of NM_001754.4, SpliceAI ≤ 0.20, not with
  PM1), CDH1 (predicted NMD, or PTC upstream of p.Glu836Ter — the most 3'
  truncating variant the CDH1 VCEP classified Pathogenic in the eRepo), ATM (PTC upstream of p.Arg3047) and PALB2 (PTC
  upstream of p.Tyr1183) at Supporting; BRCA1/BRCA2 by the exon containing the
  termination codon, from ENIGMA Specifications Table 4 V1.2 (Strong / Moderate /
  Supporting / not applicable, last-exon codon cutoffs p.1854 / p.3309), applied
  only when PVS1 is met. The PTC is the new stop codon for frameshifts. Splice
  variants (ATM/PALB2 need RNA evidence) and the SCID VCEP point system are not
  automated. `ACMG_PM5_PTC=false` restores missense-only PM5.
- **2015 classifier = Richards et al. 2015 Table 5 as written**: 1 BS + 1 BP is now
  Likely benign (was Benign); a single BS alone is VUS (was Likely benign); when a
  pathogenic and a benign rule are both met the result is VUS (pathogenic rules
  previously took precedence). BA1 still always gives Benign.
- **Bayesian classifier, BA1**: BA1 no longer overrides the point sum; it contributes
  −8 points and the class follows the final score (BA1 alone is still Benign).
- **BA1/BS1 variant exceptions**: new `ba1_exclude` column, and `bs1_exclude` now
  accepts MANE `c.` changes as well as `p.` changes. The LGMD VCEP "BS1 exception"
  alleles (DYSF c.2697+1G>A = NM_003494.3:c.2643+1G>A, ANO5 c.191dup, CAPN3
  c.1746-20C>G and c.2120A>G) receive neither BA1 nor BS1.
- **BP4 VCEP scope**: new `bp4_scope` (`not_applicable` DCLRE1C/IL2RG/IL7R/JAK3/RAG1;
  `splice_only` CDH1/PALB2; `missense_only` HRAS/KRAS/SOS1; `no_synonymous` PAH),
  `bp4_missense_splice_max` (GP9 SpliceAI = 0, also for synonymous/intronic via
  `bp4_splice_cutoff` 0; RUNX1 missense SpliceAI ≤ 0.20) and `bp4_exclude_c`
  (PTEN intron 1 donor/acceptor). BMPR2 BP4 now requires all three of CADD ≤ 22.7,
  AlphaMissense ≤ 0.169 and REVEL ≤ 0.29 (was 2 of 3).
- **BP4, all genes**: splice-prediction BP4 is no longer applied to variants in
  introns 5' of the coding start (5'UTR introns, `c.-N+M` / `c.-N-M`; e.g. GAA
  c.-32-13T>G).
- **PVS1**: applied only when LoF is an established mechanism (VCEP, ClinGen
  haploinsufficiency score 3, or null variants ≥20% of the gene's ≥5 ClinVar
  P/LP variants — `ACMG_PVS1_MIN_NULL_FRACTION` / `ACMG_PVS1_MIN_PLP`; replaces
  the earlier ≥3-null count, selected against the VCEP PVS1 decisions and HI
  scores, Reviewer 2); LOEUF is shown for reference only. NMD is predicted from the PTC position and the MANE v1.5
  exon structure (50-nt rule). The ClinVar-count Moderate caps were removed.
- **PVS1 and gene-disease validity** (ClinGen SVI PVS1 Table 1): for genes
  without a VCEP specification, PVS1 is decreased by two levels (Very Strong →
  Moderate) unless the ClinGen gene-disease validity is Definitive or Strong
  (Moderate, weaker or not curated); `ACMG_PVS1_GDV_DOWNGRADE`.
- **PM4 not with PVS1** at any strength for all genes (SVI PVS1 recommendation).
- **Gene-disease validity**: P/LP capped at VUS when all ClinGen curations for
  the gene are Limited or weaker; uncurated genes are flagged.
- **PS3** text mining counts distinct cited PMIDs and is capped at Supporting;
  **PS4** and **PP1** are applied only from curated sources (supplement or
  ClinVar expert-panel criteria); evidence text names the source.
- **BS2**: gnomAD-based BS2 only for genes whose VCEP applies BS2 with
  population data (VCEP counts, genotypes, sex restriction and strength tiers).
  Genes without a VCEP specification no longer receive gnomAD-based BS2 by any
  route (homozygotes, hemizygotes or carriers): penetrance and age of onset
  cannot be judged automatically, and healthy carriers of adult-onset or
  incompletely penetrant disorders are expected in gnomAD. For them, as for
  genes whose VCEP bars population data, BS2 comes from curated evidence or a
  ClinVar expert-panel (>=3-star) assertion. `ACMG_BS2_AD_MAX_LOEUF` (the
  interim LOEUF gate) was removed and is ignored.
- **BS2 VCEP rules checked against the CSpec text**: ADA (≥3 homozygotes Strong,
  ≥1 Supporting); ETHE1, POLG, SLC19A3, GAMT and SCN1B (≥1 homozygote, "observed
  in a healthy individual"); MTM1 (≥1 hemizygous male or ≥1 homozygous female; new
  `bs2_xl_hom` column); GALT (GALT activity confirmation) and COCH/KCNQ4 (internal
  cohort) set to not applicable. 23 genes now use gnomAD-based BS2, all with
  VCEP-stated thresholds.
- **Input variant types** (Reviewer 1): symbolic CNV/SV alleles, breakends and
  `*` are skipped with the reason in `<output>.not_annotated.tsv`; `explain`
  rejects them; indels > 50 bp are classified with a `LARGE_INDEL` warning.
- **Gene-disease validity cap**: genes with an implemented VCEP specification
  are exempt (`GENE_VALIDITY_NOTE` instead), because their ClinGen curations may
  concern other diseases (PIK3CA, POLG were wrongly capped).
- **PS4 vs PM3**: for genes without a VCEP specification, PS4 is suppressed when
  PM3 is met, with the reason in the PS4 evidence.
- **PS1/PM5 comparators** must have ClinVar review status ≥2 stars by default
  (`ACMG_PS1_MIN_STARS` / `ACMG_PM5_MIN_STARS`; 1 restores the previous policy);
  recorded in the provenance. On the eRepo benchmark with curated evidence this
  raised PS1 F1 0.311→0.384 and PM5 0.470→0.505 and the exact five-tier
  agreement by 0.5–0.8 points (Reviewer 1).
- **`--cspec GENE=CSPEC_ID`** (classify / explain): evaluate RYR1, ACTA1 or VWF
  under one of their disease-specific ClinGen specifications, as in HUVar_app
  (Reviewer 2); default remains the conservative merged specification.
- **TP53 PS3/BS3** from the VCEP functional flowchart on systematic assay data
  (`resources/shared/tp53_functional.tsv`, built by
  `scripts/build_tp53_functional.py` from the NCI TP53 Database R21: Kato 2003,
  Giacomelli 2018, Kotler 2018, Kawaguchi 2005), with the VCEP caveats (no
  SpliceAI-based PP3; PVS1). Priority: supplement > ClinVar expert panel > TP53
  assay data > text mining. **Opt-in** (`--with-tp53-functional`,
  `ACMG_USE_TP53_FUNCTIONAL`): the database terms prohibit monetization, so the
  default stays commercial-use ready; attribution added to NOTICE and README. Matches the TP53 VCEP PS3/BS3 calls for 160 of 172
  eRepo variants.
- **eRepo supplement**: BS2 added for genes whose VCEP bars gnomAD-based BS2
  (`bs2 = not_applicable`; 586 GRCh38 rows), so internal-cohort BS2 no longer
  depends on reading the ClinVar expert-panel record.
- **Expert-panel criterion import (PS3/PS4/PP1/BS2)** now reads only the
  expert-panel SCVs (not other submitters' comments on the same variant) and
  applies clause-level negation ("cannot be applied", "below the threshold",
  "not assessed", "insufficient to apply", …). The previous rule imported PS4/PP1
  that the VCEP had explicitly not applied. Requires rebuilding the ClinVar
  SQLite (`setup_data.py --only clinvar-sqlite`).
- **PM1**: the ClinVar hotspot heuristic for genes without a VCEP definition has
  been **removed** (PM1 from curated evidence only); calibrated against the
  VCEP-defined regions it reached a residue-level F1 of ~0.3.
  `ACMG_PM1_HEURISTIC`, `ACMG_PM1_WINDOW` and `ACMG_PM1_MIN_PATH_VARIANTS` are
  ignored, and the provenance no longer records them.
  IL2RG PM1 requires PM2 and no BA1/BS1/BS2 (`PM1_CO_REQUIREMENTS`).
- **PM1**: VCEP definitions added for genes of multi-gene specifications that the
  CSpec miner skipped — Hearing Loss (PM1 only for the KCNQ4 pore region; the
  other 10 genes not applicable), Brain Malformations (MTOR, AKT3, PIK3CA,
  PIK3R2 domains, Supporting), IL2RG (PM1_Strong residues; transmembrane 263-283 only when a
  charged/polar residue is introduced, new `alt_aa` column), HBA2 (AHSP-binding
  residues) and HBB (non-coding only → not applicable to missense).
- **PP2** fallback thresholds re-derived against the VCEP PP2 decisions
  (131 genes, gene-level train/test split): ≥5 P/LP missense (was 10) and a
  Z-branch benign ceiling of 30% (was 15%).
- **BS2 / PAH**: the gnomAD homozygote route is withheld for PAH (GN006 requires
  homozygosity "in a healthy adult", i.e. a normal phenylalanine level, which
  gnomAD cannot establish); BS2 for PAH comes only from curated sources.
- **PS1/PM5**: VCEP-classified comparators where the VCEP requires them
  (`vcep_comparator_rules.tsv`, eRepo snapshot 2026-05-28).
- **PM1** fallback counts only other P/LP missense; **PP2** fallback adds gnomAD
  missense meeting the gene's BS1 threshold to the benign side
  (`pp2_gene_stats.tsv`); both configurable.
- **PP3/BP4**: note when the gene's VCEP predictor is not active in the run.
- Only **Released** ClinGen specifications are used; the criteria provenance
  table (S1) reports `cspecStatus`.
- Supplement rows from the user take priority over eRepo rows.

### Added

- Data provenance: `data_manifest.json` written by `setup_data.py`,
  `<output>.provenance.json` per run, data versions in `acmg-classify status`.
- `setup_data.py --force-clingen` (ClinGen dosage / gene-validity tables).
- `--exclude-self-expert-panel` (benchmark mode).
- ClinVar SQLite columns `functional_pmids` and `ep_criteria` (rebuild the
  SQLite with `--only clinvar-sqlite` after moving the old file aside).

### Performance

- **Batched, connection-reuse annotation.** The per-variant annotation path was
  reworked to remove the per-variant database setup that dominated wall-time at
  scale (~300k-variant VCFs), especially on network-mounted data. Output is
  byte-identical to the previous per-variant path (covered by new equivalence
  tests):
  - **gnomAD** stats are fetched for the whole batch in a single DuckDB JOIN up
    front (`GnomadDB.precompute()` → `cached()`), replacing a per-variant DuckDB
    connection + query in `_annotate_one`.
  - **Tabix-backed readers** (ClinVar VCF, AlphaMissense, REVEL, RepeatMasker,
    and the opt-in BayesDel / CADD) reuse a persistent thread-local
    `pysam.TabixFile` handle across the batch instead of opening the file per
    variant. BayesDel / CADD handles are created only when their licence-gated
    opt-in is active, so no non-commercial file is opened otherwise.
  - **ESM1b** uses a persistent thread-local read-only SQLite connection for the
    batch instead of `sqlite3.connect()` per variant.
  - Reference-DB file-existence checks are **memoized** for the reader/DB
    lifetime, avoiding a `stat()` per variant per file (6+ tabix files + ESM1b) —
    a real cost on network mounts.

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

- **RYR1 multi-spec AF basis (CSpec overlay).** The per-CSpec side table
  (`disease_prevalence_multispec.tsv`, used by the interactive CSpec switch; the
  conservative batch/CLI row is unaffected) wrongly pinned the RYR1 Congenital
  Myopathy CSpecs (GN150 AD / GN179 AR) to the popmax **point** estimate
  inherited from the Malignant Hyperthermia spec (GN012). Their text specifies
  the **filtering allele frequency (FAF95)**; on few alleles the point AF ≫ FAF95
  and could wrongly fire BS1. The AF metric is now resolved per-CSpec, so both
  congenital-myopathy specs fall back to the FAF95 evaluator default.
- **gnomAD `popmax_pop` label determinism.** On a popmax-AF tie between the
  exomes and genomes rows (GRCh37), `_merge_rows` picked the population label by
  list position, so the single-variant `query()` and the batch `precompute()`
  JOIN could report different `gnomad_popmax_pop` values for the same variant. A
  deterministic tiebreak (larger AN, then population name) makes both paths — and
  repeat runs — agree. Label-only: AF / FAF95 / nhomalt and all classifications
  were already identical.
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
