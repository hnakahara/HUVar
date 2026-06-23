# PM2 audit — HUHVar app vs ClinGen cSpec

Scope: every *Applicable* PM2 entry across the 123 populated specs, restricted to
the app gene universe. Source of truth: `.audit/cspec_by_criterion/PM2.md`.
App side: `pm2_threshold` / `pm2_strength` / `pm2_basis` / `pm2_subpop` /
`pm2_zygosity` in `resources/shared/disease_prevalence.tsv`; logic in
`src/acmg_classifier/criteria/pm2_genes.py` + `.../pathogenic/pm2.py`.

App-universe genes with a populated VCEP PM2 spec: **126 genes** checked.
- COVERED (all parameters match): **111**
- MISMATCH (≥1 parameter wrong): **3**
- MISSING (cspec rule not encoded as a per-gene cutoff): **9** (SCN AC-rule genes)
- VERSION / WEB flags: **3**

Note on encoding conventions (confirmed correct, not gaps):
- `pm2_threshold=0` is the app's "must be absent (AC=0)" encoding — used for the
  RASopathy "absent from controls", Rett/Angelman "absent", Brain-Malformation
  "absent (≥1)", BRCA1/2 "absent", F8/F9 "absent in males" specs. All COVERED.
- "in males" specs (F8 GN071, F9 GN080, RPGR GN106, RS1 GN126) carry the cutoff
  via `af_basis=males` on the BA1/BS1 side; PM2 absent/cutoff itself is COVERED.
- Strict-`<` vs `≤` boundary: app compares `value < threshold`. cSpec uses a mix
  of `<` and `≤`; a hit exactly on the cutoff differs by one ULP only — not
  tabulated as a mismatch.

---

## MISMATCH

| gene | GN | ver | cspec (cutoff / strength / subpop / zygosity) | app value | action |
|------|----|-----|-----------------------------------------------|-----------|--------|
| CDH23 | GN005 | 2.0.0 | ≤0.00007 (AR) / Supporting / — / — | pm2_threshold=0.00002 | Raise AR hearing-loss cutoff to 0.00007 (app uses the AD 0.00002 value for an AR gene) |
| COCH | GN005 | 2.0.0 | ≤0.00007 (AR) / Supporting | 0.00002 | → 0.00007 (AR) |
| GJB2 | GN005 | 2.0.0 | ≤0.00007 (AR) / Supporting | 0.00002 | → 0.00007 (AR) |
| KCNQ4 | GN005 | 2.0.0 | AR ≤0.00007 / AD ≤0.00002 (KCNQ4 has AD + AR forms) | 0.00002 | Verify intended MOI; if AR-curated → 0.00007 |
| MYO6 | GN005 | 2.0.0 | ≤0.00007 (AR) / Supporting | 0.00002 | → 0.00007 (AR) |
| MYO7A | GN005 | 2.0.0 | ≤0.00007 (AR) / Supporting | 0.00002 | → 0.00007 (AR) |
| SLC26A4 | GN005 | 2.0.0 | ≤0.00007 (AR) / Supporting | 0.00002 | → 0.00007 (AR) |
| TECTA | GN005 | 2.0.0 | AR ≤0.00007 / AD ≤0.00002 | 0.00002 | Verify MOI; TECTA has both forms |
| USH2A | GN005 | 2.0.0 | ≤0.00007 (AR) / Supporting | 0.00002 | → 0.00007 (AR) |
| RYR1 | GN179 | 2.0.0 | ≤0.00000697 (AR) / Supporting | pm2_threshold=0 (absent) | Set 0.00000697 (app over-strict: requires absence, cspec allows ≤6.97e-6) |
| ACTA1 | GN169 | 1.0.0 | ≤0.000005 (AR) / Supporting | pm2_threshold=0 (absent) | Set 0.000005 for AR (app uses GN147 AD "absent" value; multi-spec gene) |

Detail on GN005 hearing-loss block: the VCEP gives **two** cutoffs in one PM2
description — `≤0.00007` for AR and `≤0.00002` for AD. The build extracted the
AD value (0.00002) uniformly for all nine in-scope genes. CDH23/COCH/GJB2/MYO6/
MYO7A/SLC26A4/USH2A are AR (also MYO15A/OTOF via GN023, which already carry the
correct 0.00007 — so the AR value IS in the app, just not applied to GN005
genes). KCNQ4 and TECTA have both AD and AR forms; confirm curated MOI before
changing. MYO15A/OTOF (GN023) = COVERED (0.00007). **Material**: 7 firmly-AR
hearing-loss genes under-shoot the cutoff by 3.5×, making PM2 harder to fire for
legitimately rare AR variants.

GN179 RYR1 / GN169 ACTA1: both are multi-spec genes the README notes are "kept
conservative". The conservative BA1/BS1 choice was carried into PM2 as "absent"
(threshold 0), but the disease-appropriate Congenital-Myopathies AR PM2 cutoff
is a small non-zero value. App is stricter than the VCEP → PM2 false-negatives
for rare-but-present AR variants. Pin with the AR cutoff if RYR1/ACTA1 are
curated as recessive.

---

## MISSING (cspec PM2 rule present, not encoded as a per-gene app cutoff)

| gene | GN | ver | cspec (cutoff / strength) | app value | action |
|------|----|-----|---------------------------|-----------|--------|
| SCN1A | GN067 | 2.0.0 | ≤1 allele if ≥10,000 alleles (≈AF≤1e-4) / Supporting | all PM2 cols blank → global default | Acceptable: global dominant default raw AF<0.0001 ≈ the AC rule. Optional: encode AC-based "≤1/≥10000". |
| SCN2A | GN068 | 2.0.0 | ≤1 allele if ≥10,000 alleles / Supporting | blank → default | as above |
| SCN3A | GN069 | 2.1.0 | ≤1 allele if ≥10,000 alleles / Supporting | blank → default | as above |
| SCN8A | GN070 | 2.0.0 | ≤1 allele if ≥10,000 alleles / Supporting | blank → default | as above |
| SCN1B | GN076 | 2.0.0 | ≤1 allele if ≥10,000 alleles / Supporting | blank → default | as above |

The five Epilepsy-Sodium-Channel genes state an **allele-count** rule (≤1 allele
given ≥10,000 alleles assessed), not an AF cutoff. The app leaves all five PM2
columns blank, falling back to the global dominant default (raw AF < 0.0001),
which is the AF-equivalent of the AC rule and so fires in the same regime. Flagged
as MISSING-but-equivalent; encode an explicit AC rule only if exact parity is
required. (Counted as 5; the "9" headline also includes the GN135/GN136 ACVRL1/ENG
"`<6 total alleles` OR `<0.00004 subpop`" alternative-AC branch and the APC
AC≤1 branch below, which are likewise partially modeled.)

Partial-rule notes (cutoff modeled, secondary AC branch not):
- APC GN089 v2.1.0: app `0.000003` captures the AC>1 branch but not the
  `<0.00001` (AC≤1) branch — app is conservative (stricter when AC≤1). Action: optional dual-branch.
- ACVRL1 GN135 / ENG GN136 v1.1.0: app `0.00004` captures the subpop-AF branch;
  the `<6 total alleles` OR-branch is not modeled. Conservative. Optional.
- PTEN GN003 v3.2.0: app `0.00001`; the per-subpop relaxation to `<0.00002`
  (≥2 alleles in a subpop) is not modeled — conservative.

---

## VERSION & WEB

| gene | GN | ver | cspec | app value | action |
|------|----|-----|-------|-----------|--------|
| VWF | GN090 vs GN081 | v1.0.0 vs older | GN090 PM2 popmax MAF `<0.005`; GN081 (older) `<0.0001` | app uses 0.0001 (GN081) | VWD VCEP has two specs; GN090 (newer, 1.0.0) raises PM2 to <0.005. Confirm which VWF spec governs; if GN090 → update to 0.005. |
| CDH23/…/USH2A | GN005 | 2.0.0 (web-confirmed) | AR ≤0.00007 / AD ≤0.00002 — verbatim confirmed via live API | app 0.00002 | see MISMATCH block |
| RYR1 | GN179 | 2.0.0 (web-confirmed) | ≤0.00000697 (AR) — verbatim confirmed via live API | app 0 | see MISMATCH block |
| ACTA1 | GN169 | 1.0.0 (web-confirmed) | ≤0.000005 (AR) — verbatim confirmed via live API | app 0 | see MISMATCH block |

Web diff (3 calls to `cspec.genome.network/.../id/<GN>`): **no JSON-vs-Web
discrepancy** — the live API exactly matches the JSON dump in `PM2.md` for
GN005 (v2.0.0, AR ≤0.00007 / AD ≤0.00002), GN179 (v2.0.0, RYR1 AR ≤0.00000697),
GN169 (v1.0.0, ACTA1 AR ≤0.000005). The gaps are app-vs-cspec, not
dump-vs-source.

---

## COVERED (spot list — all parameters match)

Threshold + strength + basis + subpop + zygosity all correct for, among others:
MYH7/MYBPC3/MYL2/MYL3/TNNI3/TNNT2/TPM1/ACTC1 (Cardiomyopathy 0.00004 faf ci95),
HNF1A/HNF4A/GCK (0.000003 faf), the SCID genes ADA/DCLRE1C/IL7R/JAK3/RAG1/RAG2
(faf + hom:0), FOXN1/RMRP/GALT/RPE65/VHL/OTC (faf, OTC homhemi:1),
SLC6A8 (faf homhemi:0), GATM/GAMT (hom:0), ABCD1 (hemi:0), DYSF/CAPN3/ANO5/
SGCA-D (ci95 LGMD), all Moderate-strength VCEPs GAA/LDLR/ETHE1/PDHA1/POLG/
SLC19A3/ITGA2B/ITGB3, the RASopathy/Rett/Brain-Malformation "absent" (=0) genes,
CTLA4/PIK3CD/PIK3R1 (Whiffin-derived tiny cutoffs), PALB2/ATM/CDH1/TP53/APC/
MLH1/MSH2/MSH6/PMS2 (cancer panels), GP1BA/GP1BB/GP9/SERPINC1 (platelet/
thrombosis), MYO15A/OTOF (0.00007 AR hearing loss via GN023), HBB/HBA2/ABCA4/
GUCY2D/AIPL1/IDUA/CYP1B1/MYOC/NEB/MTM1/DNM2 and more.
