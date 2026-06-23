# BA1 audit — HUHVar app vs ClinGen cSpec

Scope: every populated cSpec spec with an Applicable BA1 (Stand Alone) entry,
cross-checked against `resources/shared/disease_prevalence.tsv`
(`ba1_threshold`, `af_basis`). Released / Approved-For-Release specs only are
expected to carry an app cutoff; Pilot / In-Prep / Submitted specs legitimately
have no row value (flat default 0.05 applies). The app gene universe is fully
covered (every in-scope gene has a row; non-listed genes correctly absent).

Counts (released specs with a defined cSpec BA1):
- COVERED: 130
- MISMATCH: 1 (RPGR)
- MISSING: 0
- af_basis=males: 2 rows (RPGR, RS1) — both present; RS1 correct, RPGR value wrong (see MISMATCH)

A "specialization" = any BA1 cutoff other than the flat 0.05 default, or
`af_basis=males`. Of 130 covered released genes, **128** carry a non-default
(specialized) cutoff; only CYP1B1 (GN104) and RPGR (GN106) legitimately use the
0.05 figure — and GN104 is a genuine gnomAD ≥0.05 cutoff, while GN106 is the
legacy ESP/1000G/ExAC "5% in males" (kept because it is the only value).

## MISSING

| gene | GN | spec ver | cspec cutoff | app ba1_threshold | action |
|------|----|----------|--------------|-------------------|--------|
| — (none) | | | | | All released in-scope genes have a populated `ba1_threshold`. |

Genes whose only spec is Pilot/In-Prep/Submitted/Deleted (ACO2, CACNA1F,
CEP290, CHEK2, CHM, CPS1, CYP27A1, DDC, DDX41, DNAH5, F11, F7, FKRP, FOXC1,
FOXP3, G6PD, GATA4, GBA1, GDAP1, GFAP, GRIN1, GRIN2A, GRIN2B, GRIN2D, HBA1,
HNF1B, INS, LRRK2, MFN2, MUTYH, MYH6, NAGS, NDP, NF1, NKX2-5, NR2F2, OFD1,
OPA1, PAX6, PCCA, PCCB, PEX1, PKD1, PKD2, PROC, PROS1, PRPH2, RAD51C, RHO,
RYR2, SDHB, SERPING1, SOD1, SORD, SPRED1) correctly have a blank
`ba1_threshold` → flat default 0.05. Not counted as MISSING because the cSpec
has no released BA1 number yet. (GFAP/GN157 is "Approved For Release" but the
slice carries no BA1 entry; blank is acceptable.)

## MISMATCH

| gene | GN | spec ver | cspec cutoff | app ba1_threshold | action |
|------|----|----------|--------------|-------------------|--------|
| RPGR | GN106 | 1.0.0 | 0.05 (males) | 0.00005 (males) | **Fix.** cSpec BA1 = "Allele frequency **in males** is above **5%** in ESP / 1000 Genomes / ExAC in the subpopulation with the highest frequency" → per README parsing (RPGR-style legacy 5% is *kept* when it is the only value), BA1 should be **0.05** with `af_basis=males`. The app row has a `manual override` to **0.00005**, which equals the gene's BS1 male cutoff — i.e. the override appears to have copied the BS1 value into BA1. Set `ba1_threshold=0.05` (keep `af_basis=males`), or document why a non-spec stricter BA1 is intentional. Note: 0.00005 is far stricter than 0.05, so this *over-fires* BA1 (more aggressive benign calls) relative to the VCEP. |

RS1 (GN126, v1.0.0): cSpec BA1 = "≥2×10⁻⁴ in males" = 0.0002; app `ba1_threshold=0.0002`, `af_basis=males` → COVERED (the correctly handled "in males" companion to RPGR).

## VERSION & WEB

Web diff performed on the single material discrepancy (1 fetch, cap ≤3):

| gene | GN | spec ver (app/JSON) | web ver | cspec cutoff (JSON) | web cutoff | result |
|------|----|---------------------|---------|---------------------|------------|--------|
| RPGR | GN106 | 1.0.0 | 1.0.0 | 0.05 (males) | 0.05 (males) | JSON ⟷ Web **agree** (both 0.05, "in males"). The deviation is purely on the **app** side (manual override to 0.00005). No version skew. |

No other gene showed an app value implying an older spec version than the JSON
slice. All released app cutoffs match the spec version recorded in the slice
(within rounding), including parsing edge cases verified against the README
rules:
- Rett/Angelman sub-population rule: app uses 0.000083 (latest GN032-037), not the 0.0005 headline. ✔
- KCNQ1 legacy-5% dropped → 0.004 (gnomAD). ✔
- CYP1B1 genuine gnomAD ≥0.05 kept. ✔
- Hearing-loss GN005 AR genes 0.005 / AD genes 0.001 split applied per gene. ✔
- Multi-spec RASopathy / cardiomyopathy / SCID gene-specific cutoffs all match. ✔

---

### Summary

- **131** Applicable BA1 entries across 123 populated specs reviewed; **130**
  released in-scope genes matched the app within rounding (**COVERED**).
- **0 MISSING** — every released in-scope gene has a populated `ba1_threshold`;
  pilot-only genes correctly fall back to the 0.05 default.
- **1 MISMATCH — RPGR (GN106):** app `ba1_threshold=0.00005` (manual override)
  vs cSpec/Web **0.05 in males**. The override copied the BS1 male value into
  BA1; recommend resetting BA1 to 0.05 (basis=males) unless intentionally
  stricter.
- `af_basis=males` correctly set for both X-linked "in males" genes (RPGR, RS1);
  RS1 value correct, RPGR value wrong (above).
