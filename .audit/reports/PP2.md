# PP2 audit — HUHVar app vs ClinGen cSpec

Scope: PP2 (missense in a gene with low benign-missense rate where missense is a
common disease mechanism). Source of truth: `.audit/cspec_by_criterion/PP2.md`
(15 Applicable entries) + per-gene PP2 applicability flags read from every
populated `resources/clingen/cspec_json/GN*.json`. App side: `pp2` / `pp2_requires`
columns of `resources/shared/disease_prevalence.tsv`, consumed by
`src/acmg_classifier/criteria/pp2_genes.py` + `criteria/pathogenic/pp2.py`.

## Method
- Extracted PP2 `applicability` per evidence strength from each GN JSON ruleSet.
  Applicable values in the registry: `Applicable`, `Applicable with VCEP specification`.
  Not-applicable values: `Not applicable`, `Not Applicable for this VCEP`.
- Applied the README resolution rules: most-gene-specific spec wins; specificity
  tie across distinct diseases → conservative `not_applicable`; description-level
  **blanket-negation** ("Not applicable due to…") and **gene-exclusion**
  ("…but not PIK3R2") override the strength-level flag.
- Cross-checked the resolved per-gene decision against the app `pp2` column for
  every gene the app records (131 rows) and every gene any spec covers.

## Result: full coverage

**MISSING: none.** Every gene with a VCEP PP2 decision is recorded in the app.

| gene | GN | ver | cspec PP2 decision (+requires) | app value | action |
|------|----|----|-------------------------------|-----------|--------|
| — | — | — | (no gaps) | — | none |

**MISMATCH: none.** All app `pp2` values agree with the resolved cspec decision,
including the conditional and tie/negation cases.

| gene | GN | ver | cspec PP2 decision (+requires) | app value | action |
|------|----|----|-------------------------------|-----------|--------|
| — | — | — | (no disagreements) | — | none |

### Applicable genes — all COVERED
| gene | GN | ver | cspec decision | app value |
|------|----|----|----------------|-----------|
| PTEN | GN003 | 3.2.0 | Applicable (Supporting) | applicable |
| AKT3 | GN018 | 1.1.0 | Applicable (z>3.09) | applicable |
| MTOR | GN018 | 1.1.0 | Applicable (z>3.09) | applicable |
| PIK3CA | GN018 | 1.1.0 | Applicable (z>3.09) | applicable |
| FBN1 | GN022 | 1.0.0 | Applicable (caveat: needs other path. evidence) | applicable |
| PTPN11 | GN043 | 2.3.0 | Applicable (z>3.09) | applicable |
| MAP2K1 | GN045 | 2.3.0 | Applicable (z>3.09) | applicable |
| BRAF | GN049 | 2.3.0 | Applicable (z>3.09) | applicable |
| GCK | GN086 | 3.1.0 | Applicable (all missense) | applicable |
| TPM1 | GN100 | 1.0.0 | Applicable (HCM-only) | applicable |
| BMPR2 | GN125 | 2.0.0 | Applicable **requires PM2,PP3** | applicable / PM2,PP3 |
| PPP1CB | GN128 | 1.3.0 | Applicable (z>3.09) | applicable |
| DNM2 | GN148 | 1.0.0 | Applicable (z=4.87) | applicable |

`pp2_requires` is recorded only for BMPR2 (`PM2,PP3`) — matches GN125
("PM2_supporting and PP3 must be met"). No other Applicable spec imposes
co-criteria. App enforces this post-hoc (see registry suppression pass; code
in `pp2_genes.py::requires`).

## VERSION & WEB diff (≤3 material cases)

These are the cases where a naive strength-flag read of the JSON would disagree
with the app; all three were confirmed against the live web API
(`.../SequenceVariantInterpretation/id/<GN>`). **No JSON-vs-Web discrepancy.**

| gene | GN | ver | cspec PP2 decision (+requires) | app value | action |
|------|----|----|-------------------------------|-----------|--------|
| KCNQ1 | GN112 | 1.0.0 | Supporting flag = "Applicable", but description **blanket-negates**: "Not applicable due to benign variation throughout KCNQ1 (Z=1.83 < 3)". Web = same. | not_applicable | COVERED — app honours the blanket negation, not the raw flag. Correct. |
| PIK3R2 | GN018 | 1.1.0 | Spec Applicable but **excludes PIK3R2** by name ("applicable to MTOR, PIK3CA and AKT3 but not PIK3R2"). Web = same. | not_applicable | COVERED — app honours the gene exclusion. Correct. |
| ACTA1 | GN147 vs GN169 | 2.0.0 / 1.0.0 | **Specificity tie, distinct conditions**: GN147 Applicable (z=6.09) vs GN169 Not applicable. Web GN169 confirmed Not applicable. | not_applicable | COVERED — app takes the conservative tie-break (README documents `--override ACTA1:pp2=applicable` if HCM-appropriate). |

Version note: all cited specs match the GN versions in
`.audit/cspec_by_criterion/_spec_index.tsv`; no stale-version drift found for PP2.

## RASopathy panel vs single-gene specs (sanity check)
GN004 (panel of 12 RASopathy genes) marks PP2 Applicable for "all RASopathy
genes". The newer single-gene RASopathy specs (GN038–GN049, GN087, GN094,
GN127/128) supersede it. In those single-gene specs only **PTPN11 (GN043),
MAP2K1 (GN045), BRAF (GN049), PPP1CB (GN128)** keep PP2 Applicable; the rest
(SHOC2, NRAS, RAF1, SOS1, SOS2, KRAS, HRAS, RIT1, MAP2K2, MRAS, LZTR1, RRAS2)
are Not applicable in their own specs. The app's per-gene `pp2` values match this
most-specific resolution exactly — the panel's blanket "applicable" is correctly
not propagated to the declined genes.

## Counts
- Populated specs scanned: 131 ruleset rows (123 unique GN specs).
- cSpec Applicable PP2 entries: 15 (1 generic GN001 + 14 gene-bearing). All 13
  app-relevant applicable genes COVERED.
- App `pp2` populated rows: 132 (13 applicable, 119 not_applicable).
- MISSING: 0 · MISMATCH: 0 · Version drift: 0 · JSON-vs-Web discrepancies: 0.
