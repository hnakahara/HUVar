# BP7 Audit — VCEP specializations vs app (bp7_phylop / bp7_intronic + bp7.py)

Source of truth: `.audit/cspec_by_criterion/BP7.md` (160 Applicable BP7 entries, 123 populated specs).
App side: `resources/shared/disease_prevalence.tsv` columns `bp7_phylop` / `bp7_intronic`;
`src/acmg_classifier/criteria/benign/bp7.py`; `src/acmg_classifier/criteria/bp_genes.py`.
Global default: phyloP `bp7_phylop_max = 2.0`; intronic default = Walker deep-intronic (donor ≥ +7, acceptor ≤ −21).

App model has only 4 BP7 levers: numeric phyloP cutoff, `na` (conservation non-informative), `noncanonical` intronic (|dist| ≥ 3), or blank (Walker +7/−21). It has **no** lever for VCEP-specific intronic ranges between those two modes (e.g. −4/+7), and **no** per-gene SpliceAI cutoff (the splice gate is hard-coded 0.10 / SQUIRLS 0.20 globally).

---

## Summary counts

| Class | Count | Notes |
|---|---|---|
| COVERED (phyloP cutoff) | 13 | 0.1, 0.2, 1.5, 2.0, 0 cutoffs all transcribed correctly |
| COVERED (conservation `na`) | 19 | TP53, BRCA1/2, RUNX1, MYOC, LCA, SCID, GALT, ABCD1 |
| COVERED (intronic `noncanonical`) | ~25 | RASopathy/PIK3 panels, VHL, MYOC, RUNX1, DICER1, BMPR2, MLH1/MSH2/MSH6/PMS2 |
| **MISMATCH** | **9** | Cardiomyopathy −4/+7 intronic (8 genes) + SLC6A8 −4/+7 |
| **MISSING (na candidates)** | **0 hard** | InSiGHT/PTEN/CDH1 omit conservation but do not declare it non-informative → NOTE only |
| VERSION&WEB | 0 discrepancies | 3 web-checked specs match TSV/index version and parameters |

Top gap: a single systematic **intronic-range MISMATCH** affecting the 8 Cardiomyopathy genes (and SLC6A8): their VCEPs extend BP7 to the acceptor side at **−4 outward**, but the app applies the Walker **−21** default, so intronic variants at −4..−20 are wrongly excluded from BP7.

---

## MISMATCH

App value differs from cspec. `intronic(app)` = effective eligible range.

| gene | GN | ver | cspec (phyloP / intronic) | app value (phylop / intronic) | action |
|---|---|---|---|---|---|
| MYH7   | GN002 | 2.0.0 | default / **intronic −4 and +7 outward** | 2.0 default / blank → ≥+7 or ≤−21 | Extend intronic to acceptor −4; need a mode/value between default and `noncanonical` |
| MYBPC3 | GN095 | 1.0.0 | default / intronic −4/+7 outward | blank → ≥+7 or ≤−21 | same as MYH7 |
| TNNI3  | GN098 | 1.0.0 | default / intronic −4/+7 outward | blank → ≥+7 or ≤−21 | same |
| TNNT2  | GN099 | 1.0.0 | default / intronic −4/+7 outward | blank → ≥+7 or ≤−21 | same |
| TPM1   | GN100 | 1.0.0 | default / intronic −4/+7 outward | blank → ≥+7 or ≤−21 | same |
| ACTC1  | GN101 | 1.0.0 | default / intronic −4/+7 outward | blank → ≥+7 or ≤−21 | same |
| MYL2   | GN102 | 1.0.0 | default / intronic −4/+7 outward | blank → ≥+7 or ≤−21 | same |
| MYL3   | GN103 | 1.0.0 | default / intronic −4/+7 outward | blank → ≥+7 or ≤−21 | same |
| SLC6A8 | GN027 | 2.0.0 | default / **intronic beyond −4 or +7** | blank → ≥+7 or ≤−21 | Extend intronic to acceptor −4 (same fix); SpliceAI ≤0.10 already matches |

Root cause: `bp7.py` `_intronic_eligible` supports only `_DEEP_INTRONIC_ACCEPTOR_MAX = -21` (default) or `noncanonical` (|dist| ≥ 3). The −4/+7 VCEP range (donor ≥+7 = default, acceptor ≤−4) is representable by neither. `noncanonical` would be *too permissive* on the donor side (admits +3..+6). Recommend a parametric acceptor/donor cutoff column, or a new `cardio` intronic mode.

Note (not counted): SQUIRLS/global splice gate is fixed at 0.20 / SpliceAI 0.10. Several VCEPs set different SpliceAI cutoffs (DYSF/LGMD ≤0.05, GP genes / HNF / KCNQ1 / RPGR 0.2, VWF GN090 = 0). These are not expressible in the BP7 columns at all — out of scope for `bp7_phylop`/`bp7_intronic` but worth a separate splice-cutoff column. KCNQ1/HNF1A/HNF4A/GCK use SpliceAI <0.2 with phyloP <2.0; app phyloP=2.0 is correct but its splice gate 0.10 is stricter than the VCEP's 0.2 (conservative, not a false-benign risk).

---

## MISSING

No hard MISSING parameters: every explicit phyloP cutoff, `na` declaration, and `noncanonical` extension in cspec is present in the TSV.

NOTE-level (conservation omitted, not declared non-informative — left as default 2.0, defensible):

| gene | GN | ver | cspec (phyloP / intronic) | app value | action |
|---|---|---|---|---|---|
| APC | GN089 | 2.1.0 | no conservation clause / intronic +7/−21 | 2.0 default / blank | Consider `na`? cspec omits but does not declare non-informative — keep default |
| MLH1 | GN115 | 2.0.0 | no conservation clause / intronic +7/−21 | 2.0 / noncanonical | intronic covered; conservation omission only |
| MSH2/MSH6/PMS2 | GN137-139 | 2.0.0 | no conservation clause / +7/−21 | 2.0 / noncanonical | same as MLH1 |
| PTEN | GN003 | 3.2.0 | no conservation clause / +7/−21 | 2.0 / blank | omission only |
| CDH1 | GN007 | 3.1.0 | no conservation clause / +7/−21 | 2.0 / blank | omission only |

These VCEPs simply do not mention conservation; they do not state it is non-informative (unlike TP53/SCID/GALT/RUNX1 which use "no requirement"/"not required"). Marking them `na` would over-apply BP7. Recommend leaving default; flagged for curator review only.

---

## VERSION & WEB

3 most-material cases web-checked against `https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/<GN>`. No JSON-vs-Web discrepancies; all versions match `_spec_index.tsv` and TSV.

| gene | GN | ver (index/web) | cspec (phyloP / intronic) | app value | result |
|---|---|---|---|---|---|
| MYH7 | GN002 | 2.0.0 / 2.0.0 | default / intronic −4 and +7 outward (web confirmed verbatim) | 2.0 / blank (≥+7,≤−21) | version OK; confirms intronic MISMATCH above |
| VHL  | GN078 | 1.1.0 / 1.1.0 | phyloP ≤0.2 / silent or intronic if BP4 met | 0.2 / noncanonical | COVERED; web confirms 0.2 |
| FOXN1| GN113 | 2.3.0 / 2.3.0 | default / conservation **required** (web confirmed "not highly conserved" mandatory) | 2.0 / blank (no `na`) | COVERED; correctly NOT `na` (distinct from sibling SCID na genes) |
