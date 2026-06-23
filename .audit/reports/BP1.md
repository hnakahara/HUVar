# BP1 Audit — VCEP specializations vs HUHVar app

Scope: ClinGen criterion BP1 ("missense in a gene where truncating variants are the primary
mechanism → more likely benign"). Source of truth: `.audit/cspec_by_criterion/BP1.md`,
`_spec_index.tsv`, and the per-gene BP1 applicability flag in every
`resources/clingen/cspec_json/GN*.json`. App side: `bp1`, `bp1_target`, `bp1_exclude`,
`bp1_strength`, `bp1_no_splice` in `resources/shared/disease_prevalence.tsv`, evaluated by
`src/acmg_classifier/criteria/benign/bp1.py` (semantics in `criteria/bp_genes.py`).

## Result summary

- BP1 = **Applicable** in cspec for **19 genes** across 18 populated GN specs (RASopathy GoF
  genes target truncating; PALB2/APC missense; BRCA1/BRCA2 broad+Strong). GN001 is the generic
  ACMG baseline (no gene) — not app-relevant.
- All 19 applicable genes are **COVERED** in the app with matching target / strength / exclude /
  no_splice parameters.
- Every other populated spec marks BP1 **Not applicable**; the app encodes `bp1=not_applicable`
  for those genes (BP1 is also never applied by default — `bp1.py` requires an explicit VCEP
  `applicable` flag), so the "decline" decision is faithfully captured.
- **MISSING: 0. MISMATCH: 0.**
- 1 item flagged in VERSION & WEB (KCNQ1 flag-vs-description quirk — app handled correctly).
- Web diff on the 3 most material specs (GN112, GN092, GN089): **no JSON-vs-Web discrepancy.**

## MISSING (cspec BP1 specialization not encoded in app)

| gene | GN | ver | cspec (decision/strength/region) | app value | action |
|------|----|-----|----------------------------------|-----------|--------|
| _(none)_ | | | | | All applicable genes are present |

## MISMATCH (app value disagrees with cspec)

| gene | GN | ver | cspec (decision/strength/region) | app value | action |
|------|----|-----|----------------------------------|-----------|--------|
| _(none)_ | | | | | — |

## COVERED (for completeness)

| gene(s) | GN | ver | cspec (decision/strength/region) | app value |
|---------|----|-----|----------------------------------|-----------|
| SHOC2,NRAS,RAF1,SOS1,SOS2,PTPN11,KRAS,MAP2K1,HRAS,RIT1,MAP2K2,BRAF | GN004 | 1.0.0 | Applicable/Supporting, target=truncating (GoF: LoF is benign) | bp1=applicable, target=truncating |
| SHOC2 | GN038 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| NRAS | GN039 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| RAF1 | GN040 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| SOS1 | GN041 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| SOS2 | GN042 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| PTPN11 | GN043 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| KRAS | GN044 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| MAP2K1 | GN045 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| HRAS | GN046 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| RIT1 | GN047 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| MAP2K2 | GN048 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| BRAF | GN049 | 2.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| MRAS | GN087 | 1.4.0 | Applicable/Supporting, truncating | applicable / truncating |
| RRAS2 | GN127 | 1.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| PPP1CB | GN128 | 1.3.0 | Applicable/Supporting, truncating | applicable / truncating |
| PALB2 | GN077 | 1.2.0 | Applicable/Supporting, missense ("apply to all missense") | applicable / missense |
| APC | GN089 | 2.1.0 | Applicable/Supporting, missense; exclude codon 1021-1035 (β-catenin 1st repeat) | applicable / missense / exclude=1021-1035 |
| BRCA1 | GN092 | 1.2.0 | Applicable/**Strong**, broad (silent/missense/in-frame indel); exclude RING 2-101, coiled-coil 1391-1424, BRCT 1650-1857; SpliceAI≤0.1 | applicable / broad / Strong / exclude=2-101;1391-1424;1650-1857 / no_splice=yes |
| BRCA2 | GN097 | 1.2.0 | Applicable/**Strong**, broad; exclude PALB2-binding 10-40, DNA-binding 2481-3186; SpliceAI≤0.1 | applicable / broad / Strong / exclude=10-40;2481-3186 / no_splice=yes |

## VERSION & WEB

| gene | GN | ver | cspec (decision/strength/region) | app value | action |
|------|----|-----|----------------------------------|-----------|--------|
| KCNQ1 | GN112 | 1.0.0 | JSON/Web flag: Supporting=**Applicable**, but description self-declines: "Not applicable, as pathogenic KCNQ1 variants are not limited to truncating variants, but can be missense as well." | bp1=not_applicable | **No change.** App correctly follows the prose decline over the contradictory flag; classifying KCNQ1 missense as benign would be wrong. Documented quirk only. |

Web verification (≤3 cases, all current spec versions; no JSON-vs-Web discrepancy found):
- GN112 (KCNQ1) v1.0.0 — Web confirms Supporting=Applicable with the self-declining description. Matches JSON; app handling correct.
- GN092 (BRCA1) v1.2.0 — Web confirms Strong=Applicable, domains RING 2-101 / coiled-coil 1391-1424 / BRCT 1650-1857, SpliceAI≤0.1. Exactly matches app `bp1_exclude` / `bp1_strength` / `bp1_no_splice`.
- GN089 (APC) v2.1.0 — Web confirms Supporting=Applicable, exclude codon 1021-1035. Matches app `bp1_exclude`.

All spec versions in `_spec_index.tsv` match the JSON `version` fields used above; no stale-version flags.
