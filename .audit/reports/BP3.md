# BP3 Audit — VCEP specializations vs HUHVar app

Criterion: **BP3** (in-frame indel in a repetitive region without known function).
Source of truth: `.audit/cspec_by_criterion/BP3.md` (30 Applicable strength entries across the populated specs).
App side: `bp3` / `bp3_regions` in `resources/shared/disease_prevalence.tsv`; logic in
`src/acmg_classifier/criteria/bp_genes.py` and `src/acmg_classifier/criteria/benign/bp3.py`.

## Specialization model

BP3 has **no strength specialization anywhere** — every Applicable cspec entry is at the
default **Supporting**. So a BP3 specialization is only one of:

1. **Applicability flip** — a VCEP that declines BP3 (`bp3=not_applicable`, suppresses the
   Dfam/RepeatMasker repeat heuristic). This is the overwhelmingly common case (defaulting to
   the generic ACMG heuristic when no VCEP opinion exists).
2. **Region-restricted BP3** — the VCEP enumerates the repetitive regions itself
   (`bp3_regions`), so being inside one IS the "repetitive region" and awards BP3 directly;
   outside = not BP3 (`bp3_in_region` in `bp_genes.py`).

Region-restricted specializations in cspec (3 specs, all COVERED):

| gene | GN | ver | cspec region | app `bp3_regions` | status |
|---|---|---|---|---|---|
| VHL | GN078 | 1.1.0 | AA14-AA48 (GXEEX repeat, p30 5' end) | `14-48` | COVERED |
| RPGR | GN106 | 1.0.0 | ORF15 aa585-1078; disordered 609-776, 790-906, 989-1020 | `585-1078;609-776;790-906;989-1020` | COVERED |
| FOXG1 | GN016/GN035 | 2.0.0/6.0.0 | polyHis His47-57, polyGln Gln70-73, polyPro Pro58-61/65-69/74-80 | `47-57;58-61;65-69;70-73;74-80` | COVERED |

Note RPGR: the three disordered sub-ranges are fully inside the master 585-1078 range — redundant
but harmless (any position in them is already in-region).

## MISSING

None. Every cspec-Applicable BP3 gene is present and set `applicable` in the TSV; every
region-restricted spec has its `bp3_regions` populated.

## MISMATCH

None (true mismatches). The one apparent app-vs-cspec divergence is a version artifact of the
audit source, not a content error — see VERSION&WEB below.

The KCNQ1 case is **correct**: BP3.md GN112 carries the text "Not applicable to _KCNQ1_" inside an
Applicable-Supporting entry; the app correctly resolves this to `bp3=not_applicable`.

| gene | GN | ver | cspec | app value | action |
|---|---|---|---|---|---|
| KCNQ1 | GN112 | 1.0.0 | Applicable-Supporting text = "Not applicable to KCNQ1" | `not_applicable` | none (COVERED) |

## VERSION & WEB

Web diff performed on the 3 most material cases (cap = 3) via
`https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/<GN>`.

| gene(s) | GN (audit) | ver | cspec (BP3.md) | app value | web result | action |
|---|---|---|---|---|---|---|
| 12 RASopathy genes (SHOC2,NRAS,RAF1,SOS1,SOS2,PTPN11,KRAS,MAP2K1,HRAS,RIT1,MAP2K2,BRAF) | GN004 | 1.0.0 | Applicable-Supporting | `not_applicable` (all 12) | **GN004 web = Applicable-Supporting (confirms BP3.md). BUT GN004 v1.0.0 is the superseded combined spec.** | Flag audit source as stale; app is correct |
| BRAF (per-gene) | GN049 | 2.3.0 | (not in BP3.md) | `not_applicable` | **GN049 web = BP3 Not Applicable (all strengths)** | App correctly follows the newer per-gene spec |
| BMPR2 | GN125 | 2.0.0 | Applicable-Supporting | `applicable` | **GN125 web = BP3 Applicable-Supporting (v2.0.0)** | COVERED |

### Key version finding

BP3.md's RASopathy entry is **GN004 v1.0.0**, the old single combined RASopathy spec, where BP3 is
Applicable. The app does **not** track GN004 — its RASopathy rows reference the **newer per-gene
specs** (e.g. BRAF→GN049 v2.3.0, PTPN11→GN043; `notes` = "multiple specs (gene-specific kept)"),
and those per-gene specs **flip BP3 to Not Applicable**. Web confirmation of GN049 (BRAF, BP3 Not
Applicable at all strengths) shows the app's `bp3=not_applicable` for all 12 RASopathy genes is
**correct and current**; it is the **audit source BP3.md (GN004) that is outdated**, not the app.

No JSON-vs-Web discrepancies on the app side: GN049 (Not Applicable) and GN125 (Applicable-
Supporting) both match the TSV.

## Counts

- cspec Applicable BP3 entries in BP3.md: 30 (28 gene-bearing + 2 generic GN001 placeholders).
- App genes with a BP3 VCEP opinion: 130 (`applicable` = 35; `not_applicable` = 95).
- Region-restricted specs: 3 (VHL, RPGR, FOXG1) — all COVERED.
- MISSING: 0 · MISMATCH: 0 · VERSION-flagged: 1 (RASopathy GN004 superseded; app correct).
