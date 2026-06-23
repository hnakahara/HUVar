# BP4 audit — gene-specific REVEL cutoffs vs `disease_prevalence.tsv`

Scope: ClinGen cSpec BP4 specializations anchored on a **REVEL upper cutoff**
(`REVEL ≤ value` fires BP4). App columns: `revel_bp4_supporting/moderate/strong`.
Source of truth: `.audit/cspec_by_criterion/BP4.md`. Convention (per
`resources/clingen/README.md`): firing edge = **upper** bound of a band;
single-tier capping; most-specific spec wins; specs citing REVEL without a
number leave columns blank → genome-wide Bergquist 2024 default.

Result: **~80 REVEL-anchored genes COVERED**, **3 MISMATCH**, **1 MISSING**.
All four flagged genes are coagulation/platelet/thrombosis specs whose
"REVEL ≤ X AND SpliceAI ≤ 0.1" text was mis-mined — the SpliceAI 0.1 (or a
blank) landed in the REVEL column instead of the stated REVEL cutoff.

---

## MISSING

App row exists but `revel_bp4_*` is blank, while cspec states a REVEL cutoff.

| gene | GN | ver | cspec REVEL cutoff (tier) | app value | action |
|------|----|-----|---------------------------|-----------|--------|
| GP9 | GN083 | 1.1.0 | ≤0.290 (Supporting) | *(blank)* → default 0.290 | set `revel_bp4_supporting=0.29` |

Note: the Bergquist default Supporting is also 0.290, so GP9's *effective*
behaviour is unchanged — but the column should be filled for provenance and to
make it explicit/robust to default changes. Sibling platelet genes GP1BA
(GN079) and GP1BB (GN082) have identical text and **are** populated (0.29);
GP9 was simply skipped.

---

## MISMATCH

App REVEL Supporting cutoff disagrees with cspec (both numbers shown).

| gene | GN | ver | cspec REVEL cutoff (tier) | app value | action |
|------|----|-----|---------------------------|-----------|--------|
| F8 | GN071 | 2.0.0 | **≤0.3** (Supporting) | **0.1** | fix `revel_bp4_supporting` 0.1 → 0.3 |
| F9 | GN080 | 2.1.0 | **≤0.3** (Supporting) | **0.1** | fix `revel_bp4_supporting` 0.1 → 0.3 |
| SERPINC1 | GN084 | 1.1.0 | **≤0.30** (Supporting) | **0.1** | fix `revel_bp4_supporting` 0.1 → 0.3 |

Root cause: each spec reads "REVEL score of 0.3 or below **AND** a SpliceAI
score ≤ 0.1". The miner captured the **SpliceAI 0.1** as the REVEL cutoff. Net
clinical effect: BP4 is **3× too strict** for these genes (only fires at
REVEL ≤0.1 instead of ≤0.3), under-calling benign-supporting evidence and
missing true benign missense calls in REVEL (0.1, 0.3].

---

## NON-REVEL (coverage boundary — out of scope for `revel_bp4_*`)

BP4 specializations whose computational anchor is **not** REVEL. These are
correctly *not* in the REVEL columns; listed as the audit boundary. App handles
SpliceAI BP4 generically in `bp4.py` (≤0.10 Supporting), but the gene-specific
BayesDel / CADD / AlphaMissense / HCI-prior / non-0.10 SpliceAI cutoffs below
are **not** encoded per-gene.

| gene(s) | GN | ver | non-REVEL BP4 anchor |
|---------|----|-----|----------------------|
| TP53 | GN009 | 2.4.0 | BayesDel ≤ -0.008 (Mod) / < 0.16 (Sup); SpliceAI <0.2 |
| BRCA1 | GN092 | 1.2.0 | BayesDel no-AF ≤0.15 AND SpliceAI ≤0.1 (domain-gated) |
| BRCA2 | GN097 | 1.2.0 | BayesDel no-AF ≤0.18 AND SpliceAI ≤0.1 (domain-gated) |
| CDH1 | GN007 | 3.1.0 | splicing predictors only (3 in agreement) |
| RUNX1 (synon/intron) | GN008 | 3.1.0 | SpliceAI ≤0.20 (REVEL <0.5 also present → covered) |
| AKT3,MTOR,PIK3CA,PIK3R2 | GN018 | 1.1.0 | splicing tools only (2/3 of varSEAK/SpliceAI/MaxEntScan) |
| APC | GN089 | 2.1.0 | missense N/A; synon/intron splicing predictors |
| PALB2 | GN077 | 1.2.0 | missense "do not use"; SpliceAI ≤0.1 only |
| VHL | GN078 | 1.1.0 | missense predictors barred; SpliceAI ≤0.1 + VarSeak |
| MLH1 | GN115 | 2.0.0 | HCI-prior <0.11; SpliceAI ≤0.1 |
| MSH2 | GN137 | 2.0.0 | HCI-prior <0.11; SpliceAI ≤0.1 |
| MSH6 | GN138 | 2.0.0 | HCI-prior <0.11; SpliceAI ≤0.1 |
| PMS2 | GN139 | 2.0.0 | MAPP/PP2 prior <0.11; SpliceAI ≤0.1 |
| CTLA4 | GN122 | 1.0.0 | REVEL <0.25 **AND** CADD <20 (REVEL covered; CADD extra) |
| PIK3CD | GN141 | 1.0.0 | REVEL ≤0.290 **AND** CADD ≤22.7 (REVEL covered; CADD extra) |
| PIK3R1 | GN160 | 1.0.0 | REVEL ≤0.290 **AND** CADD ≤21.5 (REVEL covered; CADD extra) |
| BMPR2 | GN125 | 2.0.0 | 2-of-3 REVEL≤0.29/AlphaMissense≤0.169/CADD≤22.7 (REVEL covered) |
| ABCA4 | GN164 | 1.0.0 | REVEL covered; synon/indel via CADD ≤17.3/17.4-20 |
| HBB / HBA2 | GN170 / GN173 | 1.0.0 | REVEL <0.7 **AND** SpliceAI ≤0.3 (REVEL covered; note SpliceAI 0.3 not 0.1) |

Non-REVEL-only specializations (no REVEL component at all): **9** specs
(TP53, BRCA1, BRCA2, CDH1, GN018 brain-malformations panel, APC, PALB2, VHL,
and the 4 InSiGHT/Lynch HCI-prior genes MLH1/MSH2/MSH6/PMS2).
Specs with a REVEL component **plus** an additional non-REVEL tool (CADD /
AlphaMissense / a non-0.10 SpliceAI): CTLA4, PIK3CD, PIK3R1, BMPR2, ABCA4,
HBB, HBA2 — REVEL side is covered; the second-tool agreement requirement is
not modelled per-gene.

---

## VERSION & WEB

### Version flags
- **PAH (GN006 v2.0.0, Approved For Release — not Released).** cSpec grants BP4
  at **Strong + Moderate + Supporting** ("Applicable as described in Pejaver"),
  but only Supporting cites a REVEL band (0.183–0.290). App: `sup=0.29` only →
  PAH capped at Supporting. Defensible (no REVEL number for the higher tiers →
  default), but the multi-tier intent is not reflected. Re-mine if/when the
  Strong/Moderate tiers gain explicit REVEL numbers.
- **Rett/Angelman genes** correctly resolved to the newer v6/v7 specs
  (GN032–037: REVEL ≤0.290) over the superseded v2 GN016 (≤0.15). COVERED.
- **RYR1** tie (GN012 Malignant-Hyperthermia <0.5 vs GN150/GN179 Congenital-
  Myopathies ≤0.15) resolved conservatively to **0.5** per the documented
  most-conservative tie rule. Acceptable; flagged because the disease-specific
  myopathy value (0.15) is stricter.

### Web diff (live cSpec API vs local TSV), top 3 material cases — all confirm the gaps
Endpoint: `https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/<GN>`

| gene | GN | web (live) REVEL BP4 | local TSV | verdict |
|------|----|----------------------|-----------|---------|
| F8 | GN071 | "REVEL score of 0.3 or below AND SpliceAI ≤0.1" → **≤0.3** | 0.1 | TSV wrong (MISMATCH) |
| SERPINC1 | GN084 | "REVEL score ≤0.30 and ... SpliceAI ≤0.1" → **≤0.30** | 0.1 | TSV wrong (MISMATCH) |
| GP9 | GN083 | "REVEL score ≤ 0.290 AND SpliceAI = 0" → **≤0.290** | *(blank)* | TSV incomplete (MISSING) |

F9 (GN080) shares F8's verbatim text and the same defect; not re-fetched (web
cap ≤3) but treated as confirmed by the identical GN071 wording.
