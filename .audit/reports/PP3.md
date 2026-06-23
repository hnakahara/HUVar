# PP3 audit — gene-specific REVEL cutoffs (HUHVar vs ClinGen cSpec)

Scope: criterion **PP3** (computational evidence supporting deleterious effect).
App side encodes **REVEL-anchored** cutoffs only, in the
`revel_pp3_supporting / revel_pp3_moderate / revel_pp3_strong` columns of
`resources/shared/disease_prevalence.tsv` (consulted only when `--insilico-tool revel`;
blank → genome-wide Bergquist 2024 defaults 0.644 / 0.773 / 0.879(+3) / 0.932).
Code: `src/acmg_classifier/criteria/revel_genes.py`, `src/acmg_classifier/criteria/pathogenic/pp3.py`.
Source of truth: `.audit/cspec_by_criterion/PP3.md` (168 Applicable strength entries).

App semantics note: `_revel_pp3()` fires at the highest tier with `REVEL >= cutoff`.
A gene whose VCEP grants only Supporting fills just `*_supporting`, capping the gene
at Supporting (no promotion on a high score). README single-tier-capping +
most-specific-spec rules confirmed.

---

## MISSING — REVEL-anchored cutoffs (or tiers) the app does not represent

| gene | GN | ver | cspec REVEL cutoff (tier) | app value | action |
|------|----|-----|---------------------------|-----------|--------|
| AIPL1 | GN208 | 1.0.0 | Strong ≥0.093*; Moderate 0.774–0.092*; Supporting 0.644–0.773 | sup=0.644 only | **Add Moderate+Strong** — but see VERSION&WEB: cspec source has typos; do NOT transcribe 0.093/0.092 literally |
| GP9 | GN083 | 1.1.0 | Moderate ≥0.773; Supporting ≥0.644 | sup=0.644, mod=0.773, **bp4 blank** | PP3 side COVERED; flagged only because BP4 cols empty (out of PP3 scope) |
| AIPL1 (Moderate tier) | GN208 | 1.0.0 | Moderate ~0.774–0.932 (intended) | absent | **Add mod≈0.774** once typo resolved |

Net PP3 *gene*-level MISSING: only **AIPL1** has genuinely unrepresented higher tiers
(Moderate + Strong). Every other multi-tier REVEL spec is fully populated in the TSV.

---

## MISMATCH — app value differs from cspec REVEL cutoff

| gene | GN | ver | cspec REVEL cutoff (tier) | app value | action |
|------|----|-----|---------------------------|-----------|--------|
| HBB | GN170 | 1.0.0 | Supporting REVEL **> 0.8** (strict) | sup=`0.8` | Boundary semantics: app fires at REVEL == 0.8 (`>=`), spec needs strictly above. Set sup=0.801 or document `>=` approximation |
| HBA2 | GN173 | 1.0.0 | Supporting REVEL **> 0.8** (strict) | sup=`0.8` | Same `>` vs `>=` boundary mismatch as HBB |
| GAA | GN010 | 2.0.0 | Supporting REVEL **> 0.7** (strict) | sup=`0.7` | Minor `>` vs `>=`; app is 1 hundredth permissive at exactly 0.7 |
| RYR1 (Congenital Myop.) | GN150/GN179 | 2.0.0 | Supporting REVEL ≥0.7 | sup=`0.85` | **Strictness MISMATCH** — app uses 0.85 (the Malignant-Hyperthermia GN012 cutoff `>0.85`). Multi-spec tie kept conservative (higher = fewer PP3 calls). cspec Congenital-Myopathy value is 0.7. Confirm intended disease |
| SLC9A6 | GN033 | 6.0.0 | Supporting REVEL ≥0.664 | sup=`0.664` | COVERED (listed to show the unusual 0.664, not 0.644, is correctly transcribed) |

Note: the `>` vs `>=` rows are *pervasive* (most VCEPs write `>0.7`/`>0.75`/`>0.8`,
the app stores the bare number and compares with `>=`). This is a systematic 1-score-unit
over-permissiveness at the exact boundary, not a per-gene data error. Flagged here as the
material instances; if the SVI intends strict `>`, it affects ~25 Supporting rows.

---

## NON-REVEL specializations the app cannot represent (out of scope for REVEL columns)

These specs anchor PP3 on a numeric cutoff from a **non-REVEL** predictor. They are
correctly NOT in the `revel_pp3_*` columns, but the app cannot reproduce the VCEP rule
for these tools (SpliceAI is handled generically at Moderate via `_spliceai_pp3`, but the
gene-specific SpliceAI thresholds — 0.5, 0.38, 0.8, 0.3 — are not encoded; BayesDel / CADD /
AlphaMissense / aGVGD / HCI-prior / MAPP have no PP3 path at all).

| gene(s) | GN | predictor + cutoff | notes |
|---------|----|--------------------|-------|
| TP53 | GN009 | BayesDel ≥0.16 (+aGVGD class); SpliceAI ≥0.2 | no REVEL; pure BayesDel/aGVGD |
| BRCA1 | GN092 | BayesDel no-AF ≥0.28; SpliceAI ≥0.2 | domain-gated BayesDel |
| BRCA2 | GN097 | BayesDel no-AF ≥0.30; SpliceAI ≥0.2 | domain-gated BayesDel |
| RUNX1 | GN008 | REVEL ≥0.88 **or** SpliceAI ≥0.38 | REVEL tier IS encoded (sup=0.88); the 0.38 SpliceAI is the non-REVEL part |
| CTLA4 | GN122 | REVEL ≥0.75 **AND** CADD ≥20 | app encodes REVEL 0.75 but cannot enforce the AND-CADD requirement |
| PIK3CD | GN141 | REVEL ≥0.644 **AND** CADD ≥25.3 | REVEL encoded; AND-CADD not enforceable |
| PIK3R1 | GN160 | REVEL ≥0.644 **AND** CADD ≥26.0 | REVEL encoded; AND-CADD not enforceable |
| BMPR2 | GN125 | 2-of-3 {REVEL≥0.644, AlphaMissense≥0.792, CADD≥25.3} | REVEL encoded; 2-of-3 logic not enforceable |
| ABCA4 | GN164 | REVEL (mod>0.772, sup 0.644–0.772) **or** CADD 25.3/28.1; SpliceAI 0.2/0.8 | REVEL encoded; CADD/SpliceAI-band not |
| HBB / HBA2 | GN170/GN173 | CADD>23.5 fallback; SpliceAI>0.3 | REVEL>0.8 encoded; CADD-fallback not |
| MLH1 | GN115 | HCI prior P >0.81 (mod) / >0.68 (sup); SpliceAI ≥0.2 | **no REVEL** — HCI-prior tool; nothing encodable |
| MSH2 | GN137 | HCI prior P >0.81 / >0.68; SpliceAI ≥0.2 | no REVEL |
| MSH6 | GN138 | HCI prior P >0.81 / >0.68; SpliceAI ≥0.2 | no REVEL |
| PMS2 | GN139 | MAPP/PP2 prior P >0.81 / >0.68; SpliceAI ≥0.2 | no REVEL |
| ABCD1 | GN105 | REVEL >0.85; SpliceAI ≥0.5 | REVEL encoded (0.85); SpliceAI 0.5 not |
| GP1BA/GP1BB/GP9 | GN079/82/83 | REVEL tiers + SpliceAI ≥0.5 | REVEL encoded; SpliceAI 0.5 not |
| SCN1A/2A/3A/8A/1B | GN067-070,076 | "REVEL, follow ClinGen rec." (no number) | blank → Bergquist default applies (correct per README) |
| APC | GN089 | splice predictors only, "do not use REVEL for missense" | no numeric REVEL |
| PALB2 | GN077 | "Missense: do not use"; SpliceAI ≥0.2 | correctly blank |
| ADA/DCLRE1C/IL7R/JAK3/RAG1/RAG2/IL2RG | GN114/116/119/121/123/124/129 | SpliceAI ≥0.2 only, "do not apply to missense" | SCID panel — no REVEL by design |
| RMRP | GN088 | RNAsnp p<0.1 | exotic tool; not encodable |

Specs using a **non-REVEL predictor with a numeric cutoff** for PP3:
**~30 GN ids** (SpliceAI-gene-specific, BayesDel, CADD, AlphaMissense, aGVGD, HCI/MAPP prior,
RNAsnp). Of these, **6 panels have NO REVEL path at all** (TP53, BRCA1, BRCA2, the 4 InSiGHT
HCI/MAPP genes MLH1/MSH2/MSH6/PMS2, the SCID SpliceAI-only genes) — the app's PP3 cannot
reproduce their VCEP rule by design.

---

## VERSION & WEB

Web diff (≤3 cases) via `https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/<GN>`:

| gene | GN | JSON (PP3.md) | Web (live API) | verdict |
|------|----|---------------|----------------|---------|
| AIPL1 | GN208 | Strong ≥0.093; Mod 0.774–0.092; Sup 0.644–0.773 | **identical** — Strong "≥0.093", Mod "0.774–0.092", Sup "0.644–0.773" | Web confirms the JSON; the **typos originate in cSpec itself** (0.093/0.092 are clearly meant to be 0.932/0.932 reading "0.774–0.932"). App correctly stored only Sup=0.644 and skipped the corrupt tiers (monotonicity guard). **Do not transcribe the literal typos.** True intended: Strong ≥0.932, Mod 0.774–0.931, Sup 0.644–0.773 |
| HBB | GN170 | REVEL >0.8; SpliceAI >0.3; CADD>23.5 | **identical** — "REVEL score > 0.8 OR SpliceAI > 0.3 … CADD PHRED > 23.5" | Web confirms strict `>0.8`. App `sup=0.8` + `>=` is the boundary mismatch above |
| RPGR | GN106 | Strong 0.932; Mod 0.773–0.931; Sup 0.644–0.772 | **identical** — Strong "above 0.932", Mod "0.773–0.931", Sup "0.644–0.772" | COVERED — TSV (0.644/0.773/0.932) matches exactly |

Version flags: all three live-checked specs match their JSON snapshot; no stale-version
drift detected. RS1 (GN126) TSV stores strong=`0.931` (cspec "above 0.931") vs RPGR's
`0.932` — both faithful to their respective specs (intentional 1-unit VCEP difference,
not an app error).

---

## SUMMARY

- **REVEL-anchored PP3 specs cross-checked:** ~55 GN ids carry a numeric REVEL PP3 cutoff;
  all map to TSV rows.
- **COVERED:** the large majority — RASopathy (0.7), Cardiomyopathy (0.7), Hearing Loss (0.7),
  LGMD (0.7), Glaucoma/CYP1B1/MYOC (0.644/0.773/0.932), PAH/GAMT/OTC (0.644/0.773/0.932),
  RPGR/RS1, Rett-Angelman (0.644, SLC9A6 0.664), VHL (0.664), MYOC, etc. Multi-tier specs
  are correctly populated.
- **MISSING (genuine):** **1 gene — AIPL1 (GN208)** lacks Moderate+Strong tiers, blocked by
  corrupt source numbers (0.093/0.092). Action: add Mod≈0.774, Strong≈0.932 with a typo note,
  or leave capped at Supporting pending VCEP correction.
- **MISMATCH:** **RYR1** (app 0.85 vs Congenital-Myopathy 0.7 — multi-spec conservative tie,
  confirm disease) is the one substantive value mismatch. Plus a **systematic `>` vs `>=`
  boundary** over-permissiveness affecting ~25 Supporting rows (HBB/HBA2 >0.8, GAA >0.7, …) —
  app fires at the exact cutoff, several VCEPs write strict `>`.
- **NON-REVEL (out of scope):** **~30 GN specs** anchor PP3 on SpliceAI(gene-specific)/BayesDel/
  CADD/AlphaMissense/aGVGD/HCI-prior/MAPP/RNAsnp. **6 panels have no REVEL path at all**
  (TP53, BRCA1, BRCA2, InSiGHT MLH1/MSH2/MSH6/PMS2, SCID SpliceAI-only genes) and cannot be
  represented in the REVEL columns by design.

**Top gaps:** (1) AIPL1 missing higher tiers (source-typo-blocked); (2) RYR1 0.85-vs-0.7
disease-tie; (3) systematic strict-`>` boundary semantics; (4) AND/2-of-3 multi-tool logic
(CTLA4, PIK3CD, PIK3R1, BMPR2) un-enforceable — REVEL tier stored but co-tool requirement lost.
