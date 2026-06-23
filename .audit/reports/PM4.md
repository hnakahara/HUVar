# PM4 Audit — VCEP specializations vs HUHVar app

Source of truth: `.audit/cspec_by_criterion/PM4.md` (175 Applicable strength entries across populated specs).
App side: `pm4` column in `resources/shared/disease_prevalence.tsv` + `src/acmg_classifier/criteria/pathogenic/pm4.py`.

## Structural finding (root cause)

The app can represent PM4 in only **two** states:

- `pm4` column ∈ {`applicable`, `not_applicable`, empty} — **binary applicability only**.
- `pm4.py` consumes *only* `not_applicable` (builds a decline set). Everything else fires at **default Moderate** for any in-frame indel / stop-loss outside a repeat region.

Therefore the app **cannot encode any PM4 specialization other than a full decline**:
- No `PM4_Supporting` strength.
- No `PM4_Strong` strength.
- No gene-specific region / residue restriction (allow-list domains or deny-list regions).
- No size-based strength split (1 aa = Supporting vs ≥2 aa = Moderate).

All non-default PM4 parameters below are consequently **MISSING** (not representable), unless noted as MISMATCH.

---

## MISMATCH (app value contradicts the spec it cites)

| gene | GN | ver | cspec (strength/region) | app value (`pm4`) | action |
|------|----|----|--------------------------|-------------------|--------|
| RYR1 | GN150 / GN179 | 2.0.0 | **Applicable** Strong/Moderate/Supporting ("no change — use as originally described") | `not_applicable` | FIX: RYR1 row cites `source_vcep=Congenital Myopathies` (GN179), which makes PM4 applicable. App wrongly declines PM4 for RYR1. Set `applicable` (+ ideally Strong tier). NB: a different RYR1 spec (GN012 Malignant Hyperthermia) is not the cited source. |

(No other `not_applicable` gene is contradicted: AKT3, APC, BRCA1, BRCA2, MLH1, MSH2, MSH6, MTOR, PALB2, PIK3CA, PIK3CD, PIK3R2, PMS2, TP53 do not appear as Applicable in PM4.md → declining is COVERED.)

---

## MISSING — non-default strength (PM4_Strong)

| gene | GN | ver | cspec (strength/region) | app value | action |
|------|----|----|--------------------------|-----------|--------|
| RUNX1 | GN008 | 3.1.0 | **PM4_Strong**: in-frame indel hitting any RHD residue R107,K110,A134,R162,R166,S167,R169,G170,K194,T196,D198,R201,R204; also stop-loss → extension. PM4_Mod same residues; PM4_Supporting other RHD aa 89-204 | `applicable` | Add Strong tier + residue/region model |
| MECP2 | GN036 / GN016 | 6.0.0 / 2.0.0 | **PM4_Strong** for stop-loss; PM4_Mod excludes Pro-rich p.381-405; PM4_Supporting <3 aa | `applicable` | Add Strong (stop-loss) + deny-region |
| UBE3A | GN037 / GN016 | 7.0.0 / 2.0.0 | **PM4_Strong** for stop-loss; PM4_Mod/Supporting tiers | `applicable` | Add Strong (stop-loss) |
| RPGR | GN106 | 1.0.0 | **PM4_Strong** stop-loss at aa 1153 (38-aa extension); PM4_Mod region-restricted exons 1-14 + ORF15 aa585-1078 | `applicable` (region col present for BP3/PM2, not PM4) | Add Strong + PM4 region allow-list |
| NEB | GN146 | 1.0.0 | **PM4_Strong** in-frame del (esp. exon 55, repetitive but pathogenic); PM4_Mod default; PM4_Supporting | `applicable` | Add Strong + repeat-exception |
| MYOC | GN019 | 2.1.0 | **PM4_Mod** within olfactomedin domain AA246-502 >10%; **PM4_Supporting** ≤10% (region-gated both tiers) | `applicable` | Add region gate + Supporting |

## MISSING — non-default strength (PM4_Supporting) and/or size split

| gene | GN | ver | cspec (strength/region) | app value | action |
|------|----|----|--------------------------|-----------|--------|
| GATM | GN025 | 2.0.0 | Mod for ≥2 aa, **Supporting** for single aa | `applicable` | Add size-based Supporting |
| GAMT | GN026 | 2.0.0 | Mod for stop-loss/≥2 aa, **Supporting** single aa | `applicable` | Add size-based Supporting |
| GAA | GN010 | 2.0.0 | Mod ≥2 aa <1 exon, **Supporting** single aa | `applicable` | Add size-based Supporting |
| IDUA | GN091 | 1.2.0 | Mod stop-loss/≥2 aa, **Supporting** single aa | `applicable` | Add size-based Supporting |
| HNF1A | GN017 | 3.1.0 | Mod >1 aa, **Supporting** single aa | `applicable` | Add size-based Supporting |
| HNF4A | GN085 | 4.0.0 | Mod >1 aa, **Supporting** single aa | `applicable` | Add size-based Supporting |
| GCK | GN086 | 3.1.0 | Mod >1 aa, **Supporting** single aa | `applicable` | Add size-based Supporting |
| DICER1 | GN024 | 1.4.0 | **Mod** only in RNase IIIb domain p.Y1682-S1846; **Supporting** outside (+repeat denylist) | `applicable` | Add domain gate + Supporting |
| RPE65 | GN120 | 1.0.0 | Mod ≥2 aa w/ conserved residue (PhyloP>2); **Supporting** 1 aa | `applicable` | Add PhyloP+size Supporting |
| GUCY2D | GN167 | 1.0.0 | Mod ≥2 aa conserved; **Supporting** 1 aa | `applicable` | Add PhyloP+size Supporting |
| AIPL1 | GN208 | 1.0.0 | Mod ≥2 aa conserved; **Supporting** 1 aa | `applicable` | Add PhyloP+size Supporting |
| ABCA4 | GN164 | 1.0.0 | Mod >1 aa / stop-loss / PhyloP≥7.367; **Supporting** 1 aa | `applicable` | Add Supporting |
| CDKL5 | GN034 / GN016 | 6.0.0 / 2.0.0 | Mod (deny C-terminus exons19-21/after p.904); **Supporting** <3 aa | `applicable` | Deny-region + Supporting |
| FOXG1 | GN035 / GN016 | 6.0.0 / 2.0.0 | Mod (deny His-rich p.37-57, Pro/Gln-rich p.58-86, Pro-rich p.105-112); **Supporting** <3 aa | `applicable` | Deny-region + Supporting |
| SLC9A6 | GN033 / GN016 | 6.0.0 / 2.0.0 | Mod; **Supporting** <3 aa | `applicable` | Add Supporting |
| TCF4 | GN032 / GN016 | 6.0.0 / 2.0.0 | Mod; **Supporting** <3 aa | `applicable` | Add Supporting |

## MISSING — gene-specific region / extra-requirement restriction (Moderate-only, non-default applicability logic)

| gene | GN | ver | cspec (strength/region) | app value | action |
|------|----|----|--------------------------|-----------|--------|
| PTEN | GN003 | 3.2.0 | Restricted to catalytic-motif indels (PM1) + protein extension | `applicable` | Add region gate |
| VHL | GN078 | 1.1.0 | Restricted to B/alpha domains (β63-155, α156-192, β2 193-204); deny indels before codon 54 | `applicable` (region col=14-48 used for BP3) | Add PM4 region gate |
| KCNQ1 | GN112 | 1.0.0 | Mod for any-size in-frame indel (location > size); mutually excl. PVS1/PP3 | `applicable` | Note location-not-size logic |
| CYP1B1 | GN104 | 1.0.0 | Stop-loss NOT a disease mechanism → PM4 not applicable to stop-loss | `applicable` | Add stop-loss exclusion |
| CDH1 | GN007 | 3.1.0 | Apply **only** to stop-loss variants | `applicable` | Add indel exclusion |
| ATM | GN020 | 1.5.0 | Use **only** for stop-loss variants | `applicable` | Add indel exclusion |
| OTC | GN156 | 1.0.0 | ≥1 aa but < whole exon; ≥1 exon → defer to PVS1 | `applicable` | Add size cap |
| LDLR | GN013 | 1.2.0 | Indels < 1 exon or whole-exon dup not in PVS1; must also meet PM2 | `applicable` | Add caveat |
| RS1 | GN126 | 1.0.0 | Indels < 1 exon, non-repeat, not in PVS1; must meet PM2 | `applicable` | Add caveat |
| RUNX1→FOXN1/ADA/DCLRE1C/IL7R/JAK3/RAG1/RAG2/IL2RG (SCID) | GN113/114/116/119/121/123/124/129 | various | PM4 on deletions requires deleted region contain known P/LP (Mod) or VUS (Supporting) variant | `applicable` | Add deletion-content requirement + Supporting |
| CTLA4 | GN122 | 1.0.0 | Mod only if ≥2 aa + conserved (PhyloP≥2) + SpliceAI<0.2; mut.excl PVS1/PP3 | `applicable` | Add conservation/splice gate |
| PIK3R1 | GN160 | 1.0.0 | Mod if ≥2 aa + PhyloP≥2 (indel) or stop-loss +≥2 aa C-term; polymorphic-region exclusion | `applicable` | Add conservation gate |
| FBN1 | GN022 | 1.0.0 | Caveat: cannot apply simultaneously with PVS1 (any strength) | `applicable` | Add mutual-exclusion |
| VWF | GN081 | 1.0.0 | N/A to type 2B variants (GoF) | `applicable` | Add subtype exclusion |
| HBB | GN170 | 1.0.0 | In-frame indels only (no stop-loss listed) | `applicable` | Minor: stop-loss scope |

(Cardiomyopathy genes GN002/095/098-103 MYH7/MYBPC3/TNNI3/TNNT2/TPM1/ACTC1/MYL2/MYL3: "may require downgrading to Supporting" — discretionary, no fixed param; treated as default Moderate = not actionable.)

---

## VERSION & WEB

| gene | GN | ver | cspec (strength/region) | app value | action |
|------|----|----|--------------------------|-----------|--------|
| RUNX1 | GN008 | 3.1.0 | PM4_Strong/Mod/Supporting RHD residues — **web JSON matches PM4.md verbatim** (residue list R107…R204 confirmed) | `applicable` | No JSON/web discrepancy |
| RPGR | GN106 | 1.0.0 | PM4_Strong stop-loss aa1153; PM4_Mod region exons1-14 + ORF15 585-1078 — **web matches** | `applicable` | No JSON/web discrepancy |
| MECP2 | GN036 | 6.0.0 | PM4_Strong stop-loss; Mod deny p.381-405; Supporting <3 aa — **web matches** | `applicable` | No JSON/web discrepancy |

Version notes: app rows track the **latest** GN per gene (e.g. MECP2→GN036 v6, UBE3A→GN037 v7, TCF4→GN032 v6, SLC9A6→GN033 v6, CDKL5→GN034 v6, FOXG1→GN035 v6), superseding the older GN016 v2 multi-gene Rett spec. No stale-version mismatches found in the `pm4` column itself; the columns simply lack the strength/region fields to carry these specs.
