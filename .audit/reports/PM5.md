# PM5 VCEP Specialization Audit

Source of truth: `.audit/cspec_by_criterion/PM5.md` (282 Applicable entries across 123 populated specs).
App side: `pm5_grantham` / `pm5_excludes` / `pm5_max` / `pm5_lp` in `resources/shared/disease_prevalence.tsv`;
code `src/acmg_classifier/criteria/pm5_genes.py`, `pathogenic/pm5.py`, `grantham.py`.

App PM5 parameter semantics (from `pm5_genes.py`):
- `pm5_grantham` = Grantham gate operator: `ge` (candidate >= comparator) / `gt` (strictly >). Blank = no gate.
- `pm5_excludes` = criteria PM5 may not combine with (PM1 / PM1,PS1).
- `pm5_max` = strength ceiling (Strong / Moderate / Supporting). Blank => Moderate ceiling. PM5_Strong granted only when cap=Strong.
- `pm5_lp` = `no` => comparator must reach Pathogenic (LP-only comparator does not trigger PM5). Blank => LP comparator accepted at Supporting.

Scope of "specialization" = a Grantham/BLOSUM gate, a non-default tier ceiling (Strong / Supporting), a PM1/PS1 exclusion, an LP-acceptance flip, or a non-missense applicability flip (PTC/nonsense PM5). Generic Moderate-only with LP→Supporting is the app default and is not separately listed.

---

## MISSING (specialization in cspec, not represented in app)

| gene | GN | ver | cspec (grantham/tier/exclude) | app value | action |
|------|----|-----|-------------------------------|-----------|--------|
| PTEN | GN003 | 3.2.0 | **BLOSUM62 ≤** known variant (Moderate only; Strong/Supporting NA); P **or LP** comparator accepted | grantham=`` (none), cap=Moderate, lp=`no` | Add a BLOSUM62-based gate. App engine (`grantham.py`) only implements Grantham; **no BLOSUM62 matrix exists** — gate is silently dropped, so PTEN PM5 fires with no chemical-severity check. Web-confirmed (GN003). Also lp should be blank (cspec accepts LP comparator), but with no Supporting tier in cspec the lp=`no` is defensible; BLOSUM gap is the real miss. |
| RUNX1 / SCN1A / SCN2A / SCN3A / SCN8A / SCN1B (paralog & nonsense rules) | GN008 / GN067-070,076 | various | PM5 extended to **nonsense/frameshift** (RUNX1 downstream c.98) and **paralogous-gene** missense (SCN family NDD paralogues) | missense-only evaluator (`pm5.py` returns not-met for non-missense) | Out-of-model: app evaluates only same-gene missense. Nonsense-PM5 / paralogue-PM5 tiers are not implemented. Document as known scope limitation (affects RUNX1, SCN*, SCID point-system genes, ATM/PALB2 PTC, BRCA1/2 PTC). |
| ATM / PALB2 | GN020 / GN077 | 1.5.0 / 1.2.0 | PM5_Supporting applies **only to truncating/splice PTC** upstream of an anchor residue (no missense PM5) | cap=Supporting | Supporting ceiling captured, but the trigger is PTC, not missense; app's missense-only path will essentially never fire here. Flag as scope limitation, not a value error. |
| BRCA1 / BRCA2 | GN092 / GN097 | 1.2.0 | PM5_PTC: **protein-termination-codon** rule, per-exon Strong/Moderate/Supporting | cap=Strong | App treats as missense Strong gene; the actual PTC/per-exon logic is unmodeled. Scope limitation. |
| FOXN1/ADA/DCLRE1C/IL7R/JAK3/RAG1/RAG2/IL2RG (SCID) | GN113-129 | 2.x | **nonsense point-system** PM5 (4+→Strong, 2+→Mod, 1→Supp), PVS1 downgrade interactions | cap=Strong | Point-system + PVS1 interaction unmodeled; app applies generic missense tiers. Scope limitation. |

## MISMATCH (app value differs from cspec)

| gene | GN | ver | cspec (grantham/tier/exclude) | app value | action |
|------|----|-----|-------------------------------|-----------|--------|
| RYR1 | GN012 (MH, used) vs GN150/GN179 (CongMyo) | 2.0.0 / 2.0.0 | GN012: grantham **gt**, **Moderate** only, comparator must be P. GN150/179: **Strong** allowed, **no Grantham**, P only | grantham=`gt`, cap=**Strong**, lp=`` (blank) | Two VCEPs disagree. App blends them: keeps GN012's `gt` gate but raises cap to GN150 Strong. If RYR1 is curated under Malignant Hyperthermia (GN012), cap should be **Moderate**; if Congenital Myopathies (GN150/179), the `gt` gate should be removed. Also both specs require P comparator → `lp` should be `no`; currently blank (LP wrongly accepted at Supporting). Resolve per intended disease context and set lp=`no`. |
| ACVRL1 / ENG | GN135 / GN136 | 1.1.0 | "≥2 missense LP **or** P at codon → Strong; 1 LP **or** P → Moderate" — LP comparators **do** count | cap=Strong, lp=**`no`** | cspec explicitly counts likely-pathogenic comparators toward Moderate/Strong. App lp=`no` blocks LP-only comparators → under-calls. Set lp=blank for ACVRL1, ENG. |
| HNF1A / HNF4A / GCK | GN017 / GN085 / GN086 | 3.1.0 / 4.0.0 / 3.1.0 | Grantham `ge` gates the **Moderate** tier; **Strong** (2 P at residue) is **not** Grantham-gated | grantham=`ge` (applied to all tiers incl. Strong) | Minor over-gate: app applies the `ge` filter before counting distinct pathogenic hits, so a Strong-eligible residue could be withheld if the candidate fails Grantham vs one comparator. Low severity; note only. |

## VERSION & WEB

| gene | GN | ver | cspec (grantham/tier/exclude) | app value | action |
|------|----|-----|-------------------------------|-----------|--------|
| PIK3R1 | GN160 | 1.0.0 | Grantham **strictly higher** (GT); Moderate/Supporting | grantham=`gt`, cap=Moderate | **COVERED — web-confirmed** (GN160 says "must have a higher Grantham score"). `gt` is correct. |
| VHL | GN078 | 1.1.0 | Grantham **equal or larger** (GE); **Moderate only**, no Supporting; LP not accepted | grantham=`ge`, cap=Moderate, lp=`no` | **COVERED — web-confirmed** (GN078). All three params correct. |
| PTEN | GN003 | 3.2.0 | **BLOSUM62 ≤**, Moderate only | no gate | **MISMATCH — web-confirmed** (GN003 uses BLOSUM62, not Grantham). See MISSING table; engine cannot express a BLOSUM gate. |

Version note: RYR1 appears under two GN ids at the same version family; PTEN GN003 is v3.2.0 (current). No stale-version flags otherwise — `_spec_index.tsv` versions match the spec headers used here.

---

## Summary

- **Grantham/BLOSUM gates** — 16 genes carry `ge`/`gt` in the app (APC, CDKL5, DICER1, FOXG1, GALT, GCK, HNF1A, HNF4A, MECP2, PIK3CD `ge`; PIK3R1, RYR1 `gt`; RPGR, RS1, RUNX1, SLC9A6, TCF4, UBE3A, VHL). These match cspec, including the two strictly-`gt` cases (PIK3R1, RYR1-MH) — both verified.
- **1 hard MISSING gate**: **PTEN** requires BLOSUM62 (≤), which the engine cannot represent (Grantham-only). Top gap — PM5 currently fires for PTEN with no chemical-severity check.
- **3 value MISMATCHes**: **RYR1** (blended MH/CongMyo spec: cap vs gate conflict + lp should be `no`), **ACVRL1/ENG** (lp=`no` wrongly excludes LP comparators that cspec counts), **HNF1A/HNF4A/GCK** (Grantham over-gates the Strong tier — minor).
- **Scope limitations (MISSING, by design)**: nonsense/PTC PM5 (BRCA1/2, ATM, PALB2, RUNX1), paralogue PM5 (SCN1A/2A/3A/8A), and SCID nonsense point-systems (FOXN1, ADA, DCLRE1C, IL7R, JAK3, RAG1/2, IL2RG) are not modeled — the evaluator is missense-only. These are not value errors but unimplemented specialization classes.

**Top 3 actions**: (1) PTEN — add BLOSUM62 gate or explicitly document non-coverage; (2) ACVRL1/ENG — clear `pm5_lp` (LP comparators count per cspec); (3) RYR1 — resolve the MH-vs-CongMyo spec conflict and set `pm5_lp=no`.
