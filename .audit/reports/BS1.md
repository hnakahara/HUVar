# BS1 audit — cspec vs HUHVar app (`disease_prevalence.tsv`)

Source of truth: `.audit/cspec_by_criterion/BS1.md` (142 Applicable BS1 strength entries across populated specs).
App side: `resources/shared/disease_prevalence.tsv` columns `bs1_threshold` / `bs1_strength` / `bs1_exclude`, consumed by `src/acmg_classifier/criteria/benign/bs1.py`.
Parsing per `resources/clingen/README.md` (range → lower edge; sub-population rule; legacy 5% dropped; `*`/`Ter` equal; blank `bs1_strength` → Strong default).

Scope note: cspec BS1.md genes restricted to the app gene universe and to specs whose BS1 strength is *Applicable*. Genes whose spec is Pilot / In Prep / Submitted (no Applicable BS1) correctly have a blank `bs1_threshold` in the app (ACO2, CEP290, CHEK2, MUTYH, NF1, RHO, RYR2, SOD1, SDHB, PKD1/2, etc.) — not flagged.

## Summary counts

| Result | Cutoff | Strength | Exclusion |
|--------|--------|----------|-----------|
| COVERED | ~118 genes | all (Strong/VeryStrong) | MYOC ✓ |
| MISMATCH | 1 (RPGR) | 0 | 0 |
| MISSING | 1 disease-tier (ACTA1 AD) | 0 | 2 (RYR1 GN150, RYR1 GN179) |
| Second-tier not representable | — | 8 BS1_Supporting / 2 BS1_VeryStrong-pair lost | — |

Net: the cutoff/strength columns reproduce the cspec Strong-tier BS1 cutoff for essentially every covered gene. The real gaps are (a) one wrong cutoff (RPGR), (b) RYR1's two exclusion-variant lists never transcribed into `bs1_exclude`, and (c) an architectural inability to encode VCEPs' lower BS1 strength tiers (Supporting/Moderate) — the app supports only one BS1 strength per gene.

---

## MISMATCH

| gene | GN | ver | cspec (cutoff / strength / exclude) | app value | action |
|------|----|-----|-------------------------------------|-----------|--------|
| RPGR | GN106 | 1.0.0 | **≥8.3×10⁻⁵ = 0.000083** / Strong / — (in males) | bs1_threshold=**0.000005**, af_basis=males, Strong | **Fix cutoff → 0.000083.** App value is ~16.6× too low (likely a mis-keyed 5×10⁻⁶ vs the spec's most-frequent-pathogenic-allele 8.3×10⁻⁵). Web-confirmed (see below). Row tagged "manual override" in notes — the override is wrong. |

## MISSING

### Exclusion variants (recurrent disease alleles barred from BS1)

| gene | GN | ver | cspec exclude | app `bs1_exclude` | action |
|------|----|-----|---------------|-------------------|--------|
| RYR1 | GN150 | 2.0.0 | `p.Val4842Met` (c.14524G>A), `p.Arg109Trp` (c.325C>T) — AD Congenital Myopathies | *(empty)* | **Add to `bs1_exclude`.** Well-known pathogenic, above threshold. |
| RYR1 | GN179 | 2.0.0 | `p.Arg2241Ter` (c.6721C>T), `c.10348-6C>G` — AR Congenital Myopathies | *(empty)* | **Add to `bs1_exclude`.** Second is intronic (no p.) — current `_norm_pchange` only matches `p.` changes, so c.-level exclusions can't be encoded as-is (parser limitation). |

Only MYOC `p.Gln368Ter` (GN019) is currently populated; it is correct. No other Applicable BS1 spec in scope defines a bare-protein recurrent-allele exclusion (LGMD/DYSF and TP53/BRCA "exception lists" live in external supplementary files, not as inline protein changes, so they are out of scope for `bs1_exclude`).

### Disease-tier cutoff not selected (multi-spec conservative resolution)

| gene | GN (chosen) | other spec | cspec cutoffs | app value | action |
|------|-------------|-----------|---------------|-----------|--------|
| ACTA1 | GN169 (AR, 0.00025) | GN147 (AD, **0.00000781**) | AR 0.00025 / AD 7.81×10⁻⁶ | 0.00025 (AR, conservative) | Documented tie ("multiple specs kept conservative"). AD-disease ACTA1 variants get a 32× too-high BS1 cutoff (harder to fire BS1 → conservative/benign-safe). Acceptable per build policy; flag only. |

## VERSION & WEB

### Multi-spec / conservative-resolution flags (cutoff present, disease-ambiguous — not hard gaps)

| gene | chosen GN | competing specs | note |
|------|-----------|-----------------|------|
| RYR1 | GN179 (AR Cong. Myopathy 0.000697) | GN012 MH (0.0008), GN150 AD (4.86×10⁻⁶) | App uses 0.000697 (not the highest 0.0008). "Kept conservative" but GN012 MH cutoff 0.0008 is actually *higher* → for an MH-context variant the app fires BS1 slightly too easily. Minor; main issue is the missing exclusions above. |
| ACTA1 | GN169 AR | GN147 AD | see MISSING table |

### BS1 lower-strength tiers not representable (schema limit — one BS1 strength per gene)

cspec defines a **BS1_Supporting** (or extra) band the app cannot encode; only the Strong cutoff is kept:

| gene(s) | GN | cspec extra tier | app |
|---------|----|------------------|-----|
| ABCA4 | GN164 | BS1_Supporting >0.00163 (+ Strong 0.0163) | only Strong 0.0163 |
| PTEN | GN003 | BS1_Supporting 0.0000043–0.000043 (+ Strong from 0.000043) | only Strong 0.000043 |
| BRCA1 | GN092 | BS1_Supporting FAF>0.00002 (+ Strong 0.0001) | only Strong 0.0001 |
| BRCA2 | GN097 | BS1_Supporting FAF>0.00002 (+ Strong 0.0001) | only Strong 0.0001 |
| ACVRL1 | GN135 | BS1_Supporting >0.0008 (+ Strong >0.002) | only Strong 0.002 |
| ENG | GN136 | BS1_Supporting >0.0008 (+ Strong >0.002) | only Strong 0.002 |
| CDH23/COCH/GJB2/KCNQ4/MYO6/MYO7A/SLC26A4/TECTA/USH2A | GN005 | BS1_Supporting AR ≥0.0007 (+ Strong ≥0.003) | only Strong 0.003 |
| MYO15A / OTOF | GN023 | BS1_Supporting AR ≥0.0007 (+ VeryStrong 0.003) | only VeryStrong 0.003 |

These are correct as far as the strongest tier goes; the gap is a feature limitation (no `bs1_threshold_supporting`), not a data error. Recommend a follow-up if Supporting-tier BS1 is desired.

### Web diffs (cspec API JSON vs registry) — capped at 3

| gene | GN | JSON (BS1.md) | Web API | match? |
|------|----|--------------|---------|--------|
| RPGR | GN106 | "≥8.3×10⁻⁵ … in males … most frequent pathogenic allele" | identical verbatim; cutoff 8.3×10⁻⁵, Strong, in males | ✅ JSON==Web. App row is the outlier (0.000005), confirming the MISMATCH is an app error, not a stale slice. |
| RYR1 | GN150 | BS1 AD ≥0.00000486; exclusions p.Val4842Met, p.Arg109Trp | identical verbatim | ✅ JSON==Web. Confirms missing `bs1_exclude`. |
| ABCA4 | GN164 | BS1 Strong >0.0163, BS1_Supporting >0.00163 | identical (both tiers) | ✅ JSON==Web. Confirms Supporting tier dropped by schema, not by slice. |

No JSON-vs-Web discrepancies found in the sampled specs; the audit slice (`BS1.md`) faithfully reflects the live registry.
