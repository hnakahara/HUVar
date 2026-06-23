# PVS1 VCEP-specialization audit

Source of truth: `.audit/cspec_by_criterion/PVS1.md` (317 Applicable strength entries across 123 populated specs) + raw `resources/clingen/cspec_json/GN*.json` for PVS1 = Not Applicable specs.

App side checked:
- `resources/shared/disease_prevalence.tsv` column `pvs1` (suppression flag).
- `src/acmg_classifier/pvs1/vcep_pvs1.py` `_SPECS` (61 gene trees) + `apc.py` (APC).
- `resources/shared/vcep_pvs1_splice_exons.tsv` (per-exon canonical-splice overrides: DICER1, CDKL5, HNF4A, PAH, ABCD1, GAA, GAMT, IDUA).

**Scope of "parameterized":** a cspec PVS1 rule that gives concrete codon bands, NMD-escape boundaries, named in-frame/critical exons, or an explicit start-loss / N/A override. Generic "Applicable, follow the SVI/Tayoun decision tree" boilerplate (e.g. SCN1A, the LGMD genes, the Hearing-Loss flowchart genes, SLC6A8, the GN014 mitochondrial flowsheet genes, BMPR2, LDLR, SERPINC1, ITGA2B/ITGB3) is **not** counted as a gap.

**PVS1 = Not Applicable (LoF not the mechanism):** 33 specs scanned from raw JSON (MYH7, all RASopathy genes incl. MRAS/PPP1CB/RRAS2, AKT3/MTOR/PIK3CA/PIK3R2, MYOC, PIK3CD, all 6 Cardiomyopathy thin-filament/regulatory genes TNNI3/TNNT2/TPM1/ACTC1/MYL2/MYL3, ACTA1, DNM2, RYR1, VWF, RMRP). **All 33 are correctly `not_applicable` in `disease_prevalence.tsv`. No suppression gaps.**

---

## MISSING — cspec defines a parameterized PVS1 specialization, app has nothing

| gene | GN | spec ver | cspec value (strength) | app value | recommended action |
|------|----|----------|------------------------|-----------|--------------------|
| PIK3R1 | GN160 | 1.0.0 | NMD c.917–c.1890→PVS1; codons 631–645 / cSH2 646–718→Strong; 719–724→Moderate; splice skip exon 8/10/12→PVS1, in-frame exon 9/11(nSH2) /13/14(iSH2) /15/16(cSH2)→Strong; **start-loss N/A**; exon 1–7 / AR-agamma-only region (c.4–c.916) N/A | `pvs1=applicable`, no `_SPECS` entry → runs **generic tree only** | Add `_GeneSpec` with trunc_bands + exon-aware splice + `start_lost=NOT_MET`. Multi-transcript / disease-scope (MONDO:1060136) exclusion likely cannot be fully modelled; encode at least the codon bands + start-loss N/A. |
| GATM | GN025 | 2.0.0 | NMD-escape last exon 9 / last 50 nt exon 8 (3' of c.1109)→Strong(>10%)/Moderate(<10%); **initiator codon→PVS1_Moderate** (next Met p.130) | `pvs1=applicable`, no `_SPECS` entry → generic tree only | Add `_GeneSpec(trunc_nmd=(VERY_STRONG, None), start_lost=MODERATE, splice=VERY_STRONG, deletion=VERY_STRONG)`. Mirrors the already-encoded GAMT (same VCEP). |
| GALT | GN158 | 1.0.0 | In-frame exons 6, 7, 9 (splice/skip handled at reduced strength per attached decision tree) | `pvs1=applicable`, no `_SPECS` and no row in `vcep_pvs1_splice_exons.tsv` | Add splice-exon overrides for GALT exons 6/7/9 (Strong/Moderate per %-protein from `vcep_pvs1_exons.tsv`), plus a baseline `_GeneSpec`. |
| OTOF | GN023 | 1.0.0 | PVS1 null default, **exon 46 (c.5841–c.5994) is an exception** requiring scrutiny (high-freq / non-pathogenic LoF region) | `pvs1=applicable`, no `_SPECS` entry | Add exon-46 carve-out (suppress / downgrade) for OTOF; default PVS1 elsewhere. |
| MYO15A | GN023 | 1.0.0 | PVS1 null default, **exon 8 (c.4033–4038) and exon 26 (c.5911–c.5964) are exceptions** | `pvs1=applicable`, no `_SPECS` entry | Add exon 8 / exon 26 carve-outs for MYO15A. |
| RYR1 | GN179 | 2.0.0 | In-frame deletion / in-frame exon-skip in pore/TM region (aa4800–4950, exons 100–103)→PVS1_Strong | `pvs1=not_applicable` (Congenital Myopathies suppresses PVS1 generally) | Low priority: the app suppresses RYR1 PVS1 entirely (consistent with GN012/GN150/GN179 marking the generic null rule N/A). The Strong in-frame-pore *up-call* is the one residual case lost by blanket suppression. Document, or add a narrow exon-100–103 in-frame override. |

Notes on **partial coverage** (encoded but in-frame splice exon exceptions only approximated): DICER1, GAA, TP53, HNF4A, GCK, CDKL5, RUNX1, NEB are in `_SPECS` with a flat `splice=VERY_STRONG`; only DICER1/CDKL5/HNF4A/PAH/ABCD1/GAA/GAMT/IDUA have per-exon rows in `vcep_pvs1_splice_exons.tsv`. NEB (in-frame exons 3–180,182 → Strong; exon 55 critical) and GCK (exons 8/9 >10%→PVS1, exons 4/5 active-site→PVS1) lack exon rows — splice variants skipping in-frame exons in these genes are over-called at Very Strong vs the cspec Strong/Moderate. These are refinements, not full MISSING, and are explicitly listed as out-of-scope in the module docstrings.

---

## MISMATCH — encoded but diverges from cspec

| gene | GN | spec ver | cspec value (strength) | app value | recommended action |
|------|----|----------|------------------------|-----------|--------------------|
| — | — | — | — | — | No hard value mismatches found. All `_SPECS` codon bands / start-loss / splice defaults reviewed against the verbatim cspec text agree (including the documented transcript-offset remaps for SLC9A6 +10, UBE3A +20, MECP2 e1=e2+12, and MECP2 start-loss = N/A). |

The only divergences are the deliberately-deferred refinements noted above (in-frame splice-exon strength for genes whose splice default is left flat) and the deferred exon-level dup "proven/presumed in tandem" and RNA-evidence modifiers — all documented as out-of-scope in `vcep_pvs1.py` / `apc.py` docstrings, not silent bugs.

---

## VERSION & WEB

| gene | GN | spec ver (JSON) | app implied ver | web check | result |
|------|----|-----------------|-----------------|-----------|--------|
| PIK3R1 | GN160 | 1.0.0 | n/a (not encoded) | fetched `…/id/GN160` | Web v1.0.0 = JSON; codon bands (c.917–1890, 631–645, 646–718, 719–724), exon 8/10/12 & cSH2/nSH2/iSH2 rules, start-loss N/A all match. No JSON-vs-Web discrepancy. |
| GATM | GN025 | 2.0.0 | n/a (not encoded) | fetched `…/id/GN025` | Web v2.0.0 = JSON; NMD-escape (exon 9 / exon 8 3' of c.1109), initiator Met130→Moderate, >10%/<10% split all match. No discrepancy. |
| OTOF / MYO15A | GN023 | 1.0.0 | n/a (not encoded) | fetched `…/id/GN023` | Web v1.0.0 = JSON; OTOF exon 46 (c.5841–5994), MYO15A exon 8 (c.4033–4038) & exon 26 (c.5911–5964) exceptions match. No discrepancy. |

No app `source_vcep` / `notes` value implies an older spec version than the JSON for the encoded genes (the RASopathy and Rett/Angelman re-versioned specs — e.g. GN032–037 v6/v7, GN038–049 v2.3.0 — are the ones the app cites, and the encoded codon bands match those later versions). Web fetch hard cap (3) respected.
