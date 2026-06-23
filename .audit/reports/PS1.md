# PS1 audit — cspec vs HUHVar app

Source of truth: `.audit/cspec_by_criterion/PS1.md` (247 Applicable strength entries across 123 populated specs).
App side: `resources/shared/disease_prevalence.tsv` columns `ps1` / `ps1_splice`;
code `src/acmg_classifier/criteria/ps1_genes.py` + `src/acmg_classifier/criteria/pathogenic/ps1.py`.

## App model (what can be encoded)
The app encodes only three PS1 specializations:
1. **Applicability flip** — `ps1 = not_applicable` (gene's VCEP declines PS1, e.g. CDH1).
2. **Splice-extension state** — `ps1_splice` ∈ {`""` missense-only, `canonical`, `noncanonical`}.
3. **Comparator-based strength** — PS1 fires Strong if a same-change comparator is ClinVar
   Pathogenic, else Moderate (`_ps1_strength`). This is *generic*, not per-gene.

The app **cannot** encode: a per-gene strength *cap* (e.g. PS1 max = Supporting), or
paralog/analogous-residue PS1 routes. These are structural gaps, flagged below.

## Methodology note on strength tiering
Most VCEPs list PS1_Strong (P comparator) + PS1_Moderate (LP comparator) — this two-tier
P/LP pattern is exactly what `_ps1_strength` reproduces, so it is **COVERED** generically and
not listed as a gap. Only *deviations* from that pattern (Supporting-only caps, paralog rules,
splice-mode mismatches, missing splice extensions) are actionable.

---

## MISSING (specialization in cspec, not encoded by app)

| gene | GN | ver | cspec (strength/splice) | app value | action |
|------|----|----|--------------------------|-----------|--------|
| ABCD1 | GN105 | 1.0.0 | PS1 applicable for **canonical splice-site** variant at same position, Strong (P) / Moderate (LP) | `ps1_splice=` (missense-only → blocks all splice PS1) | Set `ps1_splice=canonical`. App currently denies PS1 for any ABCD1 splice variant; cspec explicitly grants it at canonical positions. |
| RMRP | GN088 | 1.3.0 | PS1 **downgraded to Supporting** (only); same-nucleotide-position P/LP rule | `ps1=applicable`, `ps1_splice=` | No per-gene strength cap exists. App fires PS1 Strong/Moderate via comparator → over-weights. Needs a Supporting cap (new column/handling) or documented exception. |
| HBA2 | GN173 | 1.0.0 | PS1_Moderate = same change shown pathogenic **in a paralogue gene** | `ps1=applicable`, `ps1_splice=` | Paralog-PS1 route not modeled; app only matches same-gene hgvs_p. Out-of-model gap. |
| KCNQ1 | GN112 | 1.0.0 | PS1_Moderate via **paralogous KCNQ2** corresponding variant | `ps1=applicable`, `ps1_splice=canonical` | Paralog-PS1 route not modeled (splice handled OK). Out-of-model gap. |
| SCN1A/2A/3A/8A | GN067–070 | 2.x | PS1 Strong/Mod/Supp incl. **paralogous-gene identical-AA** (SCN1A/2A/3A/8A) routes | `ps1=applicable`, `ps1_splice=canonical` | Paralog-PS1 route not modeled (missense same-gene + splice handled OK). Out-of-model gap (affects 4 genes). |
| RASopathy genes (HRAS,KRAS,NRAS,RIT1,MRAS,RRAS2; BRAF/RAF1; SOS1/SOS2; MAP2K1/MAP2K2) | GN004/038–049/087/127 | 1–2.3 | PS1 applies to **analogous residue positions across paralogous group** | `ps1=applicable`, `ps1_splice=` | Cross-gene analogous-residue route not modeled; app matches only same-gene codon. Out-of-model gap (group rule). |

## MISMATCH (app value present but disagrees with cspec)

| gene | GN | ver | cspec (strength/splice) | app value | action |
|------|----|----|--------------------------|-----------|--------|
| GP1BA / GP1BB / GP9 / ITGA2B / ITGB3 | GN079/082/083/011 | 1.1–2.1 | PS1 missense-only; **no splice extension stated** (Platelet Disorders) — Strong(P)/Moderate(LP) | `ps1_splice=` ✓ | None — agrees (listed only to confirm Platelet genes correctly left missense-only). |
| RMRP | GN088 | 1.3.0 | Supporting-only + same-position rule | strength uncapped | (see MISSING) — also a strength mismatch: app can emit Strong, cspec caps at Supporting. |

No `ps1_splice` mode (`canonical` vs `noncanonical`) disagreements were found against cspec for
genes that carry a splice extension. Spot-checked the InSiGHT MMR genes (MLH1/MSH2/MSH6/PMS2 →
`noncanonical`, cspec = "same non-canonical splice nucleotide") and DICER1 (`noncanonical`,
cspec = "non-canonical intronic… at same nucleotide") — both **COVERED & correct**.

## VERSION & WEB

| gene | GN | ver | cspec (web) | app value | action |
|------|----|----|-------------|-----------|--------|
| ABCD1 | GN105 | 1.0.0 | Web confirms: "PS1 at strong level for 1 …canonical splice site variant (at the same position)"; Supporting N/A | `ps1_splice=` | Confirmed MISSING — set `canonical`. JSON slice == web. |
| RMRP | GN088 | 1.3.0 | Web confirms: "Downgraded to PS1_Supporting"; same-nucleotide-position P/LP rule | uncapped strength | Confirmed MISMATCH — needs Supporting cap. JSON slice == web. |
| DICER1 | GN024 | 1.4.0 | Web confirms: non-canonical intronic only (canonical is PVS1) | `ps1_splice=noncanonical` | Confirmed COVERED & correct. JSON slice == web. |

No version skew detected between `.audit` slice and live cspec API for the three checked GNs
(GN105 v1.0.0, GN088 v1.3.0, GN024 v1.4.0 all match). Web diff capped at 3 per instructions.

---

## Summary counts
- **MISSING: 2 in-model** (ABCD1 splice extension; RMRP Supporting cap) **+ ~3 out-of-model classes**
  (paralog PS1: HBA2, KCNQ1, 4×SCN genes; RASopathy analogous-residue group ~13 genes).
- **MISMATCH: 1** (RMRP strength — same item, dual-classified).
- **COVERED:** splice-mode assignments (`canonical`/`noncanonical`) verified correct for the
  InSiGHT MMR genes and DICER1; generic P/LP two-tier strength matches the dominant cspec pattern;
  CDH1 `not_applicable` flip correct.
- **Web:** all 3 fetched GNs (105/088/024) — live cspec == `.audit` slice, no discrepancy.

### Top gaps (priority order)
1. **ABCD1 `ps1_splice` empty → should be `canonical`** — app actively withholds a PS1 the VCEP grants. Single-cell TSV fix.
2. **Paralog/analogous-residue PS1 unmodeled** — structural; affects RASopathy group (~13 genes), 4 SCN genes, KCNQ1, HBA2. No column exists; the AA-match query is same-gene only.
3. **RMRP Supporting-only cap unmodeled** — app may over-weight to Strong; needs a per-gene strength cap mechanism.
