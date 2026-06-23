# BS2 Audit — VCEP gene-specific specializations vs HUHVar app

Source of truth: `.audit/cspec_by_criterion/BS2.md` (178 Applicable strength entries across populated specs).
App side: `resources/shared/disease_prevalence.tsv` columns `bs2`, `bs2_count`, `bs2_female_only`, `bs2_hom_only`; code `src/acmg_classifier/criteria/bs2_genes.py` + `criteria/benign/bs2.py`.

## How the app models BS2
- `bs2=applicable` + `bs2_count` (+ `bs2_hom_only`/`bs2_female_only`) → gnomAD-count path. `bs2_count` overrides ALL inheritance-mode thresholds (defaults hom=2, hemi=2, het=3).
- `bs2=not_applicable` → gnomAD path barred; falls back to ClinVar expert-panel (>=3-star) BS2 evidence. Used both for VCEPs that bar population data (RASopathy) AND for VCEPs whose BS2 needs internal cohorts gnomAD can't supply (CDH1, TP53, SERPINC1, etc.).
- **Structural limit:** there is NO `bs2_strength` column. The gnomAD-count path always asserts BS2 at **Strong**. VCEPs that tier Strong/Moderate/Supporting on count cannot be represented; the TSV picks the lowest-tier count, so BS2 fires Strong at the most lenient (Supporting) threshold.

## Counts
- Specs with a concrete gene-specific BS2 parameter (count / zygosity / sex / strength tiers): ~30 spec families.
- COVERED (parameter faithfully represented): ~12
- MISMATCH: 8
- MISSING (gnomAD-usable threshold not encoded; barred): 2
- VERSION/structural notes: several (strength-tier loss; female_only dead column)

---

## MISSING — gnomAD-usable count threshold present in spec but app bars the gnomAD path

| gene | GN | ver | cspec (count/zygosity/sex) | app value | action |
|------|----|----|----------------------------|-----------|--------|
| DCLRE1C | GN116 | 2.2.0 | >=3 hom Strong / >=1 hom Supporting (gnomAD-eligible) | bs2=not_applicable | Set bs2=applicable, bs2_count=1, bs2_hom_only=1, inh=AR (mirror IL7R/JAK3/RAG1/RAG2 SCID treatment) |
| IL2RG | GN129 | 2.2.0 | >=3 hemizygotes Strong / >=2 hemi Supporting (explicit "in gnomAD") | bs2=not_applicable | Set bs2=applicable, bs2_count=2, inh=XL (hemi threshold). Currently no BS2 unless ClinVar EP review exists |

(Both SCID-panel genes are explicitly gnomAD-count based, so `not_applicable` wrongly routes them to the ClinVar-only fallback. Sibling SCID genes RMRP/IL7R/JAK3/RAG1/RAG2 are correctly `applicable`.)

---

## MISMATCH — parameter encoded but value/semantics diverge from spec

| gene | GN | ver | cspec (count/zygosity/sex) | app value | action |
|------|----|----|----------------------------|-----------|--------|
| GUCY2D | GN167 | 1.0.0 | gnomAD Strong >=6 hom; Supporting >=3 hom (gnomAD v4.1.0+) | bs2_count=3 (fires Strong at 3) | Raise count to 6 (gnomAD Strong threshold); 3 is the Supporting count and over-fires Strong |
| AIPL1 | GN208 | 1.0.0 | gnomAD Strong >=6 hom; Supporting >=3 hom (gnomAD v4.1.0+) | bs2_count=3 (fires Strong at 3) | Raise count to 6 (same as GUCY2D) |
| BMPR2 | GN125 | 2.0.0 | >=3 hom Strong / >=2 hom Mod / >=1 hom Supporting | bs2_count=1, hom_only=1 (fires Strong at 1) | Count=1 is the Supporting tier; gnomAD path always asserts Strong → over-strong. Needs strength tiering or count=3 for Strong |
| RMRP | GN088 | 1.3.0 | >=3 hom Strong / >=2 hom Supporting | bs2_count=2, hom_only=1 (fires Strong at 2) | count=2 is the Supporting tier; over-fires Strong. Set count=3 for Strong or add tiering |
| IL7R / JAK3 / RAG1 / RAG2 | GN119/121/123/124 | 2.2.0/2.3.0 | >=3 hom Strong / >=1 hom Supporting | bs2_count=1, hom_only=1 (fires Strong at 1) | count=1 is Supporting tier; over-fires Strong. Set count=3 for Strong or add tiering |
| OTC | GN156 | 1.0.0 | >5 female homozygotes OR 5 male hemizygotes (XL) | bs2_count=6, hom_only=1, XL | hom_only forces nhomalt; spec also allows 5 male hemizygotes. count=6 (>5) is right for hom but the hemi-path (5 hemi) is unreachable under hom_only. Verify intended; consider hemi handling |
| ABCD1 | GN105 | 1.0.0 | >10 hemizygous healthy adult (>40y) MALES (XL) | bs2_count=11, XL (no sex/age) | count=11 correct for hemi; age>40 + male restriction not encoded (gnomAD lacks age). Acceptable approximation — note only |

---

## VERSION & WEB diff (3 most material cases verified live)

Endpoint: `https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/<GN>`

| gene | GN | ver (JSON / web) | web says (verbatim) | app value | action |
|------|----|------------------|---------------------|-----------|--------|
| DCLRE1C | GN116 | 2.2.0 / 2.2.0 ✓ | "BS2_Strong ... at least 3 homozygotes"; "BS2_Supporting ... at least 1 homozygote" | bs2=not_applicable | Web confirms gnomAD-eligible hom-count rule → flip to applicable (see MISSING) |
| IL2RG | GN129 | 2.2.0 / 2.2.0 ✓ | "BS2_Strong: Observed in >=3 hemizygotes in gnomAD"; "Supporting ... at least 2 hemizygotes in gnomAD" | bs2=not_applicable | Web explicitly says "in gnomAD" → flip to applicable, count=2 hemi (see MISSING) |
| GUCY2D | GN167 | 1.0.0 / 1.0.0 ✓ | Strong: ">=6 homozygotes in gnomAD v.4.1.0+"; Supporting: ">=3 homozygotes in gnomAD" | bs2_count=3 | Web confirms Strong gnomAD threshold is 6, not 3 → raise count to 6 (see MISMATCH) |

No version skew detected for the 3 checked specs (JSON versions match live cspec).

---

## Notes / lower priority
- **Dead column:** `bs2_female_only` is implemented in code (TP53/DICER1 rationale) but set on NO gene in the TSV (TP53 GN009 and DICER1 GN024 are `not_applicable` → ClinVar fallback, bypassing the female-count path). The female-restricted specs (TP53 >=8 females >=60y; DICER1 40+ females tumor-free) are handled via `not_applicable` instead. Harmless but the `female_only` evaluator branch is unreachable.
- **Strength-tier loss (systemic):** BMPR2, RMRP, the SCID genes, GUCY2D, AIPL1, TP53, PALB2/BRCA1/BRCA2 (point-based >=4/2/1) all define Strong/Moderate/Supporting tiers. The single-threshold + always-Strong model cannot represent these; current counts are set to the lowest tier, biasing toward asserting BS2 at Strong. A `bs2_strength` (or per-tier count) column would be the principled fix.
- **Generic "Applicable, no count" specs** (SCN1A/SCN2A/SCN3A/SCN8A/SCN1B epilepsy; congenital myopathy NEB/ACTA1/DNM2/MTM1/RYR1 "no change"; GALT/GAMT/PAH/POLG/etc. enzyme-activity rules) carry no gnomAD count; mapped to applicable + default thresholds or not_applicable as appropriate — acceptable.
