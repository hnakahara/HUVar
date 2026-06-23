# PM1 Audit — VCEP hotspot specializations vs `resources/shared/pm1_hotspots.tsv`

Source of truth: `.audit/cspec_by_criterion/PM1.md` (101 Applicable PM1 strength entries).
App universe = genes present in `pm1_hotspots.tsv` (incl. `not_applicable` rows).
Scope: only specs that name a **concrete, codeable** hotspot (residues / codons / domains-with-coords / exon / non-default strength). Specs marked Applicable-default with no concrete region, or that point only to an external "Table"/"supplemental material" with no residues in the description, are flagged as **NOT-ACTIONABLE-FROM-SLICE** and not counted as gaps.

## Summary counts
- Concrete codeable cspec PM1 specializations reviewed: ~55 gene-entries
- **COVERED:** majority (PTEN, RUNX1, TP53, GAA, GP1BA, GP1BB, SERPINC1, RS1, BMPR2, F8, F9, DICER1, IDUA, GCK, HNF1A, HNF4A, ACVRL1, ENG, GUCY2D, RPE65, FOXN1, CTLA4, RAG1, RAG2, ABCD1, GALT(partial), all RASopathy genes, MYH7/MYBPC3/TNNI3/TNNT2, MECP2/CDKL5/FOXG1/TCF4/UBE3A, CYP1B1, PAH, OTC(partial), …)
- **MISSING (gene absent or concrete region not encoded):** 7
- **MISMATCH (present but region/residues/coords/strength diverge):** 2
- **NOT-ACTIONABLE-FROM-SLICE (default/Table-only, no gap):** SCN1A, SCN2A, SCN3A, SCN8A (GN067-070, "PM1 Table"); AKT3/MTOR/PIK3CA/PIK3R2 (GN018, "Table 4"); IL2RG (GN129, generic); IL7R (GN119, caveat only); GN001 (generic ACMG default).

## MISSING
| gene | GN | ver | cspec hotspot (region/strength) | app entry | action |
|------|----|-----|----------------------------------|-----------|--------|
| VHL | GN078 | 1.1.0 | Moderate: germline hotspots + key functional domains + somatic ≥10 in cancerhotspots; Supporting: somatic <10. "See table of Germline and Somatic Hotspots" | absent | Add VHL rows; resolve hotspot table (residues not in description — needs table ingestion) |
| KCNQ4 | GN005 | 2.0.0 | Moderate: pore-forming region aa **271-292** (NM_004700.4) | absent | Add KCNQ4 Moderate region 271-292 |
| PDHA1 | GN014 | 1.0.0 | Moderate: TPP-binding, αβ heterodimer, α2β2 heterotetramer, phospho-loop residues (long explicit residue list, ~70 positions) | absent | Add PDHA1 Moderate residue set |
| KCNQ1 | GN112 | 1.0.0 | Moderate: pore helix aa **300-320** (requires PM2_Supporting) | absent | Add KCNQ1 Moderate region 300-320 |
| JAK3 | GN121 | 2.3.0 | Moderate: JH2 residues **R651, C759** (R651W, C759R) | absent | Add JAK3 Moderate residues 651,759 |
| LDLR | GN013 | 1.2.0 | Moderate: missense in **exon 4**, or 60 conserved Cys residues (Supp Table 4); requires PM2 | absent | Add LDLR — needs exon-4 codon range + Cys list (Pilot) |
| FBN1 | GN022 | 1.0.0 | Strong: Cys in cbEGF-like domains; Moderate: Cys/Ca-binding/Gly motifs in EGF/TB/hybrid domains (motif-based) | absent | Motif/pattern rule — not expressible as residue ranges; defer or model as domain Cys set (Pilot) |

Non-coding-only specs (RMRP GN088 n.-32→n.4 insertion; HBB GN170 TATA/polyA; HBA2 GN173 polyA) name concrete regions but in nucleotide/promoter space, not protein residue positions the app's residue/range engine handles — flagged separately as **NON-CODEABLE-IN-CURRENT-SCHEMA**, not counted in the 7.

## MISMATCH
| gene | GN | ver | cspec hotspot (region/strength) | app entry | action |
|------|----|-----|----------------------------------|-----------|--------|
| OTC | GN156 | 1.0.0 | Moderate: 90,91,92,93,117,141,168,171,304,330 (CP-binding); 163,198,199,263,267,268 (ornithine); 302,303,304 (catalytic); 277,305,269 (conserved) — **~24 residues** | residues=`268` only | Replace with full residue set (currently encodes 1 of ~24) |
| GALT | GN158 | 1.0.0 | Moderate: active site **Phe171–Gln188** (contiguous range 171-188) | residues=`171,188`; regions empty | Encode as region `171-188` (currently only the two endpoint residues; 172-187 not awarded) |

## VERSION & WEB
| gene | GN | app→cspec ver | web diff result |
|------|----|---------------|-----------------|
| OTC | GN156 | matches v1.0.0 | Web JSON == slice verbatim; residue list confirmed (app encodes only 268 → MISMATCH stands) |
| VHL | GN078 | n/a (absent) | Web confirms Moderate+Supporting, hotspots live in external "Germline and Somatic Hotspots" table (not in JSON description) |
| KCNQ4 | GN005 | matches v2.0.0 | Web confirms aa 271-292 pore region verbatim; gene simply not in tsv |

No app entry implies an older spec version than cspec for the reviewed concrete genes (e.g. MYH7 167-931 matches GN002 v2.0.0 updated region; MECP2/CDKL5/FOXG1/TCF4/UBE3A app coords match the newer GN032-037 v6/v7 values). No version-staleness flags raised.

## Top 5 gaps
1. **OTC** (GN156) MISMATCH — encodes 1 of ~24 critical residues; 23 missing.
2. **VHL** (GN078) MISSING — entire gene absent; Moderate+Supporting hotspot table not ingested.
3. **KCNQ4** (GN005) MISSING — pore region 271-292 (NM_004700.4) absent.
4. **GALT** (GN158) MISMATCH — range 171-188 stored as two endpoints; interior residues unawarded.
5. **KCNQ1** (GN112) MISSING — pore helix 300-320 absent. (Also PDHA1 GN014 & JAK3 GN121 MISSING, concrete and directly codeable.)
