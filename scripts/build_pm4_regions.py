#!/usr/bin/env python3
"""Build pm4_regions.tsv — per-gene PM4 region / strength rules from ClinGen VCEPs.

PM4 ("protein-length change") is Moderate by default, but several VCEPs make it
gene-specific: a residue/region awards Strong, restricts where PM4 applies at
all (allow-list with an N/A default), denies repeat/low-complexity regions, or
sets a stop-loss-specific strength. These rules are transcribed verbatim from
the cspec descriptions (residue lists / domain coordinates / decision text) into
a small table the PM4 evaluator consumes — analogous to pm1_hotspots.tsv.

Row kinds (the ``strength`` column):
  * ``strong`` / ``moderate`` / ``supporting`` — an in-frame indel impacting a
    residue/region gets that strength (strongest match wins).
  * ``deny`` — an in-frame indel in this region is withheld (repeat / Pro-rich /
    low-complexity regions the VCEP excludes from PM4_Moderate).
  * ``region_default`` — the strength for an in-frame indel matching no region:
    ``not_met`` (allow-list genes), ``supporting`` (DICER1), or ``moderate``.
  * ``stoploss`` — strength for a stop-loss variant (or ``not_applicable``).
  * ``conserved_phylop`` — PhyloP cutoff; an in-frame indel fires PM4 only when
    the position is conserved (PhyloP > cutoff). Skipped if PhyloP unavailable.
  * ``deletion_content`` (``yes``) — an in-frame DELETION fires PM4 only if the
    deleted protein range contains a known ClinVar P/LP (→Moderate) or VUS
    (→Supporting) variant (the SCID panel rule).
  * ``excludes`` — criteria PM4 is mutually exclusive with (``PVS1`` / ``PVS1,PP3``);
    the registry suppresses PM4 when any listed criterion also fired.
"""
from __future__ import annotations

import argparse
import csv

# gene -> list of (strength, regions[list of (a,b)], residues[list of int]).
# Coordinates are 1-based protein positions on the VCEP's clinical transcript.
_CURATED: dict[str, list[tuple[str, list, list]]] = {
    # RUNX1 (GN008, NM_001754): PM4_Strong on the RHD residues R107,K110,A134,
    # R162,R166,S167,R169,G170,K194,T196,D198,R201,R204; PM4_Supporting on the
    # other RHD residues 89-204; outside the RHD → N/A. (Stop-loss → PVS1
    # extension, handled by the PVS1 module, so no PM4 stoploss row here.)
    "RUNX1": [
        ("strong", [], [107, 110, 134, 162, 166, 167, 169, 170, 194, 196, 198, 201, 204]),
        ("supporting", [(89, 204)], []),
        ("region_default", "not_met", None),
    ],
    # MYOC (GN019): in-frame del/ins within the conserved olfactomedin domain
    # AA 246-502 → Moderate; outside → N/A.
    "MYOC": [
        ("moderate", [(246, 502)], []),
        ("region_default", "not_met", None),
    ],
    # DICER1 (GN024): in-frame indel in the RNase IIIb domain p.Y1682-S1846 →
    # Moderate; outside → Supporting; repeat regions p.D606-609 / p.E1418-1420 /
    # p.E1422-1425 → denied.
    "DICER1": [
        ("moderate", [(1682, 1846)], []),
        ("deny", [(606, 609), (1418, 1420), (1422, 1425)], []),
        ("region_default", "supporting", None),
    ],
    # MECP2 (GN036): PM4_Moderate excludes the Pro-rich p.381-405 region (size
    # <3 aa → Supporting is handled by pm4_supporting_max_aa).
    "MECP2": [
        ("deny", [(381, 405)], []),
        ("region_default", "moderate", None),
    ],
    # CDKL5 (GN034): PM4_Moderate excludes the C-terminus after p.904 (exons
    # 19-21); the gene's MANE isoform is 960 aa.
    "CDKL5": [
        ("deny", [(905, 960)], []),
        ("region_default", "moderate", None),
    ],
    # FOXG1 (GN035): PM4_Moderate excludes the His-rich p.37-57, Pro/Gln-rich
    # p.58-86 (merged 37-86) and Pro-rich p.105-112 low-complexity regions.
    "FOXG1": [
        ("deny", [(37, 86), (105, 112)], []),
        ("region_default", "moderate", None),
    ],
    # CDH1 (GN007) / ATM (GN020): PM4 applies ONLY to stop-loss variants → an
    # in-frame indel is N/A; a stop-loss is Moderate (the default strength).
    "CDH1": [
        ("region_default", "not_met", None),
        ("stoploss", "moderate", None),
    ],
    "ATM": [
        ("region_default", "not_met", None),
        ("stoploss", "moderate", None),
    ],
    # CYP1B1 (GN104): stop-loss is NOT a disease mechanism → PM4 N/A for stop-loss
    # (in-frame indels keep the default Moderate).
    "CYP1B1": [
        ("stoploss", "not_applicable", None),
    ],
    # VHL (GN078): PM4_Moderate for in-frame indels in the beta domain (residues
    # 63-155 and 193-204) and the alpha domain (155-193) — together 63-204 of the
    # 213-aa protein; outside these domains → N/A. Stop-loss variants adding
    # significant additional amino acids → Moderate (Type 2A extensions).
    "VHL": [
        ("moderate", [(63, 204)], []),
        ("region_default", "not_met", None),
        ("stoploss", "moderate", None),
    ],
    # RPGR (GN106, ORF15 isoform MANE NM_001034853.2, 1152 aa). PM4_Moderate for
    # in-frame indels in exons 1-14 (codons 1-585) or the non-repetitive part of
    # ORF15 (aa 585-1078) — together 1-1078; the repetitive ORF15 C-terminus
    # (1079-1152) → N/A. Stop-loss at the terminal codon (aa 1153) produces a
    # 38-aa extension shown to be deleterious → PVS1... here PM4_Strong.
    "RPGR": [
        ("moderate", [(1, 1078)], []),
        ("region_default", "not_met", None),
        ("stoploss", "strong", None),
    ],
    # RPE65 (GN120): PM4 only when the indel touches a conserved residue
    # (PhyloP>2.0); ≥2 aa → Moderate, 1 aa → Supporting (size via
    # pm4_supporting_max_aa=1). CTLA4 (GN122) / PIK3R1 (GN160): conserved-nucleotide
    # gate (PhyloP>2.0) + mutually exclusive with PVS1/PP3.
    "RPE65": [
        ("conserved_phylop", "2.0", None),
    ],
    "CTLA4": [
        ("conserved_phylop", "2.0", None),
        ("excludes", "PVS1,PP3", None),
    ],
    "PIK3R1": [
        ("conserved_phylop", "2.0", None),
        ("excludes", "PVS1,PP3", None),
    ],
    # KCNQ1 (GN112): PM4 mutually exclusive with PVS1 and PP3. FBN1 (GN022): PM4
    # cannot be applied together with PVS1 at any strength.
    "KCNQ1": [
        ("excludes", "PVS1,PP3", None),
    ],
    "FBN1": [
        ("excludes", "PVS1", None),
    ],
    # SCID panel (GN113/116/119/121/123/124/129): when PM4 is applied to a
    # DELETION, the deleted protein range must contain a known ClinVar P/LP
    # (→Moderate) or VUS (→Supporting) variant not predicted to alter splicing.
    "FOXN1": [("deletion_content", "yes", None)],
    "ADA": [("deletion_content", "yes", None)],
    "DCLRE1C": [("deletion_content", "yes", None)],
    "IL7R": [("deletion_content", "yes", None)],
    "JAK3": [("deletion_content", "yes", None)],
    "RAG1": [("deletion_content", "yes", None)],
    "RAG2": [("deletion_content", "yes", None)],
    "IL2RG": [("deletion_content", "yes", None)],
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="resources/shared/pm4_regions.tsv")
    args = ap.parse_args()

    _SCALAR = ("region_default", "stoploss", "conserved_phylop",
               "deletion_content", "excludes")
    rows = []
    for gene, entries in sorted(_CURATED.items()):
        for strength, regions, residues in entries:
            if strength in _SCALAR:
                rows.append({"gene_symbol": gene, "strength": strength,
                             "regions": regions, "residues": ""})
            else:
                rows.append({
                    "gene_symbol": gene, "strength": strength,
                    "regions": ";".join(f"{a}-{b}" for a, b in regions),
                    "residues": ",".join(str(r) for r in residues),
                })
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["gene_symbol", "strength", "regions", "residues"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    print(f"PM4 region rows: {len(rows)} | genes: {len(_CURATED)} | written → {args.out}")


if __name__ == "__main__":
    main()
