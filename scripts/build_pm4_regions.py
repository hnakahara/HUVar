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
  * ``strong_deletion`` — like ``strong`` but for in-frame deletions only.
  * ``deny_below`` — in these regions PM4 needs an indel of at least N amino acids
    (N in the ``residues`` column).
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
  * ``requires`` — criteria that must also be met (``PM2``), else PM4 is suppressed.
  * ``small_outside`` — in-frame indels of at most this many amino acids are
    Supporting outside the positive regions (inside them the region strength holds).
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
    # Rett / Angelman VCEPs (MECP2 GN036, CDKL5 GN034, FOXG1 GN035, TCF4 GN032,
    # UBE3A GN037). PM4_Moderate: "in-frame deletions/insertions in a non-repeat
    # region"; PM4_Supporting: "smaller in-frame events (< 3 amino acid residues)
    # unless they occur in a functionally important region (see PM1)". So:
    #   * inside a PM1 functional region  -> Moderate at any size (``moderate`` rows);
    #   * elsewhere in a non-repeat region -> Moderate for >= 3 aa, Supporting for
    #     <= 2 aa (``region_default moderate`` + ``small_outside 2``);
    #   * regions the VCEP excludes (``deny``) -> not met.
    # MECP2 / UBE3A: "PM4_Strong is applicable to stop-loss variants ..., as several
    # stop loss variants in this gene has been described in affected individuals"
    # (``stoploss strong``); the other three keep the default Moderate.
    # Coordinates are MANE Select positions, like pm1_hotspots.tsv: the MECP2 spec
    # is written on MeCP2_e2 (NM_004992; MBD 90-162, TRD 302-306, Pro-rich
    # 381-405) and MANE is MeCP2_e1 (+12); the UBE3A spec uses the 852-aa isoform
    # (3' cysteine binding site aa 820) and MANE NM_130839.5 is +20.
    "MECP2": [
        ("moderate", [(102, 174), (314, 318)], []),
        ("deny", [(393, 417)], []),
        ("region_default", "moderate", None),
        ("small_outside", "2", None),
        ("stoploss", "strong", None),
    ],
    # CDKL5: ATP-binding region aa 19-43, TEY site aa 169-171; not for the
    # C-terminus after p.904 (exons 19-21).
    "CDKL5": [
        ("moderate", [(19, 43), (169, 171)], []),
        ("deny", [(905, 960)], []),
        ("region_default", "moderate", None),
        ("small_outside", "2", None),
    ],
    # FOXG1: Forkhead domain aa 181-275; not for the His-rich p.37-57, Pro/Gln-rich
    # p.58-86 and Pro-rich p.105-112 regions.
    "FOXG1": [
        ("moderate", [(181, 275)], []),
        ("deny", [(37, 86), (105, 112)], []),
        ("region_default", "moderate", None),
        ("small_outside", "2", None),
    ],
    # TCF4: bHLH domain aa 564-617.
    "TCF4": [
        ("moderate", [(564, 617)], []),
        ("region_default", "moderate", None),
        ("small_outside", "2", None),
    ],
    "UBE3A": [
        ("moderate", [], [840]),
        ("region_default", "moderate", None),
        ("small_outside", "2", None),
        ("stoploss", "strong", None),
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
    # PTEN (GN003): "Applies to in-frame insertions or deletions impacting at least
    # one residue in a catalytic motif (see PM1), and variants causing protein
    # extension" -> Moderate in the PM1 motifs (WPD loop 90-94, P loop 123-130, TI
    # loop 166-168), N/A elsewhere; stop-loss Moderate.
    "PTEN": [
        ("moderate", [(90, 94), (123, 130), (166, 168)], []),
        ("region_default", "not_met", None),
        ("stoploss", "moderate", None),
    ],
    # NEB (GN146): PM4_Strong "In-frame deletions ..., particularly in exon 55, are
    # deleterious and pathogenic" -> in-frame deletions in exon 55 (MANE
    # NM_001164508.2 c.7432-7536, aa 2478-2512) are Strong; other in-frame events
    # keep the default Moderate.
    "NEB": [
        ("strong_deletion", [(2478, 2512)], []),
    ],
    # RYR1 (GN179, AR congenital myopathy; general text): "Either PVS1_Strong or
    # PM4_Strong, but not both, should be used for in frame deletions in the pore
    # region (Amino acids 4800-4950)" -> in-frame deletions there are PM4_Strong
    # (PM4 is never applied together with PVS1). Not reached under the malignant
    # hyperthermia specification, where PM4 is not applicable.
    "RYR1": [
        ("strong_deletion", [(4800, 4950)], []),
    ],
    # ABCA4 (GN164): besides the size-based in-frame indel rule, PM4 also applies
    # to a synonymous or missense variant at a highly-conserved nucleotide
    # (phyloP >= 7.367): a single changed nucleotide → Supporting, >1 → Moderate.
    "ABCA4": [
        ("nt_phylop", "7.367", None),
    ],
    # RPGR (GN106, ORF15 isoform MANE NM_001034853.2, 1152 aa). PM4_Moderate for
    # in-frame indels in exons 1-14 (codons 1-584) "or, in ORF15 outside of the
    # repetitive region in amino acids 585 to 1078" (aa 1079-1152). The repetitive
    # region itself → N/A. Stop-loss at the terminal codon (aa 1153) produces a
    # 38-aa extension shown to be deleterious → PM4_Strong.
    "RPGR": [
        ("moderate", [(1, 584), (1079, 1152)], []),
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
    # AIPL1 (GN208) / GUCY2D (GN167): same LCA VCEP wording as RPE65 ("loss of at
    # least one conserved residue (PhyloP>2.0) or insertion ... adjacent to" one).
    # AIPL1 general text: "Do not apply to the region from p.Ser328 to p.His384
    # unless the indel is >= 5 amino acids" (``deny_below``: regions, min aa).
    "AIPL1": [
        ("conserved_phylop", "2.0", None),
        ("deny_below", [(328, 384)], [5]),
    ],
    "GUCY2D": [
        ("conserved_phylop", "2.0", None),
    ],
    # LDLR (GN013) / RS1 (GN126): "Variant must also meet PM2".
    "LDLR": [
        ("requires", "PM2", None),
    ],
    "RS1": [
        ("requires", "PM2", None),
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
               "deletion_content", "excludes", "nt_phylop", "requires",
               "small_outside")
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
            delimiter="\t", lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)
    print(f"PM4 region rows: {len(rows)} | genes: {len(_CURATED)} | written → {args.out}")


if __name__ == "__main__":
    main()
