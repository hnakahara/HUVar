"""Gene-specific PVS1 decision trees from ClinGen VCEP specifications.

Several VCEPs publish a gene-specific PVS1 decision tree that deviates from the
generic ClinGen SVI 2019 tree (Abou Tayoun et al. 2018). The deviations that
cause the generic tree to *under-call* (PVS1 false-negatives) are:

* **Critical-region / codon-range truncation gates** — for genes where NMD is
  not predicted (single-exon genes such as RAG1/GP9) or where the truncation
  escapes NMD in the last exon, the VCEP still grants PVS1/Strong/Moderate based
  on whether the truncation removes a critical domain (rather than withholding
  PVS1 as the generic tree does when no protein domain is annotated).
* **Initiation-codon strength overrides** — the generic tree fixes start-loss at
  Moderate, but VCEPs assign Strong (RPE65), Supporting (GCK) or N/A (VHL, where
  the downstream p19 start codon at Met54 rescues a Met1 loss).
* **Canonical ±1,2 splice / whole-gene deletion** — the VCEPs grant PVS1
  (Very Strong) for these prediction-independently; the generic tree can downgrade
  them via the ClinVar-count strength caps or when no splice predictor is loaded.

Because the handler returns a *final* strength (like :mod:`acmg_classifier.pvs1.apc`),
it runs BEFORE the generic tree and its strength caps, which is what removes the
false-negatives.

Each gene's transcript and protein length (MANE Select) is recorded in the spec
for traceability; codon bands are 1-based and inclusive on the spec transcript.

Out of scope (deferred to the generic tree or returned as ``None``): exon-level
duplication "proven/presumed in tandem" distinctions and RNA-evidence strength
modifiers, which require data not available at this layer.
"""
from __future__ import annotations

from dataclasses import dataclass

from acmg_classifier.models.enums import ConsequenceType, CriterionStrength

_S = CriterionStrength

# Truncating consequences gated by codon range.
_TRUNC = (ConsequenceType.STOP_GAINED, ConsequenceType.FRAMESHIFT)
_SPLICE = (ConsequenceType.SPLICE_DONOR, ConsequenceType.SPLICE_ACCEPTOR)

_LABEL = {
    _S.VERY_STRONG: "PVS1",
    _S.STRONG: "PVS1_Strong",
    _S.MODERATE: "PVS1_Moderate",
    _S.SUPPORTING: "PVS1_Supporting",
    _S.NOT_MET: "N/A",
}


_Band = tuple[int, int, CriterionStrength]


@dataclass(frozen=True)
class _GeneSpec:
    """One VCEP's PVS1 deviations, in a form the dispatcher can evaluate.

    Truncating (nonsense / frameshift) variants are handled by exactly one of:

    * ``trunc_bands`` — ordered ``(lo_codon, hi_codon, strength)`` tuples
      (inclusive). A variant whose ``protein_position`` falls in a band gets that
      strength; a codon matching no band → ``NOT_MET``. ``trunc_bands_fs``, when
      set, is used for frameshift variants instead (HNF1A scores nonsense and
      frameshift differently in the TAD/exon-10 region).
    * ``trunc_nmd`` — ``(nmd_predicted_strength, nmd_escape_strength)``. Used when
      the cutoff is the exon-based NMD boundary rather than a fixed codon (FBN1:
      NMD → Very Strong; last-exon / penultimate-55nt escape → Strong because the
      C-terminus is proven critical). When the escape element is ``None`` the
      escape strength follows the generic 10%-of-protein rule (Strong if the
      truncation removes >10% of the protein, else Moderate) — ACADVL, GAMT.

    ``start_lost`` / ``splice`` / ``deletion``: strength for that consequence, or
    ``None`` to defer to the generic decision tree.

    ``splice_exclude_donor_introns``: donor (``+1/+2``) splice sites of these
    introns do NOT get PVS1 (ACADVL intron 8 begins with GC, not GT, so its
    impact is not well understood and the VCEP withholds PVS1).
    """

    gene: str
    transcript: str
    aa_len: int
    start_lost: CriterionStrength | None
    splice: CriterionStrength | None
    deletion: CriterionStrength | None
    trunc_bands: tuple[_Band, ...] | None = None
    trunc_bands_fs: tuple[_Band, ...] | None = None
    trunc_nmd: tuple[CriterionStrength, CriterionStrength | None] | None = None
    splice_exclude_donor_introns: frozenset[int] = frozenset()


_SPECS: dict[str, _GeneSpec] = {
    # RPE65 (Leber congenital amaurosis VCEP, Walker 2023-adapted tree).
    # All 14 exons critical; Met1 loss → Strong (2nd Met at 93, no alt start
    # evidence); nonsense/frameshift p.Ser2–p.Gly528 → PVS1, p.Leu529–p.Ser533
    # (last residues) → PVS1_Strong.
    "RPE65": _GeneSpec(
        gene="RPE65", transcript="NM_000329.3", aa_len=533,
        trunc_bands=((2, 528, _S.VERY_STRONG), (529, 533, _S.STRONG)),
        start_lost=_S.STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # CYP1B1 (Glaucoma VCEP, Abou Tayoun-adapted tree). NMD before aa330; the
    # haem-binding domain (aa460-493) is vital, so truncations through aa493 are
    # null (PVS1); aa494-Ter do not remove the domain → PVS1_Moderate. Initiation
    # codon → PVS1_Moderate (no known alt start, pathogenic variants upstream).
    "CYP1B1": _GeneSpec(
        gene="CYP1B1", transcript="NM_000104.4", aa_len=543,
        trunc_bands=((1, 493, _S.VERY_STRONG), (494, 543, _S.MODERATE)),
        start_lost=_S.MODERATE, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # VHL VCEP. No PVS1 for truncations prior to codon 54 (downstream p19 start
    # at Met54 rescues). After Met54: NMD region / beta+alpha critical domains
    # (aa54-204) → PVS1; outside the 2nd beta domain (aa205-213) → PVS1_Moderate.
    # Any canonical exon skip or exon deletion → PVS1. Met1 start-loss → N/A
    # (p19 still produced). The 10% rule does not apply (small protein).
    "VHL": _GeneSpec(
        gene="VHL", transcript="NM_000551.4", aa_len=213,
        trunc_bands=((54, 204, _S.VERY_STRONG), (205, 213, _S.MODERATE)),
        start_lost=_S.NOT_MET, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # GCK (MDEP). PTCs throughout, including the last exon (exon 10) and the
    # last 55 nt of exon 9 that escape NMD, cause MODY → PVS1 (Very Strong).
    # Initiation codon → PVS1_Supporting (single VUS reviewed; next Met at 8).
    "GCK": _GeneSpec(
        gene="GCK", transcript="NM_000162.5", aa_len=465,
        trunc_bands=((1, 465, _S.VERY_STRONG),),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # RAG1 (SCID VCEP). Single-exon gene → NMD never predicted, so PVS1 only for
    # full-gene deletion or removing/altering a critical domain. The critical
    # region (NBD 394-460, DDBD 461-517, core 387-1011) extends to aa1011; any
    # truncation at/before 1011 removes part of it → PVS1. After 1011 only the
    # last <10% remains → N/A. No canonical splice (single exon).
    "RAG1": _GeneSpec(
        gene="RAG1", transcript="NM_000448.3", aa_len=1043,
        trunc_bands=((1, 1011, _S.VERY_STRONG), (1012, 1043, _S.NOT_MET)),
        start_lost=None, splice=None, deletion=_S.VERY_STRONG,
    ),
    # ATM (Hereditary breast cancer / ENIGMA). The C-terminal FATKIN
    # (FAT+kinase+FATC) domain is critical; a truncation anywhere up to the most
    # 3' pathogenic residue p.R3047 removes the kinase domain → PVS1. After R3047
    # (last <0.3%) no pathogenic truncations are known → N/A.
    "ATM": _GeneSpec(
        gene="ATM", transcript="NM_000051.4", aa_len=3056,
        trunc_bands=((1, 3047, _S.VERY_STRONG), (3048, 3056, _S.NOT_MET)),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # GP9 (Platelet disorder VCEP). Single coding exon → not NMD-subject. The
    # transmembrane domain (aa148-169) is critical; a truncation at/before 169
    # removes it or >10% of protein → PVS1_Strong; at/after p.170 only <10%
    # remains with unknown role → PVS1_Moderate. Initiation codon → PVS1_Moderate
    # (pathogenic c.70T>C upstream of the Met32 alt start).
    "GP9": _GeneSpec(
        gene="GP9", transcript="NM_000174.5", aa_len=177,
        trunc_bands=((1, 169, _S.STRONG), (170, 177, _S.MODERATE)),
        start_lost=_S.MODERATE, splice=None, deletion=_S.VERY_STRONG,
    ),
    # IDUA (Lysosomal diseases VCEP). NMD predicted for PTCs 5' of c.1778
    # (≈ codon 593); at/after c.1778 (last 50 nt of exon 13 + exon 14) NMD is
    # escaped and only the last <10% with unknown role remains → PVS1_Moderate.
    "IDUA": _GeneSpec(
        gene="IDUA", transcript="NM_000203.5", aa_len=653,
        trunc_bands=((1, 592, _S.VERY_STRONG), (593, 653, _S.MODERATE)),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # ACVRL1 (Hereditary haemorrhagic telangiectasia VCEP). NMD predicted
    # ≤ codon 442 → PVS1; NMD escape but critical region ≤ codon 490 → PVS1_Strong;
    # > codon 490 (unknown role, <10%) → PVS1_Moderate. Initiation codon →
    # PVS1_Moderate (no known alt start, unknown role).
    "ACVRL1": _GeneSpec(
        gene="ACVRL1", transcript="NM_000020.3", aa_len=503,
        trunc_bands=((1, 442, _S.VERY_STRONG), (443, 490, _S.STRONG),
                     (491, 503, _S.MODERATE)),
        start_lost=_S.MODERATE, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # PAH (PAH VCEP, Tayoun / Walker). Nonsense/frameshift upstream of c.1285
    # (codon ≤428) → PVS1; downstream (codon ≥429) → PVS1_Strong. Initiation
    # codon → PVS1_Strong.
    "PAH": _GeneSpec(
        gene="PAH", transcript="NM_000277.3", aa_len=452,
        trunc_bands=((1, 428, _S.VERY_STRONG), (429, 452, _S.STRONG)),
        start_lost=_S.STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # HNF1A (MDEP). Nonsense/frameshift 5' of c.1768 (codon ≤589) → PVS1
    # regardless of where the PTC lands (TAD overlaps the last 55 nt of exon 9 +
    # exon 10). 3' of c.1768 the nonsense vs frameshift cutoffs differ: nonsense
    # ≤p.601 → Strong, >p.601 → Supporting; frameshift ≤p.618 → Strong, >p.618 →
    # Supporting (added residues from a frameshift disrupt the TAD more than a
    # short C-terminal nonsense truncation). Initiation codon → PVS1.
    "HNF1A": _GeneSpec(
        gene="HNF1A", transcript="NM_000545.8", aa_len=631,
        trunc_bands=((1, 589, _S.VERY_STRONG), (590, 601, _S.STRONG),
                     (602, 631, _S.SUPPORTING)),
        trunc_bands_fs=((1, 589, _S.VERY_STRONG), (590, 618, _S.STRONG),
                        (619, 631, _S.SUPPORTING)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # GJB2 (Hearing Loss VCEP). Connexin-26; established LoF. The ORF is in a
    # single coding exon (exon 2 = last exon) so NMD is never predicted — the
    # generic tree under-calls. Truncations removing >10% of the protein
    # (codon ≤204) → PVS1; the last <10% (role of C-terminus less established) →
    # PVS1_Moderate.
    "GJB2": _GeneSpec(
        gene="GJB2", transcript="NM_004004.6", aa_len=226,
        trunc_bands=((1, 204, _S.VERY_STRONG), (205, 226, _S.MODERATE)),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # FOXG1 (Rett/Angelman VCEP). Intronless (single-exon) gene. Null variants up
    # to p.S468 → PVS1; truncating p.G469–p.Q480 → PVS1_Strong; distal of p.Q480
    # → PVS1_Moderate. Initiation codon → PVS1_Supporting. Full gene deletion →
    # PVS1. No canonical splice (single exon).
    "FOXG1": _GeneSpec(
        gene="FOXG1", transcript="NM_005249.5", aa_len=489,
        trunc_bands=((1, 468, _S.VERY_STRONG), (469, 480, _S.STRONG),
                     (481, 489, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=None, deletion=_S.VERY_STRONG,
    ),
    # DICER1 VCEP. NMD cutoff at p.Pro1850: PTC ≤1850 → PVS1; truncation 3' of it
    # → PVS1_Moderate. Initiation-codon (p.M1?) → no criterion (M1 not conserved,
    # three in-frame alt starts at Met11/17/24, unaffected cases). Canonical
    # splice → PVS1 (in-frame exon-skip exceptions to Strong/Moderate are not
    # modelled at this layer).
    "DICER1": _GeneSpec(
        gene="DICER1", transcript="NM_177438.3", aa_len=1922,
        trunc_bands=((1, 1850, _S.VERY_STRONG), (1851, 1922, _S.MODERATE)),
        start_lost=_S.NOT_MET, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # PALB2 (HBOP VCEP). The WD40 beta-propeller (C-terminal, through the very
    # last residues — Y1183/H1184/Y1185/S1186 form the "molecular Velcro" that
    # seals the toroid) and the coiled-coil are indispensable, so even NMD-escape
    # truncations that adversely affect WD40 are PVS1 (not the Strong baseline).
    # Because WD40 reaches the C-terminus, any truncation → PVS1. Canonical
    # splice → PVS1.
    "PALB2": _GeneSpec(
        gene="PALB2", transcript="NM_024675.4", aa_len=1186,
        trunc_bands=((1, 1186, _S.VERY_STRONG),),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # FBN1 (Marfan VCEP). NMD-based: nonsense/frameshift predicted to undergo NMD
    # → PVS1; predicted to escape NMD (last exon / last 55 nt of penultimate exon)
    # → PVS1_Strong (the C-terminal region is proven critical, so escape variants
    # are not withheld). Initiation codon → PVS1_Moderate. Frame-dependent splice
    # is left to the generic tree (many in-frame cbEGF exons).
    "FBN1": _GeneSpec(
        gene="FBN1", transcript="NM_000138.5", aa_len=2871,
        trunc_nmd=(_S.VERY_STRONG, _S.STRONG),
        start_lost=_S.MODERATE, splice=None, deletion=_S.VERY_STRONG,
    ),
    # GP1BA (Platelet disorder VCEP). Single coding exon → not NMD-subject (cf.
    # GP9). Transmembrane domain aa532-553 critical; a truncation at/before p.587
    # removes the TM domain or >10% of protein → PVS1_Strong; at/after p.588 only
    # <10% with unknown role remains → PVS1_Moderate. Initiation codon →
    # PVS1_Moderate (closest alt start Met68).
    "GP1BA": _GeneSpec(
        gene="GP1BA", transcript="NM_000173.7", aa_len=652,
        trunc_bands=((1, 587, _S.STRONG), (588, 652, _S.MODERATE)),
        start_lost=_S.MODERATE, splice=None, deletion=_S.VERY_STRONG,
    ),
    # CDH1 (HDGC VCEP, v3). Nonsense/frameshift predicted to undergo NMD → PVS1;
    # NMD escape (last exon 16 / last 50 nt of exon 15) → PVS1_Strong (the
    # cytoplasmic C-terminal region is critical). Canonical ±1,2 splice DEFAULT is
    # PVS1_Strong (per the v3 caveat / splicing table), not Very Strong.
    "CDH1": _GeneSpec(
        gene="CDH1", transcript="NM_004360.5", aa_len=882,
        trunc_nmd=(_S.VERY_STRONG, _S.STRONG),
        start_lost=None, splice=_S.STRONG, deletion=_S.VERY_STRONG,
    ),
    # AIPL1 (LCA VCEP, Walker 2023-adapted). All 6 exons critical. Nonsense and
    # frameshift score differently in the C-terminus (frameshift adds residues
    # via long downstream ORFs): nonsense p.Thr2–p.Ser328 → PVS1, p.Glu329–p.Ser346
    # → Strong, p.Ser347–p.His384 → Moderate; frameshift p.Thr2–p.Glu337 → PVS1,
    # p.Pro338–p.His384 → Strong. Initiation codon → PVS1 (2nd Met at 40, region
    # functionally important). Canonical ±1,2 splice and exon deletions → PVS1.
    "AIPL1": _GeneSpec(
        gene="AIPL1", transcript="NM_014336.5", aa_len=384,
        trunc_bands=((2, 328, _S.VERY_STRONG), (329, 346, _S.STRONG),
                     (347, 384, _S.MODERATE)),
        trunc_bands_fs=((2, 337, _S.VERY_STRONG), (338, 384, _S.STRONG)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # ACADVL (VLCAD VCEP). NMD-based: nonsense/frameshift NMD → PVS1; escape
    # (last exon 20 / last 50 nt of exon 19) → Strong if >10% of protein removed
    # else Moderate. Initiation codon → PVS1_Strong (leader sequence aa1-40 is
    # important; next Met at 6). Canonical ±1,2 splice → PVS1, EXCEPT the donor
    # of intron 8 (a GC donor, impact unclear → no PVS1).
    "ACADVL": _GeneSpec(
        gene="ACADVL", transcript="NM_000018.4", aa_len=655,
        trunc_nmd=(_S.VERY_STRONG, None),
        start_lost=_S.STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
        splice_exclude_donor_introns=frozenset({8}),
    ),
    # TP53 VCEP. NMD cutoff at p.Lys351: nonsense/frameshift PTC ≤350 → PVS1;
    # p.351–p.355 → PVS1_Strong; p.356–p.393 → PVS1_Moderate. Initiation codon →
    # PVS1. Canonical ±1,2 splice → PVS1 (E10-donor / E11-acceptor C-terminal-
    # shortening Moderate exceptions not modelled). Full gene deletion → PVS1.
    "TP53": _GeneSpec(
        gene="TP53", transcript="NM_000546.6", aa_len=393,
        trunc_bands=((1, 350, _S.VERY_STRONG), (351, 355, _S.STRONG),
                     (356, 393, _S.MODERATE)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # GAA (Pompe VCEP). Any nonsense/frameshift/splice PTC before codon 916 →
    # PVS1; PTC 3' of codon 916 (NMD escape) → PVS1_Moderate. Initiation codon →
    # PVS1_Strong. Canonical ±1,2 splice → PVS1 (in-frame exon-skip exceptions
    # to Strong/Moderate not modelled). Full gene deletion → PVS1.
    "GAA": _GeneSpec(
        gene="GAA", transcript="NM_000152.5", aa_len=952,
        trunc_bands=((1, 915, _S.VERY_STRONG), (916, 952, _S.MODERATE)),
        start_lost=_S.STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # GAMT (Cerebral creatine deficiency VCEP). NMD-based: nonsense/frameshift NMD
    # → PVS1; escape (last exon 6 / last 50 nt of exon 5, c.520) → Strong if >10%
    # removed else Moderate. Initiation codon → PVS1_Moderate (next Met at 42).
    # Canonical ±1,2 splice → PVS1; full gene deletion → PVS1.
    "GAMT": _GeneSpec(
        gene="GAMT", transcript="NM_000156.6", aa_len=236,
        trunc_nmd=(_S.VERY_STRONG, None),
        start_lost=_S.MODERATE, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # HNF4A (MDEP). PTCs in exon 10 + last 55 nt of exon 9 (c.1162-1216) escape
    # NMD, but the collective evidence supports PVS1 for any nonsense/frameshift
    # at codon 419 (c.1257) and 5'; at p.Gly420 (c.1258) and 3' → PVS1_Supporting.
    # Initiation codon → PVS1_Strong (P/LP variants upstream of the next Met71).
    # Canonical ±1,2 splice → PVS1 (in-frame exon-skip Strong exceptions for
    # exons 5/7/8/9/10 not modelled). Full gene deletion → PVS1.
    # MANE Select is NM_175914.5 (452 aa); the MDEP c./p. numbering (c.1217 starts
    # exon 10, last 55 nt of exon 9 = c.1162-1216, p.S419 = c.1255-1257) matches
    # this transcript exactly — verified against the generated exon table.
    "HNF4A": _GeneSpec(
        gene="HNF4A", transcript="NM_175914.5", aa_len=452,
        trunc_bands=((1, 419, _S.VERY_STRONG), (420, 452, _S.SUPPORTING)),
        start_lost=_S.STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # RUNX1 (Myeloid Malignancy VCEP, isoform c NM_001754.4). NMD-based:
    # nonsense/frameshift predicted to undergo NMD → PVS1; C-terminal truncating
    # variants NOT predicted to undergo NMD → PVS1_Strong. Canonical ±1,2 splice
    # → PVS1. Full gene deletion → PVS1 (exon 1-3 isoform-c-only deletions are
    # PVS1_Moderate but require exon-level CNV detail not modelled here).
    "RUNX1": _GeneSpec(
        gene="RUNX1", transcript="NM_001754.5", aa_len=480,
        trunc_nmd=(_S.VERY_STRONG, _S.STRONG),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # CDKL5 (Rett/Angelman VCEP, major brain isoform NM_001323289.2). Null
    # variants up to p.R948 → PVS1; any truncating variant distal of p.R948 →
    # PVS1_Moderate. Initiation codon → PVS1_Supporting. Canonical ±1,2 splice →
    # PVS1 (exon-18-flanking Strong / exon-17-flanking Moderate exceptions not
    # modelled). Full gene deletion → PVS1.
    "CDKL5": _GeneSpec(
        gene="CDKL5", transcript="NM_001323289.2", aa_len=960,
        trunc_bands=((1, 948, _S.VERY_STRONG), (949, 960, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # RPGR (eyeGENE / LCA-RP VCEP). MANE Select is the retinal ORF15 isoform
    # NM_001034853.2 (1152 aa, 15 coding exons), which is the VCEP's default. Its
    # terminal exon 15 (ORF15) is ~49% of the protein and its glutamylation
    # function is critical, so an NMD-escape truncation there disrupts
    # glutamylation → PVS1; NMD-predicted truncations (exons 1-14) → PVS1. Hence
    # every nonsense/frameshift → PVS1. Initiation codon → PVS1_Moderate (RCC1
    # domain pathogenic variants upstream). Canonical ±1,2 splice / full gene
    # deletion → PVS1.
    "RPGR": _GeneSpec(
        gene="RPGR", transcript="NM_001034853.2", aa_len=1152,
        trunc_bands=((1, 1152, _S.VERY_STRONG),),
        start_lost=_S.MODERATE, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # IL2RG (SCID VCEP, X-linked). The transmembrane domain begins at aa255 and
    # the distal cytoplasmic domain is critical, so an NMD-escape truncation that
    # reaches them is PVS1 (not the Strong baseline). Because the NMD-escape zone
    # (last exon 8 / last 50 nt of exon 7, codon ≥292) lies entirely within /
    # distal to the TM domain, every nonsense/frameshift truncation → PVS1.
    "IL2RG": _GeneSpec(
        gene="IL2RG", transcript="NM_000206.3", aa_len=369,
        trunc_bands=((1, 369, _S.VERY_STRONG),),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # MECP2 (Rett/Angelman VCEP). MANE Select is the brain-predominant MeCP2_e1
    # isoform NM_001110792.2 (498 aa); the VCEP rule is written in MeCP2_e2
    # (NM_004992) numbering (null up to p.E472, distal → Moderate). The shared
    # C-terminus is identical with a +12-residue offset (e1 = e2 + 12), so e2
    # p.E472 maps to e1 p.484. Initiation codon → N/A (MeCP2_e1 / e2 use
    # different start codons). Canonical ±1,2 splice → PVS1; full gene deletion
    # → PVS1.
    "MECP2": _GeneSpec(
        gene="MECP2", transcript="NM_001110792.2", aa_len=498,
        trunc_bands=((1, 484, _S.VERY_STRONG), (485, 498, _S.MODERATE)),
        start_lost=_S.NOT_MET, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # F9 (Coagulation factor / Hemophilia B VCEP, X-linked). NMD-escape
    # truncations remain critical (many severe hemophilia B cases in the EAHAD
    # database across the escape region), so every nonsense/frameshift → PVS1.
    # Initiation codon → PVS1_Supporting (in-frame ATGs at codons 6/8, no
    # pathogenic variants upstream). Canonical ±1,2 splice → PVS1; full gene
    # deletion → PVS1.
    "F9": _GeneSpec(
        gene="F9", transcript="NM_000133.4", aa_len=461,
        trunc_bands=((1, 461, _S.VERY_STRONG),),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # ABCD1 (X-linked adrenoleukodystrophy VCEP). NMD-based with the generic 10%
    # escape rule: nonsense/frameshift NMD (nonsense ≤c.1941 ≈ codon 647) → PVS1;
    # NMD escape → Strong if >10% removed else Moderate. Initiation codon →
    # PVS1_Moderate (≥1 pathogenic variant upstream of the closest in-frame
    # start). Canonical ±1,2 splice → PVS1; full gene deletion → PVS1.
    "ABCD1": _GeneSpec(
        gene="ABCD1", transcript="NM_000033.4", aa_len=745,
        trunc_nmd=(_S.VERY_STRONG, None),
        start_lost=_S.MODERATE, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),

    # ===================================================================
    # Batch from ClinGen cspec re-examination (genes with explicit codon /
    # domain / NMD-escape PVS1 rules not previously covered).
    # ===================================================================

    # ADA, DCLRE1C, JAK3 (SCID VCEP): standard SVI tree with the NMD-escape >10%
    # rule (Strong if >10% removed — VCEP also requires a downstream pathogenic
    # variant, approximated here — else Moderate).
    "ADA": _GeneSpec(
        gene="ADA", transcript="NM_000022.4", aa_len=363,
        trunc_nmd=(_S.VERY_STRONG, None),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    "DCLRE1C": _GeneSpec(
        gene="DCLRE1C", transcript="NM_001033855.3", aa_len=692,
        trunc_nmd=(_S.VERY_STRONG, None),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    "JAK3": _GeneSpec(
        gene="JAK3", transcript="NM_000215.4", aa_len=1124,
        trunc_nmd=(_S.VERY_STRONG, None),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # IL7R (SCID VCEP): like IL2RG, the transmembrane domain (begins aa240) and
    # distal cytoplasmic region are critical, and the NMD-escape zone (last exon 8
    # / last 50 nt of exon 7, codon ≥276) lies entirely within/distal to the TM
    # domain, so every NMD-escape truncation → Strong.
    "IL7R": _GeneSpec(
        gene="IL7R", transcript="NM_002185.5", aa_len=459,
        trunc_nmd=(_S.VERY_STRONG, _S.STRONG),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # FOXN1 (SCID VCEP): NMD (≤ codon 526, c.1577) → PVS1; NMD escape removing the
    # transactivation domain (aa511-563) or >10% of protein → Strong; distal of
    # that (≤10%) → Moderate. (The forkhead domain aa270-367 lies in the
    # NMD-predicted region, so it is already PVS1.)
    "FOXN1": _GeneSpec(
        gene="FOXN1", transcript="NM_001369369.1", aa_len=648,
        trunc_bands=((1, 526, _S.VERY_STRONG), (527, 584, _S.STRONG),
                     (585, 648, _S.MODERATE)),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # RAG2 (SCID VCEP): single-exon gene → NMD never predicted. PVS1 only for full
    # gene deletion or removing/altering the critical core (aa1-383) + PHD
    # (aa414-487) domains; a truncation at/before aa487 removes part of them →
    # PVS1; distal (the last <10%) → N/A. No canonical splice (single exon).
    "RAG2": _GeneSpec(
        gene="RAG2", transcript="NM_000536.4", aa_len=527,
        trunc_bands=((1, 487, _S.VERY_STRONG), (488, 527, _S.NOT_MET)),
        start_lost=None, splice=None, deletion=_S.VERY_STRONG,
    ),
    # CTLA4 (CHAI VCEP): truncations in exon 3 — NMD ≤ codon 172 → PVS1; escape
    # codon 173-201 → Strong; codon 202-223 (<10% C-terminus) → Moderate.
    # Initiation codon → Supporting (closest in-frame start at codon 38).
    "CTLA4": _GeneSpec(
        gene="CTLA4", transcript="NM_005214.5", aa_len=223,
        trunc_bands=((1, 172, _S.VERY_STRONG), (173, 201, _S.STRONG),
                     (202, 223, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # KCNQ1 (LQTS VCEP): truncations codons 1-581 → PVS1 (NMD); 582-620 →
    # PVS1_Moderate (NMD escape but removes the subunit-assembly domain SAD
    # 589-620); 621-676 → PVS1_Supporting (SAD retained).
    "KCNQ1": _GeneSpec(
        gene="KCNQ1", transcript="NM_000218.3", aa_len=676,
        trunc_bands=((1, 581, _S.VERY_STRONG), (582, 620, _S.MODERATE),
                     (621, 676, _S.SUPPORTING)),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # Lynch-syndrome MMR genes (InSiGHT VCEP): nonsense/frameshift PTC at/before
    # the cutoff codon → PVS1; a short window beyond it → PVS1_Moderate. Large
    # single/multi-exon deletions → PVS1.
    "MLH1": _GeneSpec(
        gene="MLH1", transcript="NM_000249.4", aa_len=756,
        trunc_bands=((1, 753, _S.VERY_STRONG), (754, 756, _S.MODERATE)),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    "MSH2": _GeneSpec(
        gene="MSH2", transcript="NM_000251.3", aa_len=934,
        trunc_bands=((1, 891, _S.VERY_STRONG), (892, 934, _S.MODERATE)),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    "MSH6": _GeneSpec(
        gene="MSH6", transcript="NM_000179.3", aa_len=1360,
        trunc_bands=((1, 1341, _S.VERY_STRONG), (1342, 1360, _S.MODERATE)),
        start_lost=_S.STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    "PMS2": _GeneSpec(
        gene="PMS2", transcript="NM_000535.7", aa_len=862,
        trunc_bands=((1, 798, _S.VERY_STRONG), (799, 862, _S.MODERATE)),
        start_lost=_S.STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # OTC (Urea cycle VCEP, X-linked): NMD truncations (before c.1033 ≈ codon 345)
    # → PVS1; nonsense/frameshift/splice at c.1033 (codon 345) and downstream →
    # PVS1_Strong (the C-terminal region is critical, with reported pathogenic
    # frameshift/stop-loss/extension variants).
    "OTC": _GeneSpec(
        gene="OTC", transcript="NM_000531.6", aa_len=354,
        trunc_bands=((1, 344, _S.VERY_STRONG), (345, 354, _S.STRONG)),
        start_lost=None, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # SLC9A6 (Christianson syndrome VCEP, X-linked). The VCEP residue numbering
    # (null ≤p.A563, p.C564-T601 → Strong, p.Y602-A669 → Moderate) is in the
    # 669-aa transcript; MANE Select NM_001379110.1 is 679 aa (N-terminal isoform
    # difference), so the cutoffs are shifted +10. Initiation codon → Supporting.
    "SLC9A6": _GeneSpec(
        gene="SLC9A6", transcript="NM_001379110.1", aa_len=679,
        trunc_bands=((1, 573, _S.VERY_STRONG), (574, 611, _S.STRONG),
                     (612, 679, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # TCF4 (Pitt-Hopkins VCEP): null variants up to p.E643 → PVS1; any truncating
    # variant distal of p.E643 → PVS1_Moderate. Initiation codon → Supporting.
    "TCF4": _GeneSpec(
        gene="TCF4", transcript="NM_001083962.2", aa_len=671,
        trunc_bands=((1, 643, _S.VERY_STRONG), (644, 671, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # UBE3A (Angelman VCEP, maternally expressed). The VCEP numbering (PVS1 up to
    # p.K841, p.A842-G850 → Strong, distal of p.G850 → Moderate) is in the 852-aa
    # isoform; MANE Select NM_130839.5 is 872 aa (N-terminal isoform difference),
    # so the cutoffs are shifted +20. Initiation codon → PVS1.
    "UBE3A": _GeneSpec(
        gene="UBE3A", transcript="NM_130839.5", aa_len=872,
        trunc_bands=((1, 861, _S.VERY_STRONG), (862, 870, _S.STRONG),
                     (871, 872, _S.MODERATE)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # GUCY2D (LCA VCEP, Walker-adapted): all exons critical. Nonsense/frameshift
    # up to p.Pro1068 → PVS1; p.Pro1069-p.Ser1103 → PVS1_Strong. Canonical ±1,2
    # splice (exons 2-20) and exon deletions → PVS1. Initiation codon → PVS1
    # (2nd Met at residue 218).
    "GUCY2D": _GeneSpec(
        gene="GUCY2D", transcript="NM_000180.4", aa_len=1103,
        trunc_bands=((1, 1068, _S.VERY_STRONG), (1069, 1103, _S.STRONG)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # RS1 (X-linked retinoschisis VCEP): nonsense/frameshift/splice/deletion in
    # c.1A-c.671 (p.Met1-Cys223) → PVS1; the final residues c.672-677
    # (p.Asp224-*225) → PVS1_Strong. The discoidin structure / disulfide bonds
    # make the whole protein critical. Initiation codon (c.1A) → PVS1.
    "RS1": _GeneSpec(
        gene="RS1", transcript="NM_000330.4", aa_len=224,
        trunc_bands=((1, 223, _S.VERY_STRONG), (224, 224, _S.STRONG)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),

    # ===================================================================
    # Batch from supplied VCEP PVS1 decision-tree files.
    # ===================================================================

    # ENG (HHT VCEP): NMD ≤ codon 601 → PVS1; NMD escape (role unknown, <10%) →
    # PVS1_Moderate. Initiation codon → PVS1_Strong.
    "ENG": _GeneSpec(
        gene="ENG", transcript="NM_001114753.3", aa_len=658,
        trunc_bands=((1, 601, _S.VERY_STRONG), (602, 658, _S.MODERATE)),
        start_lost=_S.STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # GP1BB (Platelet VCEP): coding starts in the last 50 nt of the penultimate
    # exon → not NMD-subject. TM domain aa148-173 critical; truncation before
    # p.186 (>10%) or in the TM domain → PVS1_Strong; at/after p.186 (<10%) →
    # PVS1_Moderate. Initiation codon → PVS1_Supporting (no alt start codon).
    "GP1BB": _GeneSpec(
        gene="GP1BB", transcript="NM_000407.5", aa_len=206,
        trunc_bands=((1, 185, _S.STRONG), (186, 206, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=None, deletion=_S.VERY_STRONG,
    ),
    # SCN1B (Epilepsy VCEP): NMD (stop 5' of p.Thr204) → PVS1; escape (3' of
    # p.Thr204, <10% C-terminus) → PVS1_Moderate. Initiation codon → Supporting.
    "SCN1B": _GeneSpec(
        gene="SCN1B", transcript="NM_001037.5", aa_len=218,
        trunc_bands=((1, 204, _S.VERY_STRONG), (205, 218, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # SCN2A / SCN3A / SCN8A (Epilepsy VCEP): NMD (stop 5' of the boundary codon)
    # → PVS1; NMD escape → Strong if >10% (≥~200 aa) removed else Moderate.
    "SCN2A": _GeneSpec(
        gene="SCN2A", transcript="NM_001371246.1", aa_len=2005,
        trunc_bands=((1, 1591, _S.VERY_STRONG), (1592, 1805, _S.STRONG),
                     (1806, 2005, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    "SCN3A": _GeneSpec(
        gene="SCN3A", transcript="NM_006922.4", aa_len=2000,
        trunc_bands=((1, 1586, _S.VERY_STRONG), (1587, 1801, _S.STRONG),
                     (1802, 2000, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    "SCN8A": _GeneSpec(
        gene="SCN8A", transcript="NM_014191.4", aa_len=1980,
        trunc_bands=((1, 1582, _S.VERY_STRONG), (1583, 1783, _S.STRONG),
                     (1784, 1980, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # NEB (nemaline myopathy): NMD (PTC ≤ p.8452, c.25354) → PVS1; NMD escape
    # (the last <10%) → PVS1_Moderate. Initiation codon → Supporting.
    "NEB": _GeneSpec(
        gene="NEB", transcript="NM_001164508.2", aa_len=8525,
        trunc_bands=((1, 8452, _S.VERY_STRONG), (8453, 8525, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # F8 (Coagulation Factor / Hemophilia A VCEP, X-linked): NMD (nonsense ≤
    # c.6851 ≈ codon 2283) → PVS1; NMD escape (C-terminal C2 domain critical) →
    # PVS1_Strong. Initiation codon → PVS1_Moderate (≥1 pathogenic upstream of
    # the c.96/codon-32 alt start).
    "F8": _GeneSpec(
        gene="F8", transcript="NM_000132.4", aa_len=2351,
        trunc_bands=((1, 2283, _S.VERY_STRONG), (2284, 2351, _S.STRONG)),
        start_lost=_S.MODERATE, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # PTEN VCEP: NMD (stop at/5' to p.D375, c.1121) → PVS1; NMD escape (role
    # unknown, <10%) → PVS1_Moderate. Initiation codon → PVS1.
    "PTEN": _GeneSpec(
        gene="PTEN", transcript="NM_000314.8", aa_len=403,
        trunc_bands=((1, 375, _S.VERY_STRONG), (376, 403, _S.MODERATE)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # MYBPC3 (Cardiomyopathy VCEP): LoF is established for MYBPC3. NMD (prior to
    # p.1254) → PVS1; NMD escape (role unknown, <10% C-terminus) →
    # PVS1_Moderate. Initiation codon → PVS1_Supporting.
    "MYBPC3": _GeneSpec(
        gene="MYBPC3", transcript="NM_000256.3", aa_len=1274,
        trunc_bands=((1, 1253, _S.VERY_STRONG), (1254, 1274, _S.MODERATE)),
        start_lost=_S.SUPPORTING, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # HBB (Hemoglobinopathy VCEP, β-globin). NMD only for PTCs in codons 24-87
    # (Peixeiro 2011); PTCs 5' of p.Glu23 or 3' of p.Thr88 escape NMD. Escape
    # truncations removing >10% → Strong, <10% → Moderate (so early PTCs that
    # remove most of the protein → Strong). Initiation codon → PVS1.
    "HBB": _GeneSpec(
        gene="HBB", transcript="NM_000518.5", aa_len=147,
        trunc_bands=((1, 23, _S.STRONG), (24, 87, _S.VERY_STRONG),
                     (88, 133, _S.STRONG), (134, 147, _S.MODERATE)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
    # HBA2 (Hemoglobinopathy VCEP, α-globin). NMD for PTCs 5' of p.Leu84 (c.251);
    # at/3' of p.Leu84 escape NMD → Strong if >10% removed else Moderate.
    # Initiation codon → PVS1.
    "HBA2": _GeneSpec(
        gene="HBA2", transcript="NM_000517.6", aa_len=142,
        trunc_bands=((1, 83, _S.VERY_STRONG), (84, 128, _S.STRONG),
                     (129, 142, _S.MODERATE)),
        start_lost=_S.VERY_STRONG, splice=_S.VERY_STRONG, deletion=_S.VERY_STRONG,
    ),
}


def _intron_number(intron: str | None) -> int | None:
    """Parse the intron index from VEP's ``n/N`` field (e.g. "8/19" -> 8)."""
    if not intron:
        return None
    try:
        return int(intron.split("/")[0])
    except (ValueError, IndexError):
        return None


def _escape_10pct(codon: int | None, aa_len: int) -> CriterionStrength:
    """NMD-escape strength by the generic 10%-of-protein rule: Strong if the
    truncation removes >10% of the protein, else Moderate. A missing codon is
    treated conservatively as a large (>10%) loss → Strong."""
    if codon is None:
        return _S.STRONG
    removed = aa_len - codon + 1
    return _S.STRONG if removed > 0.10 * aa_len else _S.MODERATE


def evaluate_vcep_pvs1(pc, splice_overrides=None) -> tuple[CriterionStrength, str] | None:
    """Gene-specific PVS1 strength, or ``None`` when this gene has no VCEP
    deviation for the variant's consequence (the caller then falls back to the
    generic tree). ``CriterionStrength.NOT_MET`` is an explicit VCEP N/A.

    ``splice_overrides`` is an optional
    :class:`~acmg_classifier.pvs1.vcep_pvs1_exons.SpliceExonOverrides`; when it
    has a reviewer-supplied strength for the exon a canonical splice variant
    skips, that strength replaces the flat per-gene splice default."""
    spec = _SPECS.get(pc.gene_symbol)
    if spec is None:
        return None

    cq = pc.consequence
    if cq in _TRUNC:
        # NMD-based truncation rule (FBN1, CDH1, ACADVL, GAMT): the cutoff is the
        # exon-based NMD boundary, not a fixed codon.
        if spec.trunc_nmd is not None:
            from acmg_classifier.pvs1.nmd_predictor import predicts_nmd
            nmd = predicts_nmd(pc)
            if nmd:
                strength = spec.trunc_nmd[0]
                note = "NMD predicted"
            elif spec.trunc_nmd[1] is not None:
                strength = spec.trunc_nmd[1]
                note = "NMD escape (critical C-terminus)"
            else:
                # Escape strength by the generic 10%-of-protein rule.
                strength = _escape_10pct(pc.protein_position, spec.aa_len)
                note = "NMD escape (10% rule)"
            return strength, (
                f"{spec.gene} VCEP: {cq.value}; {note} -> {_LABEL[strength]}"
            )

        codon = pc.protein_position
        if codon is None:
            return None  # cannot place the truncation → defer to generic tree
        # Frameshift uses its own bands when the VCEP scores it differently.
        bands = spec.trunc_bands
        if cq == ConsequenceType.FRAMESHIFT and spec.trunc_bands_fs is not None:
            bands = spec.trunc_bands_fs
        for lo, hi, strength in bands or ():
            if lo <= codon <= hi:
                return strength, (
                    f"{spec.gene} VCEP: truncation at codon {codon} "
                    f"({spec.transcript}) -> {_LABEL[strength]}"
                )
        return _S.NOT_MET, (
            f"{spec.gene} VCEP: truncation at codon {codon} outside PVS1 "
            f"range -> N/A"
        )

    if cq == ConsequenceType.START_LOST:
        if spec.start_lost is None:
            return None
        return spec.start_lost, (
            f"{spec.gene} VCEP: initiation-codon variant -> {_LABEL[spec.start_lost]}"
        )

    if cq in _SPLICE:
        if spec.splice is None:
            return None
        # Donor splice sites of specific introns may be excluded (ACADVL intron 8
        # is a GC donor of unclear impact → no PVS1).
        if cq == ConsequenceType.SPLICE_DONOR and spec.splice_exclude_donor_introns:
            intron_no = _intron_number(pc.intron)
            if intron_no is not None and intron_no in spec.splice_exclude_donor_introns:
                return _S.NOT_MET, (
                    f"{spec.gene} VCEP: donor of intron {intron_no} excluded "
                    f"(GC donor / unclear impact) -> N/A"
                )
        # Optional exon-aware refinement: if a reviewer supplied a strength for
        # the exon this variant skips, it overrides the flat default.
        if splice_overrides:
            from acmg_classifier.pvs1.vcep_pvs1_exons import skipped_exon
            exon = skipped_exon(cq, pc.intron)
            override = splice_overrides.lookup(spec.gene, exon)
            if override is not None:
                return override, (
                    f"{spec.gene} VCEP: canonical +/-1,2 splice skipping exon "
                    f"{exon} -> {_LABEL[override]} (exon-aware)"
                )
        return spec.splice, (
            f"{spec.gene} VCEP: canonical +/-1,2 splice variant -> {_LABEL[spec.splice]}"
        )

    if cq == ConsequenceType.TRANSCRIPT_ABLATION:
        if spec.deletion is None:
            return None
        return spec.deletion, (
            f"{spec.gene} VCEP: whole-gene / exon deletion -> {_LABEL[spec.deletion]}"
        )

    return None
