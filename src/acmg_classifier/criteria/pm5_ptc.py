"""PM5 for protein-truncating variants (PM5_PTC), from ClinGen VCEP specifications.

Several VCEPs extend PM5 ("novel change at a residue where a different
pathogenic change has been seen") from missense to nonsense / frameshift
variants: a premature termination codon (PTC) located where other pathogenic
PTC variants are established is given additional weight. HUVar implements the
VCEPs whose rule can be evaluated from the variant itself (position of the
variant or of its PTC, and NMD prediction):

Group A (position rules, PM5_Supporting)
  * RUNX1 (Myeloid Malignancy VCEP, GN008 v3.1.0): nonsense/frameshift variants
    downstream of c.98 (NM_001754.4). The variant must not impact splicing
    (SpliceAI <= 0.20); PM5 cannot be used with PM1 (``pm5_excludes`` column of
    ``disease_prevalence.tsv``, applied by the registry).
  * CDH1 (HDGC VCEP, GN007 v3.1): nonsense/frameshift variants predicted to
    undergo NMD, or whose PTC is located upstream of the last known pathogenic
    truncating variant. The latter is the most 3' truncating variant the CDH1
    VCEP classified Pathogenic in the ClinGen Evidence Repository,
    NM_004360.5:c.2506G>T p.(Glu836Ter); a PTC at codon < 836 qualifies. The
    site-specific canonical-splice recommendations are not automated.
  * ATM (HBOP VCEP, GN020 v1.5.0): PTC upstream of p.Arg3047.
  * PALB2 (HBOP VCEP, GN077 v1.2.0): PTC upstream of p.Tyr1183.
    For ATM and PALB2, splice variants qualify only with RNA evidence
    (PVS1_VS(RNA)) and are not automated.

Group B (exon table)
  * BRCA1 / BRCA2 (ENIGMA VCEP, GN092 / GN097 v1.2.0): PM5_PTC strength by the
    exon in which the termination codon occurs (Specifications Table 4,
    V1.2 2024-11-18; MANE Select NM_007294.4 / NM_000059.4 exon numbering).
    Applied only to genomic PTC changes (nonsense / frameshift, not splicing)
    and only to PTC variants annotated as PVS1 — the registry suppresses
    PM5_PTC for BRCA1/2 when PVS1 is not met.

Not implemented: the SCID VCEP point system (ADA, DCLRE1C, FOXN1, IL2RG, IL7R,
JAK3, RAG1, RAG2), which scores VCEP-classified downstream nonsense variants.

The PTC codon is the variant codon for stop-gained variants and the new stop
codon for frameshifts (``p.Arg100GlyfsTer12`` -> codon 111), as computed by
:func:`acmg_classifier.pvs1.nmd_predictor.ptc_codon`. Frameshifts without a
stop in the transcript (``fsTer?``) are not evaluated.
"""
from __future__ import annotations

import re
from typing import Optional

from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType, CriterionStrength

_S = CriterionStrength
PTC_CONSEQUENCES = {ConsequenceType.STOP_GAINED, ConsequenceType.FRAMESHIFT}

# Transcript the VCEP rule is written for (version-stripped RefSeq accession).
_TRANSCRIPT = {
    "RUNX1": "NM_001754",
    "CDH1": "NM_004360",
    "ATM": "NM_000051",
    "PALB2": "NM_024675",
    "BRCA1": "NM_007294",
    "BRCA2": "NM_000059",
}

# Group A: PTC codon must be strictly upstream of this residue.
_PTC_UPSTREAM_OF = {"ATM": 3047, "PALB2": 1183}
# CDH1: most 3' Pathogenic truncating variant classified by the CDH1 VCEP in the
# ClinGen Evidence Repository (c.2506G>T, p.Glu836Ter). No Pathogenic frameshift
# in the eRepo has a PTC further 3' (latest: c.2430del p.Phe810LeufsTer6, PTC 815).
_CDH1_LAST_PATHOGENIC_PTC = 836
_RUNX1_MIN_CPOS = 98          # variant must lie downstream of c.98
_RUNX1_MAX_SPLICE = 0.20      # SpliceAI max delta must be <= 0.20

# Group B: ENIGMA Specifications Table 4 (V1.2, 2024-11-18), PM5_PTC per exon.
# (exon, last coding c. position of the exon, strength); None = not applicable.
# The last exon is split at a codon cutoff (see _BRCA_LAST_EXON).
_BRCA_EXONS: dict[str, tuple[tuple[int, int, Optional[CriterionStrength]], ...]] = {
    "BRCA1": (
        (2, 80, _S.STRONG), (3, 134, _S.STRONG), (4, 212, _S.STRONG),
        (5, 301, _S.STRONG), (6, 441, _S.STRONG), (7, 547, _S.SUPPORTING),
        (8, 593, None), (9, 670, None), (10, 4096, _S.STRONG),
        (11, 4185, _S.STRONG), (12, 4357, _S.STRONG), (13, 4484, _S.STRONG),
        (14, 4675, _S.STRONG), (15, 4986, _S.STRONG), (16, 5074, _S.STRONG),
        (17, 5152, _S.STRONG), (18, 5193, _S.STRONG), (19, 5277, _S.STRONG),
        (20, 5332, _S.STRONG), (21, 5406, _S.STRONG), (22, 5467, _S.STRONG),
    ),
    "BRCA2": (
        (2, 67, _S.STRONG), (3, 316, _S.STRONG), (4, 425, _S.MODERATE),
        (5, 475, _S.STRONG), (6, 516, None), (7, 631, _S.STRONG),
        (8, 681, _S.STRONG), (9, 793, _S.STRONG), (10, 1909, _S.STRONG),
        (11, 6841, _S.STRONG), (12, 6937, None), (13, 7007, _S.STRONG),
        (14, 7435, _S.STRONG), (15, 7617, _S.STRONG), (16, 7805, _S.STRONG),
        (17, 7976, _S.STRONG), (18, 8331, _S.STRONG), (19, 8487, _S.STRONG),
        (20, 8632, _S.STRONG), (21, 8754, _S.SUPPORTING), (22, 8953, _S.STRONG),
        (23, 9117, _S.STRONG), (24, 9256, _S.STRONG), (25, 9501, _S.STRONG),
        (26, 9648, _S.SUPPORTING),
    ),
}
# Last exon: (exon number, last codon given PM5_PTC_Strong); PTCs at later
# codons are not applicable (BRCA1 exon 23: PTC < p.I1855; BRCA2 exon 27:
# PTC < p.T3310).
_BRCA_LAST_EXON = {"BRCA1": (23, 1854), "BRCA2": (27, 3309)}

GENES = frozenset(_TRANSCRIPT)
REQUIRES_PVS1 = frozenset({"BRCA1", "BRCA2"})

_CPOS = re.compile(r"c\.(\d+)")


def _base(acc: Optional[str]) -> str:
    return (acc or "").split(".")[0]


def requires_pvs1(gene: Optional[str]) -> bool:
    """True when the VCEP restricts PM5_PTC to variants annotated as PVS1."""
    return bool(gene) and gene in REQUIRES_PVS1


def is_cdh1_canonical_splice(pc: Optional[ConsequenceInfo]) -> bool:
    return (pc is not None and pc.gene_symbol == "CDH1" and pc.consequence in (
        ConsequenceType.SPLICE_DONOR, ConsequenceType.SPLICE_ACCEPTOR))


def is_ptc_variant(pc: Optional[ConsequenceInfo]) -> bool:
    return pc is not None and pc.consequence in PTC_CONSEQUENCES


def _variant_cpos(pc: ConsequenceInfo) -> Optional[int]:
    m = _CPOS.search((pc.hgvs_c or "").split(":")[-1])
    if m:
        return int(m.group(1))
    if pc.protein_position:
        return pc.protein_position * 3 - 2
    return None


def brca_exon_strength(gene: str, codon: int) -> tuple[int, Optional[CriterionStrength]]:
    """``(exon, strength)`` for a PTC at *codon* in BRCA1/BRCA2 (strength None =
    PM5_PTC not applicable). The exon is the one holding the first nucleotide
    of the termination codon."""
    cpos = codon * 3 - 2
    for exon, end, strength in _BRCA_EXONS[gene]:
        if cpos <= end:
            return exon, strength
    last_exon, last_codon = _BRCA_LAST_EXON[gene]
    return last_exon, (_S.STRONG if codon <= last_codon else None)


def evaluate_pm5_ptc(
    pc: ConsequenceInfo, annotation: AnnotationData,
) -> Optional[tuple[CriterionStrength, str]]:
    """PM5_PTC for a nonsense/frameshift variant.

    Returns ``None`` when the gene has no PTC-based PM5 rule (the caller keeps
    its usual missense-only behaviour), otherwise ``(strength, evidence)`` with
    ``CriterionStrength.NOT_MET`` when the rule does not apply."""
    gene = pc.gene_symbol
    cdh1_splice = gene == "CDH1" and is_cdh1_canonical_splice(pc)
    if gene not in GENES or not (is_ptc_variant(pc) or cdh1_splice):
        return None
    tx = _TRANSCRIPT[gene]
    if _base(pc.transcript_id) != tx:
        return _S.NOT_MET, f"PM5_PTC ({gene}): rule defined on {tx}, variant on {pc.transcript_id}"
    if cdh1_splice:
        # CDH1 v3.1: PM5_Supporting also for canonical splice variants (site-
        # specific splicing table). In the eRepo it accompanies PVS1_Strong
        # sites and is withheld at the in-frame PVS1_Moderate sites; the
        # registry suppresses it unless PVS1 is met at Strong or higher.
        return _S.SUPPORTING, ("PM5_PTC (CDH1, HDGC VCEP): canonical +/-1,2 splice "
                               "variant predicted to truncate -> PM5_Supporting")

    from acmg_classifier.pvs1.nmd_predictor import nmd_status, ptc_codon

    if gene == "RUNX1":
        cpos = _variant_cpos(pc)
        if cpos is None:
            return _S.NOT_MET, "PM5_PTC (RUNX1): variant c. position unknown"
        if cpos <= _RUNX1_MIN_CPOS:
            return _S.NOT_MET, (f"PM5_PTC (RUNX1): c.{cpos} not downstream of c.98 "
                                "(affects only one functional transcript)")
        sp = annotation.splice
        if sp is not None and sp.is_available and sp.max_delta is not None \
                and sp.max_delta > _RUNX1_MAX_SPLICE:
            return _S.NOT_MET, (f"PM5_PTC (RUNX1): predicted splice impact "
                                f"({sp.tool} max_delta={sp.max_delta:.2f} > 0.20)")
        return _S.SUPPORTING, (f"PM5_PTC (RUNX1, Myeloid Malignancy VCEP): "
                               f"nonsense/frameshift downstream of c.98 (c.{cpos})")

    if gene == "CDH1":
        nmd, why = nmd_status(pc)
        if nmd:
            return _S.SUPPORTING, (f"PM5_PTC (CDH1, HDGC VCEP): nonsense/frameshift "
                                   f"predicted NMD ({why})")
        codon = ptc_codon(pc)
        if codon is not None and codon < _CDH1_LAST_PATHOGENIC_PTC:
            return _S.SUPPORTING, (f"PM5_PTC (CDH1, HDGC VCEP): NMD escape but PTC codon "
                                   f"{codon} upstream of the last known pathogenic "
                                   "truncating variant p.Glu836Ter")
        return _S.NOT_MET, (f"PM5_PTC (CDH1): not predicted to undergo NMD ({why}) and "
                            "PTC not upstream of p.Glu836Ter")

    codon = ptc_codon(pc)
    if codon is None:
        return _S.NOT_MET, f"PM5_PTC ({gene}): PTC position unknown"

    if gene in _PTC_UPSTREAM_OF:
        limit = _PTC_UPSTREAM_OF[gene]
        aa = "p.Arg3047" if gene == "ATM" else "p.Tyr1183"
        if codon >= limit:
            return _S.NOT_MET, f"PM5_PTC ({gene}): PTC codon {codon} not upstream of {aa}"
        return _S.SUPPORTING, (f"PM5_PTC ({gene}, HBOP VCEP): PTC codon {codon} "
                               f"upstream of {aa}")

    # BRCA1 / BRCA2
    exon, strength = brca_exon_strength(gene, codon)
    if strength is None:
        return _S.NOT_MET, (f"PM5_PTC ({gene}): PTC codon {codon} in exon {exon} "
                            "— not applicable (ENIGMA Table 4)")
    return strength, (f"PM5_PTC ({gene}, ENIGMA Table 4): PTC codon {codon} in exon "
                      f"{exon} -> {strength.value}")
