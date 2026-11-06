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

Group C (point table)
  * SCID VCEP (ADA, DCLRE1C, FOXN1, IL2RG, IL7R, JAK3, RAG1, RAG2): nonsense
    variants scored from informative VCEP-classified truncating variants
    (:func:`evaluate_scid_pm5`): NMD-predicted VUA -> each P/LP truncating variant
    in the same exon predicted to lead to NMD +1, each B/LB one -2; VUA with a PTC
    in the final exon (no NMD) -> each P/LP truncating variant with its PTC in the
    same exon downstream of the VUA +1, each B/LB one upstream -2. >= 4 points
    Strong, >= 2 Moderate, 1 Supporting; <= 0 not applied. The registry withholds
    it without PM2 and downgrades it with PVS1 (Strong -> Moderate with any PVS1,
    Moderate -> Supporting with PVS1_VeryStrong).

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
SCID_GENES = frozenset({"ADA", "DCLRE1C", "FOXN1", "IL2RG", "IL7R", "JAK3", "RAG1", "RAG2"})
SCID_EVIDENCE_PREFIX = "PM5_PTC (SCID"
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


# ---------------------------------------------------------------------------
# Group C: SCID VCEP point table (nonsense variants)
# ---------------------------------------------------------------------------

_CPOS_ANY = re.compile(r"c\.(\d+)")
_ACC = re.compile(r"^(N[MR]_\d+|ENST\d+)")


def _record_codons(hgvs_p: Optional[str]) -> tuple[Optional[int], Optional[int]]:
    """(first altered codon, PTC codon) of a ClinVar nonsense / frameshift record."""
    from acmg_classifier.pvs1.nmd_predictor import _FS_RE, _STOP_RE
    hp = (hgvs_p or "").split(":")[-1].strip()
    m = _FS_RE.search(hp)
    if m:
        return int(m.group(1)), int(m.group(1)) + int(m.group(2)) - 1
    m = _STOP_RE.search(hp)
    if m:
        return int(m.group(1)), int(m.group(1))
    return None, None


def scid_points(pc: ConsequenceInfo, records, st, exclude=None) -> tuple[int, list[str], str]:
    """(points, per-variant notes, rule) of the SCID PM5 point table for the
    nonsense VUA ``pc`` given ClinVar truncating ``records`` and the transcript
    structure ``st``; rule is "nmd" / "final_exon" / "" (not covered)."""
    from acmg_classifier.pvs1.nmd_predictor import NMD_JUNCTION_NT, ptc_codon
    codon = ptc_codon(pc)
    if codon is None:
        return 0, [], ""
    vua_end = st.cds_start - 1 + codon * 3
    vua_exon = st.exon_of(vua_end)
    vua_nmd = st.n_exons > 1 and st.last_junction - vua_end > NMD_JUNCTION_NT
    if not vua_nmd and vua_exon != st.n_exons:
        return 0, [], ""            # NMD escape in the penultimate exon: not in the table
    rule = "nmd" if vua_nmd else "final_exon"
    pts, notes, seen = 0, [], set()
    for h in records:
        acc = _ACC.match(h.hgvs_c or "")
        if acc and _base(acc.group(1)) != _base(pc.transcript_id):
            continue                # named on another transcript: cannot place it
        first, ptc = _record_codons(h.hgvs_p)
        key = (h.chrom, h.pos, h.ref, h.alt)
        if ptc is None or key in seen or (exclude is not None and key == exclude):
            continue
        seen.add(key)
        if ptc == codon and (h.hgvs_p or "").split(":")[-1].endswith(("Ter", "*", "X")):
            continue                # the VUA's own change
        sig = (h.clinical_significance or "").lower()
        plp = "pathogenic" in sig
        ptc_end = st.cds_start - 1 + ptc * 3
        if rule == "nmd":
            cm = _CPOS_ANY.search((h.hgvs_c or "").split(":")[-1])
            change_pos = (st.cds_start - 1 + int(cm.group(1))) if cm else \
                (st.cds_start - 1 + first * 3 - 2)
            same_exon = st.exon_of(change_pos) == vua_exon
            nmd = st.last_junction - ptc_end > NMD_JUNCTION_NT
            if not (same_exon and nmd):
                continue
            score = 1 if plp else -2
        else:
            if st.exon_of(ptc_end) != vua_exon:
                continue
            if plp and ptc > codon:
                score = 1
            elif not plp and ptc < codon:
                score = -2
            else:
                continue
        pts += score
        notes.append(f"{h.variation_id or h.hgvs_p} {h.clinical_significance} {score:+d}")
    return pts, notes, rule


def evaluate_scid_pm5(pc: ConsequenceInfo, cfg, variant
                      ) -> Optional[tuple[CriterionStrength, str]]:
    """SCID VCEP PM5 point table for a nonsense variant, or None when the gene /
    consequence is not covered."""
    if pc.gene_symbol not in SCID_GENES or pc.consequence != ConsequenceType.STOP_GAINED:
        return None
    from acmg_classifier.criteria.vcep_comparators import filter_comparators
    from acmg_classifier.local_db.clinvar_sqlite import query_truncating_variants
    from acmg_classifier.pvs1.nmd_predictor import transcript_structure
    gene = pc.gene_symbol
    tag = f"{SCID_EVIDENCE_PREFIX} {gene} point table)"
    st = transcript_structure(pc.transcript_id)
    if st is None:
        return _S.NOT_MET, f"{tag}: transcript structure of {pc.transcript_id} unavailable"
    min_stars = getattr(cfg, "pm5_min_stars", 2)
    if not isinstance(min_stars, int):
        min_stars = 2
    records = query_truncating_variants(cfg.clinvar_sqlite, gene, min_stars)
    # Informative variants must be classified by the SCID VCEP (eRepo /
    # expert-panel ClinVar record), as for the missense PM5 comparators.
    records, vcep_note = filter_comparators(gene, "PM5", records, variant.assembly.value)
    chrom = variant.chrom[3:] if variant.chrom.lower().startswith("chr") else variant.chrom
    pts, notes, rule = scid_points(pc, records, st,
                                   exclude=(chrom, variant.pos, variant.ref, variant.alt))
    if not rule:
        return _S.NOT_MET, f"{tag}: NMD escape outside the final exon - not in the point table"
    detail = (f"{'; '.join(notes[:6])}" if notes else "no informative variant") + \
        (f" [{vcep_note}]" if vcep_note else "")
    where = "NMD-predicted exon" if rule == "nmd" else "final exon (no NMD)"
    if pts <= 0:
        extra = " - negative points: reconsider PVS1" if pts < 0 else ""
        return _S.NOT_MET, f"{tag}: {pts} point(s), {where}{extra} ({detail})"
    strength = _S.STRONG if pts >= 4 else _S.MODERATE if pts >= 2 else _S.SUPPORTING
    return strength, (f"{tag}: {pts} point(s), {where} -> {strength.value} ({detail}); "
                      "informative truncating variants must be P/LP without PM5 and not "
                      "only PVS1 + PM2 - review")
