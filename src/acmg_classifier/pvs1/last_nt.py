"""VCEP "G to non-G at the last nucleotide of an exon" PVS1 rules.

Exonic changes at the last base of an exon weaken the 5' splice site; three
VCEPs score them under PVS1 even though VEP calls them missense / synonymous:

* **ATM** (HBOP VCEP, Specifications v1.6, PVS1 guide): a G>non-G change at the
  last nucleotide of an exon, when the adjacent intronic sequence is not
  ``gtrrgt`` (r = purine), gets the weight the decision tree gives the
  canonical +1/+2 donor of that intron, *notched one level down*. The donor
  weight follows the tree: exon skipping that disrupts the reading frame ->
  PVS1 (lists A / D); in-frame skipping in the N-terminal HEAT-repeat solenoid
  (p.1-1892, exons 2-38) -> PVS1_Strong (list B); in-frame skipping in FATKIN
  (p.1893-3056) -> PVS1 (list E). Hence last-nt changes are PVS1_Strong, or
  PVS1_Moderate for an in-frame HEAT-repeat exon.
* **MLH1** (InSiGHT VCEP): "G>non-G at last base of exon if first 6 bases of the
  intron are not GTRRGT" -> PVS1_Strong (PVS1 only with a proven splice defect,
  which is RNA evidence and not modelled here).
* **CDH1** (CDH1 VCEP v3.1): G to non-G at the last nucleotide of an exon ->
  PVS1_Moderate.

The exon boundary comes from the packaged MANE exon structures; the intronic
+1..+6 sequence is read from the reference FASTA when one is available (the
``gtrrgt`` exclusion is applied only for ATM and MLH1, whose VCEPs state it).
Without a FASTA the rule is applied and the evidence says the intron was not
checked.
"""
from __future__ import annotations

import itertools
import re
from typing import Callable, Optional

from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import CriterionStrength as _S

_COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
_SNV_RE = re.compile(r"c\.(\d+)([ACGT])>([ACGT])$")
_GTRRGT = re.compile(r"^GT[AG][AG]GT$")

# gene -> (rule, check intronic gtrrgt)
_GENES: dict[str, tuple[str, bool]] = {
    "ATM": ("atm", True),
    "MLH1": ("strong", True),
    "CDH1": ("moderate", False),
}
_ATM_HEAT_MAX = 1892  # last residue of the N-terminal solenoid (HEAT repeats)

_LABEL = {_S.VERY_STRONG: "PVS1", _S.STRONG: "PVS1_Strong",
          _S.MODERATE: "PVS1_Moderate", _S.SUPPORTING: "PVS1_Supporting"}

# (chrom, start, end) 1-based inclusive -> genomic + strand sequence
SeqFetcher = Callable[[str, int, int], Optional[str]]


def _exon_bounds(transcript_id: str) -> Optional[list[tuple[int, int]]]:
    """Coding-exon spans as (first c., last c.) for every exon but the last."""
    from acmg_classifier.pvs1.nmd_predictor import transcript_structure
    st = transcript_structure(transcript_id)
    if st is None:
        return None
    ends = list(itertools.accumulate(st.exon_lengths))
    out, start = [], 1
    for e in ends:
        out.append((start - st.cds_start + 1, e - st.cds_start + 1))
        start = e + 1
    return out[:-1]  # the last exon has no donor


def evaluate_last_nt(
    pc: ConsequenceInfo,
    chrom: Optional[str] = None,
    pos: Optional[int] = None,
    vcf_ref: Optional[str] = None,
    fetch: Optional[SeqFetcher] = None,
) -> Optional[tuple[_S, str]]:
    """``(strength, evidence)`` for a VCEP last-nucleotide G>non-G variant, or
    ``None`` when the rule does not apply (other gene / change / position, or
    a ``gtrrgt`` intron)."""
    rule = _GENES.get(pc.gene_symbol)
    if rule is None or not pc.hgvs_c:
        return None
    m = _SNV_RE.search(pc.hgvs_c.split(":")[-1])
    if not m or m.group(2) != "G" or m.group(3) == "G":
        return None
    cpos = int(m.group(1))
    bounds = _exon_bounds(pc.transcript_id)
    if not bounds:
        return None
    exon = next(((i, s, e) for i, (s, e) in enumerate(bounds, start=1) if e == cpos), None)
    if exon is None:
        return None
    exon_no, ex_start, ex_end = exon
    kind, check_intron = rule

    intron_note = "intron +1..+6 not checked (no reference FASTA)"
    if check_intron and fetch is not None and chrom and pos and vcf_ref:
        if vcf_ref.upper() == "G":          # transcript on the + strand
            seq = fetch(chrom, pos + 1, pos + 6)
        elif vcf_ref.upper() == "C":        # transcript on the - strand
            seq = fetch(chrom, pos - 6, pos - 1)
            seq = seq.translate(_COMP)[::-1] if seq else seq
        else:
            seq = None
        if seq:
            seq = seq.upper()
            if _GTRRGT.match(seq):
                return None  # strong donor consensus: VCEP rule not met
            intron_note = f"intron +1..+6 {seq} (not GTRRGT)"
    elif not check_intron:
        intron_note = ""

    if kind == "atm":
        in_frame = (ex_end - ex_start + 1) % 3 == 0 and ex_start >= 1
        heat = (ex_end + 2) // 3 <= _ATM_HEAT_MAX
        if in_frame and heat:
            strength = _S.MODERATE
            donor = "in-frame HEAT-repeat exon skipping (donor PVS1_Strong)"
        else:
            strength = _S.STRONG
            donor = ("in-frame FATKIN exon skipping (donor PVS1)" if in_frame
                     else "frameshifting exon skipping (donor PVS1)")
        why = f"{donor}, one level down"
    elif kind == "strong":
        strength, why = _S.STRONG, "G>non-G at the last exon base"
    else:
        strength, why = _S.MODERATE, "G>non-G at the last exon base"

    note = f"; {intron_note}" if intron_note else ""
    return strength, (
        f"{pc.gene_symbol} VCEP: c.{cpos}G>{m.group(3)} at the last nucleotide of "
        f"exon {exon_no}; {why}{note} -> {_LABEL[strength]}"
    )
