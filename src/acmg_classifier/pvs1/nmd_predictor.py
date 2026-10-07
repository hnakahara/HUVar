"""NMD (nonsense-mediated decay) prediction for the ClinGen PVS1 decision tree.

Two levels of precision:

1. **PTC-based (preferred).** When the variant's transcript is a MANE Select /
   MANE Plus Clinical transcript with a known exon structure
   (``resources/shared/mane_exons.tsv``, built from the MANE v1.5 GFF by
   ``scripts/build_mane_exons.py``), the position of the premature termination
   codon (PTC) is located in cDNA coordinates and NMD is predicted when the PTC
   lies more than 50 nucleotides upstream of the last exon-exon junction
   (Abou Tayoun et al. 2018). The PTC is the variant codon for stop-gained
   variants and, for frameshifts, the new stop codon given by the HGVS protein
   notation (``p.Arg100GlyfsTer12`` -> codon 100 + 12 - 1 = 111), so a
   frameshift whose new stop falls in the last exon is correctly predicted to
   escape NMD even though the variant itself lies in an upstream exon.

2. **Exon-number fallback.** When the transcript structure or the PTC position
   is unavailable, the VEP exon number of the variant is used: the last exon
   (or a single-exon transcript) escapes NMD, any other exon is predicted to
   undergo NMD. Missing/malformed exon information defaults to "NMD predicted".
"""
from __future__ import annotations

import csv
import re
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple, Optional

from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType

# Distance (nt) upstream of the last exon-exon junction beyond which a PTC
# triggers NMD (the "50-55 nt rule"; Abou Tayoun et al. 2018).
NMD_JUNCTION_NT = 50

# parents: [0]=pvs1 [1]=acmg_classifier [2]=src [3]=repo
_PACKAGED_TSV = Path(__file__).resolve().parents[3] / "resources" / "shared" / "mane_exons.tsv"


class TxStructure(NamedTuple):
    accession: str
    cds_start: int              # 1-based cDNA position of the first coding nt
    cds_end: int                # 1-based cDNA position of the last coding nt (stop incl.)
    exon_lengths: tuple[int, ...]

    @property
    def n_exons(self) -> int:
        return len(self.exon_lengths)

    @property
    def last_junction(self) -> int:
        """cDNA position of the last nucleotide before the last exon-exon junction."""
        return sum(self.exon_lengths[:-1])

    def exon_of(self, cdna_pos: int) -> Optional[int]:
        acc = 0
        for i, ln in enumerate(self.exon_lengths, start=1):
            acc += ln
            if cdna_pos <= acc:
                return i
        return None


def _base(acc: Optional[str]) -> str:
    return (acc or "").split(".")[0]


@lru_cache(maxsize=4)
def load_structures(tsv_path: Path = _PACKAGED_TSV) -> dict[str, TxStructure]:
    """``{version-stripped accession: TxStructure}`` keyed by both RefSeq and Ensembl IDs."""
    out: dict[str, TxStructure] = {}
    if not tsv_path.exists():
        return out
    with tsv_path.open(encoding="utf-8") as fh:
        rows = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for r in rows:
            try:
                st = TxStructure(
                    accession=r["refseq"],
                    cds_start=int(r["cds_start"]),
                    cds_end=int(r["cds_end"]),
                    exon_lengths=tuple(int(x) for x in r["exon_lengths"].split(",") if x),
                )
            except (KeyError, ValueError):
                continue
            for acc in (r.get("refseq"), r.get("ensembl")):
                if acc:
                    out.setdefault(_base(acc), st)
    return out


def transcript_structure(transcript_id: Optional[str]) -> Optional[TxStructure]:
    return load_structures().get(_base(transcript_id)) if transcript_id else None


# p.Arg100GlyfsTer12 / p.(Arg100Glyfs*12) / p.R100Gfs*12 / p.Arg100Ter / p.(Arg100*)
_FS_RE = re.compile(r"p\.\(?[A-Z][a-z]{0,2}(\d+)(?:[A-Z][a-z]{0,2})?fs(?:Ter|\*|X)(\d+)")
_STOP_RE = re.compile(r"p\.\(?[A-Z][a-z]{0,2}(\d+)(?:Ter|\*|X)\)?$")


def ptc_codon(pc: ConsequenceInfo) -> Optional[int]:
    """Codon number of the premature termination codon, or None when unknown.

    Stop-gained: the variant codon. Frameshift: the new stop codon from HGVS
    ``fsTerN`` (N counts the first altered residue as 1). ``fsTer?`` (no stop in
    the transcript) returns None."""
    hp = (pc.hgvs_p or "").split(":")[-1].replace("%3D", "=").strip()
    if pc.consequence == ConsequenceType.FRAMESHIFT:
        m = _FS_RE.search(hp)
        if m:
            return int(m.group(1)) + int(m.group(2)) - 1
        return None
    if pc.consequence == ConsequenceType.STOP_GAINED:
        m = _STOP_RE.search(hp)
        if m:
            return int(m.group(1))
        return pc.protein_position
    return None


def nmd_status(pc: ConsequenceInfo) -> tuple[bool, str]:
    """``(nmd_predicted, explanation)`` for a truncating variant."""
    st = transcript_structure(pc.transcript_id)
    codon = ptc_codon(pc)
    if st is not None and codon is not None and st.n_exons >= 1:
        # cDNA position of the last nucleotide of the PTC.
        ptc_end = st.cds_start - 1 + codon * 3
        if ptc_end <= st.cds_end and st.n_exons > 1:
            dist = st.last_junction - ptc_end
            exon = st.exon_of(ptc_end)
            where = f"PTC codon {codon} in exon {exon}/{st.n_exons} ({st.accession})"
            if dist > NMD_JUNCTION_NT:
                return True, f"{where}, {dist} nt upstream of the last exon-exon junction -> NMD predicted"
            if dist > 0:
                return False, (f"{where}, within the 3'-most {NMD_JUNCTION_NT} nt of the "
                               "penultimate exon -> NMD escape")
            return False, f"{where} (last exon) -> NMD escape"
        if st.n_exons == 1:
            return False, f"single-exon transcript ({st.accession}) -> NMD escape"
    return _exon_rule(pc)


def _exon_rule(pc: ConsequenceInfo) -> tuple[bool, str]:
    exon_str = pc.exon
    if exon_str is None:
        return True, "exon unknown -> NMD assumed"
    parts = exon_str.split("/")
    if len(parts) != 2:
        return True, "exon unknown -> NMD assumed"
    try:
        exon_num, total = int(parts[0]), int(parts[1])
    except ValueError:
        return True, "exon unknown -> NMD assumed"
    if total == 1:
        return False, "single-exon transcript -> NMD escape (exon-number rule)"
    if exon_num == total:
        return False, f"variant in last exon {exon_str} -> NMD escape (exon-number rule)"
    return True, f"variant in exon {exon_str} -> NMD predicted (exon-number rule)"


def predicts_nmd(consequence: ConsequenceInfo) -> bool:
    """Return True if the variant is predicted to trigger NMD (see module docstring)."""
    return nmd_status(consequence)[0]


def is_last_exon(consequence: ConsequenceInfo) -> bool:
    """Strict "exon == last exon" test parsed from VEP's "n/N" field."""
    exon_str = consequence.exon
    if exon_str is None:
        return False
    parts = exon_str.split("/")
    if len(parts) != 2:
        return False
    try:
        return int(parts[0]) == int(parts[1])
    except ValueError:
        return False


def is_penultimate_exon(consequence: ConsequenceInfo) -> bool:
    """Penultimate-exon test (exon == last-1) on the VEP "n/N" field."""
    exon_str = consequence.exon
    if exon_str is None:
        return False
    parts = exon_str.split("/")
    if len(parts) != 2:
        return False
    try:
        return int(parts[0]) == int(parts[1]) - 1
    except ValueError:
        return False
