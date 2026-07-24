"""MANE Select fallback for caches without the flag (GRCh37).

VEP only emits ``mane_select`` in the GRCh38 cache. On GRCh37 every transcript
comes back unflagged, so the MANE-first ordering in ``_parse_vep_record`` would
otherwise pick a non-MANE isoform (e.g. PTEN NM_001304717 vs MANE NM_000314),
shifting HGVS numbering and breaking VCEP codon-range criteria. The fallback
recovers the flag by matching the gene's MANE base accession.
"""
from acmg_classifier.local_db.vep_runner import _apply_mane_fallback, _parse_vep_record
from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType

# PTEN: MANE Select NM_000314 / ENST00000371953; NM_001304717 is the longer isoform.
_PTEN_MANE = {"PTEN": ("NM_000314", "ENST00000371953")}


def _cons(tx: str, mane: bool = False, canonical: bool = False) -> ConsequenceInfo:
    return ConsequenceInfo(
        transcript_id=tx, gene_id="g", gene_symbol="PTEN",
        consequence=ConsequenceType.STOP_GAINED, biotype="protein_coding",
        is_mane_select=mane, is_canonical=canonical,
    )


def test_fallback_flags_mane_transcript_by_refseq_base():
    cons = [_cons("NM_001304717.5", canonical=True), _cons("NM_000314.8")]
    _apply_mane_fallback(cons, _PTEN_MANE)
    flagged = [c.transcript_id for c in cons if c.is_mane_select]
    assert flagged == ["NM_000314.8"]


def test_fallback_matches_ensembl_base():
    cons = [_cons("ENST00000371953.8"), _cons("NM_001304717.5", canonical=True)]
    _apply_mane_fallback(cons, _PTEN_MANE)
    assert [c.transcript_id for c in cons if c.is_mane_select] == ["ENST00000371953.8"]


def test_fallback_noop_when_real_flag_present():
    # A genuine GRCh38 flag must never be overridden by the accession map.
    cons = [_cons("NM_001304717.5", mane=True), _cons("NM_000314.8")]
    _apply_mane_fallback(cons, _PTEN_MANE)
    assert [c.transcript_id for c in cons if c.is_mane_select] == ["NM_001304717.5"]


def test_fallback_noop_without_map():
    cons = [_cons("NM_001304717.5"), _cons("NM_000314.8")]
    _apply_mane_fallback(cons, None)
    assert not any(c.is_mane_select for c in cons)


def test_parse_record_selects_mane_first_on_grch37():
    """End-to-end: a GRCh37-style record (no mane_select) sorts MANE to front."""
    record = {
        "id": "chr10:89717708:C:T",
        "transcript_consequences": [
            {"transcript_id": "NM_001304717.5", "gene_symbol": "PTEN",
             "consequence_terms": ["stop_gained"], "canonical": 1, "biotype": "protein_coding"},
            {"transcript_id": "NM_000314.8", "gene_symbol": "PTEN",
             "consequence_terms": ["stop_gained"], "biotype": "protein_coding"},
        ],
    }
    _key, cons = _parse_vep_record(record, _PTEN_MANE)
    assert cons[0].transcript_id == "NM_000314.8"
    assert cons[0].is_mane_select
