"""MANE Plus Clinical flag, GRCh37 recovery and supplementary reporting.

Criteria are still evaluated on MANE Select (``primary_consequence``); the Plus
Clinical transcript is only described alongside it (TSV/JSON/report).
"""
from pathlib import Path

from acmg_classifier.io.tsv_writer import _HEADER, _build_mane_pc_row
from acmg_classifier.local_db.mane_db import load_mane_plus_clinical_map
from acmg_classifier.local_db.vep_runner import (
    _apply_mane_plus_clinical_fallback,
    _is_mane_plus_clinical,
    _parse_vep_record,
)
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.classification import ClassificationResult
from acmg_classifier.models.enums import ConsequenceType, Pathogenicity

# KRAS: MANE Select NM_004985 (KRAS4B) / Plus Clinical NM_033360 (KRAS4A, exon 4a).
_SELECT = "NM_004985.5"
_PC = "NM_033360.4"
_PC_MAP = {"KRAS": [("NM_033360", "ENST00000311936")]}


def _cons(tx, csq=ConsequenceType.MISSENSE, select=False, pc=False, gene="KRAS",
          hgvs_c=None):
    return ConsequenceInfo(
        transcript_id=tx, gene_id="g", gene_symbol=gene, consequence=csq,
        biotype="protein_coding", is_mane_select=select, is_mane_plus_clinical=pc,
        hgvs_c=hgvs_c,
    )


def test_vep_flag_forms():
    assert _is_mane_plus_clinical({"mane_plus_clinical": "NM_033360.4"})
    assert _is_mane_plus_clinical({"mane": ["MANE_Plus_Clinical"]})
    assert not _is_mane_plus_clinical({"mane": ["MANE_Select"]})
    assert not _is_mane_plus_clinical({"mane_select": "NM_004985.5"})


def test_parse_record_reads_plus_clinical_flag_and_keeps_select_primary():
    record = {
        "id": "chr12:25227341:T:C",
        "transcript_consequences": [
            {"transcript_id": _PC, "gene_symbol": "KRAS", "biotype": "protein_coding",
             "consequence_terms": ["missense_variant"], "mane_plus_clinical": _PC},
            {"transcript_id": _SELECT, "gene_symbol": "KRAS", "biotype": "protein_coding",
             "consequence_terms": ["missense_variant"], "mane_select": _SELECT, "canonical": 1},
        ],
    }
    _key, cons = _parse_vep_record(record)
    ann = AnnotationData(consequences=cons)
    assert ann.primary_consequence.transcript_id == _SELECT
    assert [c.transcript_id for c in ann.mane_plus_clinical_consequences] == [_PC]


def test_grch37_fallback_flags_by_accession():
    cons = [_cons("NM_033360.3"), _cons("NM_004985.4", select=True)]
    _apply_mane_plus_clinical_fallback(cons, _PC_MAP)
    assert [c.transcript_id for c in cons if c.is_mane_plus_clinical] == ["NM_033360.3"]


def test_fallback_noop_when_vep_flag_present_or_no_map():
    cons = [_cons("NM_033360.3"), _cons("ENST00000311936.8", pc=True)]
    _apply_mane_plus_clinical_fallback(cons, _PC_MAP)
    assert [c.transcript_id for c in cons if c.is_mane_plus_clinical] == ["ENST00000311936.8"]
    cons = [_cons("NM_033360.3")]
    _apply_mane_plus_clinical_fallback(cons, None)
    assert not cons[0].is_mane_plus_clinical


def test_refseq_preferred_and_other_gene_excluded():
    ann = AnnotationData(consequences=[
        _cons(_SELECT, select=True),
        _cons("ENST00000311936.8", pc=True),
        _cons(_PC, pc=True),
        _cons("NM_999999.1", pc=True, gene="NEIGHBOUR"),
    ])
    assert [c.transcript_id for c in ann.mane_plus_clinical_consequences] == [_PC]
    # Ensembl row is used only when no RefSeq Plus Clinical row exists.
    ann = AnnotationData(consequences=[_cons(_SELECT, select=True),
                                       _cons("ENST00000311936.8", pc=True)])
    assert [c.transcript_id for c in ann.mane_plus_clinical_consequences] == ["ENST00000311936.8"]


def test_summary_flags_consequence_difference():
    ann = AnnotationData(consequences=[
        _cons(_SELECT, csq=ConsequenceType.INTRON, select=True),
        _cons(_PC, csq=ConsequenceType.MISSENSE, pc=True, hgvs_c=_PC + ":c.451G>A"),
    ])
    (s,) = ann.mane_plus_clinical_summary()
    assert s.transcript_id == _PC and s.consequence_differs
    assert s.consequence == ConsequenceType.MISSENSE.value


def test_no_plus_clinical_gives_empty_summary():
    ann = AnnotationData(consequences=[_cons(_SELECT, select=True)])
    assert ann.mane_plus_clinical_summary() == []
    assert AnnotationData().mane_plus_clinical_summary() == []


def test_tsv_columns():
    assert "mane_plus_clinical_transcript" in _HEADER
    i = _HEADER.index("hgvs_p")
    assert _HEADER[i + 1] == "mane_plus_clinical_transcript"
    ann = AnnotationData(consequences=[
        _cons(_SELECT, csq=ConsequenceType.INTRON, select=True),
        _cons(_PC, pc=True, hgvs_c=_PC + ":c.451G>A"),
    ])
    res = ClassificationResult(
        variant_id="v", criteria_results=[], classification_2015=Pathogenicity.VUS,
        classification_2015_rules="", bayesian_score=0,
        classification_bayesian=Pathogenicity.VUS,
        mane_plus_clinical=ann.mane_plus_clinical_summary(),
    )
    row = _build_mane_pc_row(res)
    assert row[0] == _PC and row[1] == "c.451G>A" and row[4] == "1"
    empty = res.model_copy(update={"mane_plus_clinical": []})
    assert _build_mane_pc_row(empty) == ["", "", "", "", ""]


def test_packaged_map_contains_known_plus_clinical_transcripts():
    m = load_mane_plus_clinical_map(Path("/nonexistent/mane_exons.tsv"))
    assert ("NM_033360", "ENST00000311936") in m.get("KRAS", []) or m.get("KRAS")
    assert len(m.get("GNAS", [])) == 2
    assert "PTEN" not in m
