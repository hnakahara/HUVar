"""Unit tests for ESM1b UniProt resolution (RefSeq→Ensembl fallback).

VEP --merged attaches swissprot/trembl xrefs only to Ensembl (ENST)
transcripts. Since primary_consequence prefers the RefSeq MANE transcript,
the primary often has no uniprot_id even for scorable missense variants
(e.g. TP53 R248W). _resolve_uniprot_id must fall back to a sibling Ensembl
transcript describing the same substitution.
"""
from acmg_classifier.annotation.orchestrator import _resolve_uniprot_id
from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType


def _refseq(uniprot=None, ctype=ConsequenceType.MISSENSE,
            protein_pos=248, amino_acids="R/W"):
    return ConsequenceInfo(
        transcript_id="NM_000546.6",
        gene_id="ENSG00000141510",
        gene_symbol="TP53",
        consequence=ctype,
        biotype="protein_coding",
        is_mane_select=True,
        protein_position=protein_pos,
        amino_acids=amino_acids,
        uniprot_id=uniprot,
    )


def _ensembl(uniprot="P04637", ctype=ConsequenceType.MISSENSE,
             protein_pos=248, amino_acids="R/W"):
    return ConsequenceInfo(
        transcript_id="ENST00000269305",
        gene_id="ENSG00000141510",
        gene_symbol="TP53",
        consequence=ctype,
        biotype="protein_coding",
        is_mane_select=True,
        protein_position=protein_pos,
        amino_acids=amino_acids,
        uniprot_id=uniprot,
    )


class TestResolveUniprotId:
    def test_primary_has_uniprot_returns_it_directly(self):
        primary = _refseq(uniprot="P04637")
        assert _resolve_uniprot_id(primary, [primary]) == "P04637"

    def test_falls_back_to_sibling_ensembl(self):
        """The TP53 R248W case: RefSeq primary has no uniprot, Ensembl does."""
        primary = _refseq(uniprot=None)
        sibling = _ensembl(uniprot="P04637")
        assert _resolve_uniprot_id(primary, [primary, sibling]) == "P04637"

    def test_no_sibling_with_uniprot_returns_none(self):
        primary = _refseq(uniprot=None)
        sibling = _ensembl(uniprot=None)
        assert _resolve_uniprot_id(primary, [primary, sibling]) is None

    def test_sibling_different_position_is_ignored(self):
        primary = _refseq(uniprot=None, protein_pos=248)
        sibling = _ensembl(uniprot="P04637", protein_pos=249)
        assert _resolve_uniprot_id(primary, [primary, sibling]) is None

    def test_sibling_different_substitution_is_ignored(self):
        primary = _refseq(uniprot=None, amino_acids="R/W")
        sibling = _ensembl(uniprot="P04637", amino_acids="R/Q")
        assert _resolve_uniprot_id(primary, [primary, sibling]) is None

    def test_non_missense_sibling_is_ignored(self):
        primary = _refseq(uniprot=None)
        sibling = _ensembl(uniprot="P04637", ctype=ConsequenceType.SYNONYMOUS)
        assert _resolve_uniprot_id(primary, [primary, sibling]) is None

    def test_none_primary_returns_none(self):
        assert _resolve_uniprot_id(None, []) is None

    def test_empty_consequences_with_unscorable_primary(self):
        primary = _refseq(uniprot=None)
        assert _resolve_uniprot_id(primary, []) is None

    def test_first_matching_sibling_wins(self):
        primary = _refseq(uniprot=None)
        first = _ensembl(uniprot="P04637")
        second = _ensembl(uniprot="Q99999")
        assert _resolve_uniprot_id(primary, [primary, first, second]) == "P04637"
