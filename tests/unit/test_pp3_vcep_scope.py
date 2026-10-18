"""VCEP PP3 scope (pp3_scope) and the extended BP4 scope (Genome Medicine
revision, PP3 discordance audit):

* SCID VCEP (ADA, DCLRE1C, IL2RG, IL7R, JAK3, RAG1, RAG2): PP3 only for
  synonymous/intronic splice impact; "Do not apply to missense variants".
* APC, PALB2, CDH1: no protein-level predictor; a missense meets PP3 only via
  predicted splice impact.
* PIK3/mTOR GOF panel (AKT3, MTOR, PIK3CA, PIK3R2): PP3 not applicable.
* BP4: ADA/RAG2 not applicable; APC and the PIK3/mTOR panel splice-only.
"""
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, RevelData, SpliceScore,
)
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, CriterionStrength, InSilicoTool,
)
from acmg_classifier.models.variant import VariantRecord

_ROOT = Path(__file__).resolve().parents[2]


def _cfg():
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _ROOT / "resources" / "shared" / "disease_prevalence.tsv"
    cfg.insilico_tool = InSilicoTool.REVEL
    cfg.use_cadd = False
    cfg.use_bayesdel = False
    cfg.tp53_codes_tsv = _ROOT / "resources" / "shared" / "tp53_pp3_bp4_codes.tsv"
    return cfg


def _snv():
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


def _ann(gene, consequence, revel=None, splice=None):
    return AnnotationData(
        revel=RevelData(score=revel) if revel is not None else None,
        splice=(SpliceScore(tool="spliceai", is_available=True, max_delta=splice)
                if splice is not None else None),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=consequence, biotype="protein_coding",
            hgvs_c="c.100A>G", hgvs_p="p.Lys34Glu",
        )],
    )


def _pp3(gene, cons, **kw):
    return PP3Evaluator(_cfg()).evaluate(_snv(), _ann(gene, cons, **kw))


def _bp4(gene, cons, **kw):
    return BP4Evaluator(_cfg()).evaluate(_snv(), _ann(gene, cons, **kw))


M = ConsequenceType.MISSENSE


@pytest.mark.parametrize("gene", ["ADA", "DCLRE1C", "IL2RG", "IL7R", "JAK3", "RAG1", "RAG2"])
def test_scid_missense_never_pp3(gene):
    r = _pp3(gene, M, revel=0.99, splice=0.9)
    assert not r.triggered
    assert "missense" in r.evidence


def test_scid_intronic_splice_still_pp3():
    r = _pp3("RAG2", ConsequenceType.INTRON, splice=0.5)
    assert r.triggered and r.strength == CriterionStrength.SUPPORTING


@pytest.mark.parametrize("gene", ["APC", "PALB2", "CDH1"])
def test_splice_only_missense(gene):
    assert not _pp3(gene, M, revel=0.99, splice=0.0).triggered
    r = _pp3(gene, M, revel=0.99, splice=0.6)
    assert r.triggered and "max_delta" in r.evidence


@pytest.mark.parametrize("gene", ["AKT3", "MTOR", "PIK3CA", "PIK3R2"])
def test_gof_panel_pp3_not_applicable(gene):
    assert not _pp3(gene, M, revel=0.99).triggered
    assert not _pp3(gene, ConsequenceType.INTRON, splice=0.9).triggered


def test_unrestricted_gene_unchanged():
    r = _pp3("MECP2", M, revel=0.95)
    assert r.triggered and "REVEL" in r.evidence


@pytest.mark.parametrize("gene", ["ADA", "RAG2"])
def test_bp4_not_applicable_scid(gene):
    assert not _bp4(gene, M, revel=0.01).triggered
    assert not _bp4(gene, ConsequenceType.INTRON, splice=0.0).triggered


@pytest.mark.parametrize("gene", ["APC", "MTOR", "PIK3CA"])
def test_bp4_splice_only_missense(gene):
    assert not _bp4(gene, M, revel=0.01, splice=0.0).triggered
