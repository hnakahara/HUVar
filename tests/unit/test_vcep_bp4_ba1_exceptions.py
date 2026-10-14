"""VCEP-specific BP4 scope rules and BA1/BS1 variant exceptions (Genome Medicine
revision, Reviewer 3 Comment 15 follow-up).

* LGMD VCEP "BS1 exception" variants (DYSF c.2643+1G>A on NM_003494.3 = MANE
  c.2697+1G>A, ANO5 c.191dup, CAPN3 c.1746-20C>G and c.2120A>G): neither BA1
  nor BS1 is applied.
* BP4: not applicable (SCID genes), splice-only (CDH1, PALB2), missense-only
  (HRAS, KRAS, SOS1), no synonymous (PAH), missense splice requirement (GP9 = 0,
  RUNX1 <= 0.20), PTEN intron 1 exclusion, BMPR2 all-three predictors.
"""
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.allele_frequency import variant_listed
from acmg_classifier.criteria.benign.ba1 import BA1Evaluator
from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
from acmg_classifier.criteria.benign.bs1 import BS1Evaluator
from acmg_classifier.criteria.insilico_genes import InSilicoGeneSpec, InSilicoScores
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, GnomADData, RevelData, SpliceScore,
)
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, CriterionStrength, InSilicoTool,
)
from acmg_classifier.models.variant import VariantRecord

_ROOT = Path(__file__).resolve().parents[2]
_TSV = _ROOT / "resources" / "shared" / "disease_prevalence.tsv"


def _cfg():
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _TSV
    cfg.insilico_tool = InSilicoTool.REVEL
    cfg.use_cadd = False
    cfg.use_bayesdel = False
    cfg.popmax_af_basis = True
    cfg.tp53_codes_tsv = _ROOT / "resources" / "shared" / "tp53_pp3_bp4_codes.tsv"
    return cfg


def _snv():
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


def _ann(gene, consequence, hgvs_c=None, hgvs_p=None, revel=None, splice=None, gd=None):
    return AnnotationData(
        revel=RevelData(score=revel) if revel is not None else None,
        splice=(SpliceScore(tool="spliceai", is_available=True, max_delta=splice)
                if splice is not None else None),
        gnomad=gd,
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=consequence, biotype="protein_coding",
            hgvs_c=hgvs_c, hgvs_p=hgvs_p,
        )],
    )


# --- BA1 / BS1 variant exceptions ----------------------------------------------

class TestVariantListed:
    def test_matching(self):
        assert variant_listed("c.191dup", "c.191dupA", None) == "c.191dup"
        assert variant_listed("c.1746-20C>G,c.2120A>G", "NM_000070.3:c.2120A>G", None) == "c.2120A>G"
        assert variant_listed("p.Gln368Ter", None, "p.(Gln368*)") == "p.Gln368Ter"
        assert variant_listed("c.2697+1G>A", "c.2697+2T>C", None) == ""
        assert variant_listed("", "c.1A>G", None) == ""


@pytest.mark.parametrize("gene,hgvs_c,cons", [
    ("DYSF", "c.2697+1G>A", ConsequenceType.SPLICE_DONOR),
    ("ANO5", "c.191dupA", ConsequenceType.FRAMESHIFT),
    ("CAPN3", "c.1746-20C>G", ConsequenceType.INTRON),
    ("CAPN3", "c.2120A>G", ConsequenceType.MISSENSE),
])
def test_lgmd_exceptions_block_ba1_and_bs1(gene, hgvs_c, cons):
    gd = GnomADData(faf95_popmax=0.02, popmax_af=0.02, af=0.02, filter_pass=True)
    ann = _ann(gene, cons, hgvs_c=hgvs_c, gd=gd)
    ba1 = BA1Evaluator(_cfg()).evaluate(_snv(), ann)
    bs1 = BS1Evaluator(_cfg()).evaluate(_snv(), ann)
    assert not ba1.triggered and "excludes BA1/BS1" in ba1.evidence
    assert not bs1.triggered and "excludes BS1" in bs1.evidence


def test_other_capn3_variant_still_gets_ba1():
    gd = GnomADData(faf95_popmax=0.02, popmax_af=0.02, af=0.02, filter_pass=True)
    ann = _ann("CAPN3", ConsequenceType.MISSENSE, hgvs_c="c.2121A>G", gd=gd)
    assert BA1Evaluator(_cfg()).evaluate(_snv(), ann).triggered


# --- BP4 scope ----------------------------------------------------------------

def _bp4(ann):
    return BP4Evaluator(_cfg()).evaluate(_snv(), ann)


@pytest.mark.parametrize("gene", ["DCLRE1C", "IL2RG", "IL7R", "JAK3", "RAG1"])
def test_bp4_not_applicable_scid(gene):
    r = _bp4(_ann(gene, ConsequenceType.MISSENSE, hgvs_c="c.100A>G", revel=0.01, splice=0.0))
    assert not r.triggered and "not applicable" in r.evidence
    r = _bp4(_ann(gene, ConsequenceType.SYNONYMOUS, hgvs_c="c.99A>G", splice=0.0))
    assert not r.triggered


@pytest.mark.parametrize("gene", ["CDH1", "PALB2"])
def test_bp4_splice_only(gene):
    r = _bp4(_ann(gene, ConsequenceType.MISSENSE, hgvs_c="c.100A>G", revel=0.01, splice=0.0))
    assert not r.triggered and "not applicable to missense" in r.evidence
    r = _bp4(_ann(gene, ConsequenceType.INTRON, hgvs_c="c.48+20C>T", splice=0.01))
    assert r.triggered and r.strength == CriterionStrength.SUPPORTING


@pytest.mark.parametrize("gene", ["HRAS", "KRAS", "SOS1"])
def test_bp4_missense_only(gene):
    r = _bp4(_ann(gene, ConsequenceType.MISSENSE, hgvs_c="c.100A>G", revel=0.20, splice=0.0))
    assert r.triggered
    r = _bp4(_ann(gene, ConsequenceType.SYNONYMOUS, hgvs_c="c.90C>T", splice=0.0))
    assert not r.triggered and "missense variants only" in r.evidence


def test_bp4_pah_no_synonymous():
    assert not _bp4(_ann("PAH", ConsequenceType.SYNONYMOUS, hgvs_c="c.1155C>G", splice=0.0)).triggered
    assert _bp4(_ann("PAH", ConsequenceType.INTRON, hgvs_c="c.1199+50G>A", splice=0.0)).triggered


def test_bp4_gp9_requires_spliceai_zero():
    assert _bp4(_ann("GP9", ConsequenceType.MISSENSE, hgvs_c="c.100A>G", revel=0.2, splice=0.0)).triggered
    r = _bp4(_ann("GP9", ConsequenceType.MISSENSE, hgvs_c="c.100A>G", revel=0.2, splice=0.01))
    assert not r.triggered
    r = _bp4(_ann("GP9", ConsequenceType.MISSENSE, hgvs_c="c.100A>G", revel=0.2))
    assert not r.triggered and "no splice prediction" in r.evidence
    assert not _bp4(_ann("GP9", ConsequenceType.SYNONYMOUS, hgvs_c="c.99A>G", splice=0.01)).triggered
    assert _bp4(_ann("GP9", ConsequenceType.SYNONYMOUS, hgvs_c="c.99A>G", splice=0.0)).triggered


def test_bp4_runx1_missense_splice_020():
    assert _bp4(_ann("RUNX1", ConsequenceType.MISSENSE, hgvs_c="c.100A>G", revel=0.3, splice=0.20)).triggered
    assert not _bp4(_ann("RUNX1", ConsequenceType.MISSENSE, hgvs_c="c.100A>G", revel=0.3, splice=0.21)).triggered


@pytest.mark.parametrize("hgvs_c,blocked", [
    ("c.79+7A>G", True), ("c.80-3C>T", True), ("c.78C>T", True),
    ("c.42G>A", False), ("c.790C>T", False), ("c.165-20T>C", False),
])
def test_bp4_pten_intron1(hgvs_c, blocked):
    cons = ConsequenceType.INTRON if ("+" in hgvs_c or "-" in hgvs_c) else ConsequenceType.SYNONYMOUS
    r = _bp4(_ann("PTEN", cons, hgvs_c=hgvs_c, splice=0.0))
    assert r.triggered is (not blocked)


def test_bmpr2_bp4_requires_all_three():
    spec = InSilicoGeneSpec()
    sc = lambda **k: InSilicoScores(consequence=ConsequenceType.MISSENSE, protein_position=None,
                                    revel=k.get("revel"), alphamissense=k.get("am"), cadd=k.get("cadd"),
                                    bayesdel=None, splice_max=0.0)
    assert spec.bp4("BMPR2", sc(cadd=15.6, am=0.07, revel=0.27))[0] == CriterionStrength.SUPPORTING
    assert spec.bp4("BMPR2", sc(cadd=22.1, am=0.087, revel=0.788))[0] is None
    assert spec.bp4("BMPR2", sc(cadd=21.0, am=0.11))[0] is None


@pytest.mark.parametrize("gene,hgvs_c,blocked", [
    ("GAA", "c.-32-13T>G", True), ("RUNX1", "c.-59-9G>A", True),
    ("SLC9A6", "c.-57+27T>G", True), ("NOSPEC", "NM_1.1:c.-20+5G>A", True),
    ("NOSPEC", "c.100+30G>A", False), ("NOSPEC", "c.-30G>A", False),
])
def test_bp4_not_applied_to_5utr_introns(gene, hgvs_c, blocked):
    r = _bp4(_ann(gene, ConsequenceType.INTRON, hgvs_c=hgvs_c, splice=0.0))
    assert r.triggered is (not blocked)
    if blocked:
        assert "5'UTR intron" in r.evidence
