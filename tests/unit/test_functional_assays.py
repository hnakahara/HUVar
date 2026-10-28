"""PS3/BS3 from calibrated functional assays: MaveDB clinical calibrations (Brnich
2019 OddsPath strengths, per-gene VCEP rules) and the ENIGMA BRCA1/2 Table 9.
Uses the packaged resources (MaveDB dump 2026-06-24, Table 9 v1.2)."""
from pathlib import Path

import pytest

from acmg_classifier.config import Config
from acmg_classifier.criteria.benign.manual import ManualBenignEvaluator
from acmg_classifier.criteria.functional_assays import (
    SRC_ENIGMA_T9, SRC_MAVEDB, FunctionalAssays, oddspath_strength,
)
from acmg_classifier.criteria.pathogenic.ps3 import PS3Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, SpliceScore
from acmg_classifier.models.enums import ACMGCriterion, Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

S = CriterionStrength
CT = ConsequenceType


def _cfg(**kw):
    return Config(data_dir=Path("./nonexistent-data"), **kw)


def _ann(gene, tx, hgvs_p=None, csq=CT.MISSENSE, hgvs_c=None, splice=0.01):
    c = ConsequenceInfo(transcript_id=tx, gene_id="g", gene_symbol=gene, consequence=csq,
                        biotype="protein_coding", is_mane_select=True,
                        hgvs_p=hgvs_p, hgvs_c=hgvs_c)
    sp = None if splice is None else SpliceScore(tool="openspliceai", max_delta=splice, is_available=True)
    return AnnotationData(consequences=[c], splice=sp)


def _v(chrom="chr2", pos=1, ref="A", alt="G"):
    return VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=Assembly.GRCH38)


@pytest.mark.parametrize("crit,odds,exp", [
    ("PS3", 2.0, None), ("PS3", 2.2, S.SUPPORTING), ("PS3", 9.09, S.MODERATE),
    ("PS3", 24.9, S.STRONG), ("PS3", 400, S.VERY_STRONG),
    ("BS3", 0.5, None), ("BS3", 0.47, S.SUPPORTING), ("BS3", 0.101, S.MODERATE), ("BS3", 0.043, S.STRONG),
])
def test_brnich_table3(crit, odds, exp):
    assert oddspath_strength(crit, odds) == exp


def test_msh2_primary_vs_conservative():
    # MSH2 p.Thr564Pro: primary (Scott, OddsPath 24.9) Strong; IGVF all-variants 9.09 Moderate.
    ann = _ann("MSH2", "NM_000251.3", "NP_000242.1:p.Thr564Pro")
    fa = FunctionalAssays(_cfg()).call(_v(), ann)
    assert fa.criterion == ACMGCriterion.PS3 and fa.strength == S.STRONG
    assert fa.evidence.startswith(SRC_MAVEDB) and "urn:mavedb:00000050-a-1" in fa.evidence
    fa = FunctionalAssays(_cfg(mavedb_calibration="conservative")).call(_v(), ann)
    assert fa.strength == S.MODERATE


def test_strength_recomputed_from_oddspath_not_label():
    # TP53 p.Ala138Gly primary IGVF calibration: OddsPath 10.03 -> Moderate.
    fa = FunctionalAssays(_cfg()).call(_v("chr17"), _ann("TP53", "NM_000546.6", "NP_000537.3:p.Ala138Gly"))
    assert fa.criterion == ACMGCriterion.PS3 and fa.strength == S.MODERATE


def test_bs3_needs_low_splice_score_for_protein_level_assay():
    ann = _ann("MSH2", "NM_000251.3", "NP_000242.1:p.Tyr238Cys", splice=0.05)
    fa = FunctionalAssays(_cfg()).call(_v(), ann)
    assert fa.criterion == ACMGCriterion.BS3 and fa.strength == S.STRONG
    for sp in (0.35, None):
        fa = FunctionalAssays(_cfg()).call(_v(), _ann("MSH2", "NM_000251.3", "NP_000242.1:p.Tyr238Cys", splice=sp))
        assert fa.criterion is None and "BS3 withheld" in fa.evidence


def test_vcep_assay_list_gene_not_used():
    # GCK: VCEP decision tree needs enzyme kinetics -> MaveDB not used.
    assert FunctionalAssays(_cfg()).call(_v(), _ann("GCK", "NM_000162.5", "NP_000153.1:p.Phe269Gln")) is None


def test_null_consequence_skipped_and_no_data_is_none():
    fa = FunctionalAssays(_cfg())
    assert fa.call(_v(), _ann("MSH2", "NM_000251.3", "NP_000242.1:p.Thr564Ter", csq=CT.STOP_GAINED)) is None
    assert fa.call(_v(), _ann("MSH2", "NM_000251.3", "NP_000242.1:p.Thr9999Pro")) is None
    assert FunctionalAssays(_cfg(functional_assays=False)).call(
        _v(), _ann("MSH2", "NM_000251.3", "NP_000242.1:p.Thr564Pro")) is None


def test_genomic_level_sge_key():
    # BRCA2 is Table 9-governed; the genomic matching is checked on the raw lookup.
    fa = FunctionalAssays(_cfg())
    keys = fa._keys(_v("chr13", 32356514, "G", "C"), _ann("BRCA2", "NM_000059.4", None), _ann(
        "BRCA2", "NM_000059.4").primary_consequence)
    assert ("genomic", "chr13:32356514:G:C") in keys


def test_enigma_table9_codes():
    fa = FunctionalAssays(_cfg())
    ps3 = fa.call(_v("chr17"), _ann("BRCA1", "NM_007294.4", hgvs_c="NM_007294.4:c.100C>T"))
    assert ps3.criterion == ACMGCriterion.PS3 and ps3.strength == S.STRONG
    assert ps3.evidence.startswith(SRC_ENIGMA_T9) and "Findlay 2018" in ps3.evidence
    bs3 = fa.call(_v("chr17"), _ann("BRCA1", "NM_007294.4", hgvs_c="NM_007294.4:c.100C>G"))
    assert bs3.criterion == ACMGCriterion.BS3 and bs3.strength == S.STRONG
    none = fa.call(_v("chr17"), _ann("BRCA1", "NM_007294.4", hgvs_c="NM_007294.4:c.100C>A"))
    assert none.criterion is None and "Table 9: None" in none.evidence
    # Table 9 also codes canonical-splice variants (SGE measures mRNA): applied as listed.
    sp = fa.call(_v("chr17"), _ann("BRCA1", "NM_007294.4", csq=CT.SPLICE_DONOR,
                                   hgvs_c="NM_007294.4:c.301+1G>A"))
    assert sp.criterion == ACMGCriterion.BS3
    # Not in Table 9 -> no functional code (MaveDB not consulted for BRCA1/2).
    assert fa.call(_v("chr17"), _ann("BRCA1", "NM_007294.4", "NP_009225.1:p.Pro34Ser",
                                     hgvs_c="NM_007294.4:c.99999A>G")) is None


def test_ps3_and_bs3_evaluators_use_the_call():
    cfg = _cfg()
    v = _v()
    ps3 = PS3Evaluator(cfg).evaluate(v, _ann("MSH2", "NM_000251.3", "NP_000242.1:p.Thr564Pro"))
    assert ps3.triggered and ps3.strength == S.STRONG
    # A normal-function variant: PS3 not met with the MaveDB evidence (no text mining).
    ann = _ann("MSH2", "NM_000251.3", "NP_000242.1:p.Tyr238Cys", splice=0.05)
    ps3 = PS3Evaluator(cfg).evaluate(v, ann)
    assert not ps3.triggered and ps3.evidence.startswith(SRC_MAVEDB)
    bs3 = [r for r in ManualBenignEvaluator(cfg).evaluate(v, ann) if r.criterion == ACMGCriterion.BS3][0]
    assert bs3.triggered and bs3.strength == S.STRONG


def test_tp53_uses_mavedb_when_nci_data_off():
    ann = _ann("TP53", "NM_000546.6", "NP_000537.3:p.Ala129Ser", splice=0.01)
    bs3 = [r for r in ManualBenignEvaluator(_cfg()).evaluate(_v("chr17"), ann)
           if r.criterion == ACMGCriterion.BS3][0]
    assert bs3.triggered and bs3.evidence.startswith(SRC_MAVEDB)
