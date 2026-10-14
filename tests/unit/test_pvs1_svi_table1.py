"""ClinGen SVI PVS1 recommendation (Abou Tayoun 2018): PM4 never with PVS1, and
PVS1 two levels down for genes without a Definitive/Strong gene-disease
relationship (no VCEP)."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion as C, CriterionStrength as S, ConsequenceType


def _gdv(tmp_path):
    p = tmp_path / "gdv.csv"
    p.write_text('"GENE SYMBOL","CLASSIFICATION"\n"DEFG","Definitive"\n"STRG","Strong"\n'
                 '"MODG","Moderate"\n"MODG","Limited"\n"LIMG","Limited"\n"VCEPG","Moderate"\n',
                 encoding="utf-8")
    dp = tmp_path / "dp.tsv"
    dp.write_text("gene_symbol\tinheritance\tsource_vcep\nVCEPG\tAD\tX VCEP\n", encoding="utf-8")
    return p, dp


def _ev(tmp_path, on=True):
    from acmg_classifier.classification.gene_validity import _vcep_genes
    from acmg_classifier.local_db.clingen_gene_db import load_gene_validity
    from acmg_classifier.criteria.pathogenic.pvs1 import PVS1Evaluator
    _vcep_genes.cache_clear(); load_gene_validity.cache_clear()
    p, dp = _gdv(tmp_path)
    ev = PVS1Evaluator.__new__(PVS1Evaluator)
    ev._cfg = SimpleNamespace(clingen_gene_validity_csv=p, disease_prevalence_tsv=dp,
                              pvs1_gdv_downgrade=on)
    return ev


@pytest.mark.parametrize("gene,start,expect", [
    ("DEFG", S.VERY_STRONG, S.VERY_STRONG),
    ("STRG", S.VERY_STRONG, S.VERY_STRONG),
    ("MODG", S.VERY_STRONG, S.MODERATE),
    ("NOTCURATED", S.VERY_STRONG, S.MODERATE),
    ("LIMG", S.STRONG, S.SUPPORTING),
    ("MODG", S.MODERATE, S.NOT_MET),
    ("VCEPG", S.VERY_STRONG, S.VERY_STRONG),   # VCEP gene: follows the VCEP
])
def test_gdv_downgrade(tmp_path, gene, start, expect):
    st, ev = _ev(tmp_path)._gene_validity_downgrade(gene, start, "x")
    assert st == expect
    if st != start:
        assert "decreased two levels" in ev


def test_gdv_downgrade_disabled(tmp_path):
    assert _ev(tmp_path, on=False)._gene_validity_downgrade("MODG", S.VERY_STRONG, "x")[0] == S.VERY_STRONG


def test_pm4_not_with_pvs1():
    from acmg_classifier.criteria.registry import CriteriaRegistry
    reg = CriteriaRegistry.__new__(CriteriaRegistry)
    reg._pm4_regions = SimpleNamespace(excludes=lambda g: [])
    ann = AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_x", gene_id="E", gene_symbol="ANY",
        consequence=ConsequenceType.INFRAME_DELETION, biotype="protein_coding")])
    res = [CriteriaResult.met(C.PVS1, S.MODERATE, "x"), CriteriaResult.met(C.PM4, S.MODERATE, "y")]
    reg._apply_pm4_exclusions(res, ann)
    assert res[1].suppressed and "PM4 not with PVS1" in res[1].evidence
    res = [CriteriaResult.not_met(C.PVS1, "x"), CriteriaResult.met(C.PM4, S.MODERATE, "y")]
    reg._apply_pm4_exclusions(res, ann)
    assert not res[1].suppressed
