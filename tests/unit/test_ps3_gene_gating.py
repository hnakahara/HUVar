"""PS3 per-gene gating: suppression (no PS3 / in-vitro N/A) and Supporting cap.

The text-mined PS3 fallback is gene-agnostic by count, so per-gene VCEP limits
are enforced in the evaluator: PALB2/PDHA1/POLG (no PS3 code) and CAPN3/ANO5
(PS3 not applicable for in-vitro assays) are suppressed; RPE65/SERPINC1/... cap
PS3 at Supporting so a >=3-SCV count cannot reach Moderate.
"""
from unittest.mock import MagicMock

import acmg_classifier.local_db.clinvar_sqlite as clinvar_sqlite
from acmg_classifier.criteria.pathogenic.ps3 import PS3Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.enums import (
    ACMGCriterion, Assembly, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.models.variant import VariantRecord


def _cfg(tmp_path):
    cfg = MagicMock()
    cfg.clinvar_sqlite = tmp_path / "clinvar.sqlite"
    return cfg


def _ann(gene):
    return AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
    )])


def _var():
    return VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)


def _patch_count(monkeypatch, n):
    monkeypatch.setattr(clinvar_sqlite, "query_functional_evidence", lambda *a, **k: n)


class TestSuppression:
    def test_no_ps3_gene_suppressed(self, tmp_path, monkeypatch):
        _patch_count(monkeypatch, 3)
        r = PS3Evaluator(_cfg(tmp_path)).evaluate(_var(), _ann("PALB2"))
        assert not r.triggered and "does not permit" in r.evidence

    def test_in_vitro_na_gene_suppressed(self, tmp_path, monkeypatch):
        _patch_count(monkeypatch, 5)
        r = PS3Evaluator(_cfg(tmp_path)).evaluate(_var(), _ann("CAPN3"))
        assert not r.triggered

    def test_supplement_overrides_suppression(self, tmp_path, monkeypatch):
        # A curated animal-model PS3 still applies for a suppressed gene.
        _patch_count(monkeypatch, 0)
        sup = [SupplementEntry(
            variant_id="chr1:100:C:T", criterion=ACMGCriterion.PS3,
            strength=CriterionStrength.MODERATE, evidence="variant-specific mouse model",
        )]
        r = PS3Evaluator(_cfg(tmp_path)).evaluate(_var(), _ann("CAPN3"), sup)
        assert r.triggered and r.strength == CriterionStrength.MODERATE


class TestSupportingCap:
    def test_rpe65_capped_at_supporting(self, tmp_path, monkeypatch):
        _patch_count(monkeypatch, 3)   # would be Moderate without the cap
        r = PS3Evaluator(_cfg(tmp_path)).evaluate(_var(), _ann("RPE65"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_serpinc1_capped(self, tmp_path, monkeypatch):
        _patch_count(monkeypatch, 4)
        r = PS3Evaluator(_cfg(tmp_path)).evaluate(_var(), _ann("SERPINC1"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING


class TestUncappedGene:
    def test_normal_gene_reaches_moderate(self, tmp_path, monkeypatch):
        _patch_count(monkeypatch, 3)
        r = PS3Evaluator(_cfg(tmp_path)).evaluate(_var(), _ann("GAA"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_normal_gene_single_supporting(self, tmp_path, monkeypatch):
        _patch_count(monkeypatch, 1)
        r = PS3Evaluator(_cfg(tmp_path)).evaluate(_var(), _ann("GAA"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING
