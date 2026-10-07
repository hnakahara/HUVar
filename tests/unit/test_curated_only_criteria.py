"""PS4 / PP1 are curated-evidence-only (GM revision) + expert-panel mining."""
from unittest.mock import MagicMock

import acmg_classifier.local_db.clinvar_sqlite as clinvar_sqlite
from acmg_classifier.criteria.pathogenic.pp1 import PP1Evaluator
from acmg_classifier.criteria.pathogenic.ps4 import PS4Evaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.enums import ACMGCriterion, Assembly, CriterionStrength
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.setup.clinvar_builder import _mine_ep_criteria


def _cfg(tmp_path):
    cfg = MagicMock()
    cfg.clinvar_sqlite = tmp_path / "none.sqlite"
    cfg.exclude_self_expert_panel = False
    return cfg


V = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)


def test_ps4_not_automated_from_counts(tmp_path, monkeypatch):
    monkeypatch.setattr(clinvar_sqlite, "query_affected_cases", lambda *a, **k: 25)
    r = PS4Evaluator(_cfg(tmp_path)).evaluate(V, AnnotationData())
    assert not r.triggered and "25 P/LP ClinVar submission" in r.evidence


def test_ps4_from_supplement(tmp_path):
    sup = [SupplementEntry(variant_id=V.key, criterion=ACMGCriterion.PS4,
                           strength=CriterionStrength.MODERATE, evidence="eRepo")]
    r = PS4Evaluator(_cfg(tmp_path)).evaluate(V, AnnotationData(), sup)
    assert r.triggered and r.strength == CriterionStrength.MODERATE


def test_pp1_from_expert_panel(tmp_path, monkeypatch):
    monkeypatch.setattr(clinvar_sqlite, "query_expert_panel_criteria",
                        lambda *a, **k: {"PP1": "Strong"})
    r = PP1Evaluator(_cfg(tmp_path)).evaluate(V, AnnotationData())
    assert r.triggered and r.strength == CriterionStrength.STRONG


def test_pp1_text_mining_info_only(tmp_path, monkeypatch):
    monkeypatch.setattr(clinvar_sqlite, "query_segregation_evidence", lambda *a, **k: 3)
    r = PP1Evaluator(_cfg(tmp_path)).evaluate(V, AnnotationData())
    assert not r.triggered and "mention co-segregation" in r.evidence


def test_mine_ep_criteria():
    t = ("The c.100C>T variant ... criteria applied: PS3_Moderate, PM2_Supporting, "
         "PP1_Strong, PS4.")
    assert _mine_ep_criteria(t) == "PP1:Strong;PS3:Moderate;PS4:Strong"
    assert _mine_ep_criteria("PS3 not met. PM2_Supporting") is None
    assert _mine_ep_criteria("Criteria not met: PP1. PS3_Supporting applied") == "PS3:Supporting"
    assert _mine_ep_criteria("") is None
