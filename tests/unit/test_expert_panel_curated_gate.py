"""Expert-panel (>=3-star) ClinVar criteria (PS3/PS4/PP1/BS2) are curated
evidence: they are imported only when curated evidence is in use for the run
(``cfg.curated_evidence``, set by the pipeline when a supplement is supplied)."""
from unittest.mock import MagicMock, patch

import pytest

import acmg_classifier.local_db.clinvar_sqlite as clinvar_sqlite
from acmg_classifier.config import Config
from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
from acmg_classifier.criteria.curated_sources import expert_panel_enabled
from acmg_classifier.criteria.pathogenic.pp1 import PP1Evaluator
from acmg_classifier.criteria.pathogenic.ps3 import PS3Evaluator
from acmg_classifier.criteria.pathogenic.ps4 import PS4Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import ACMGCriterion, Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.models.variant import VariantRecord

V = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)


def _cfg(tmp_path, curated):
    dp = tmp_path / "disease_prevalence.tsv"
    dp.write_text("gene_symbol\tbs2\n", encoding="utf-8")
    cfg = MagicMock()
    cfg.clinvar_sqlite = tmp_path / "none.sqlite"
    cfg.disease_prevalence_tsv = dp
    cfg.exclude_self_expert_panel = False
    cfg.curated_evidence = curated
    cfg.use_tp53_functional = False
    return cfg


def test_default_config_is_automated_only():
    assert Config().curated_evidence is False
    assert not expert_panel_enabled(Config())
    assert expert_panel_enabled(Config(curated_evidence=True))
    assert not expert_panel_enabled(Config(curated_evidence=True, exclude_self_expert_panel=True))


@pytest.mark.parametrize("ev_cls,crit", [(PP1Evaluator, "PP1"), (PS4Evaluator, "PS4")])
def test_pp1_ps4_expert_panel_only_with_curated_evidence(tmp_path, monkeypatch, ev_cls, crit):
    monkeypatch.setattr(clinvar_sqlite, "query_expert_panel_criteria",
                        lambda *a, **k: {crit: "Strong"})
    assert not ev_cls(_cfg(tmp_path, False)).evaluate(V, AnnotationData()).triggered
    r = ev_cls(_cfg(tmp_path, True)).evaluate(V, AnnotationData())
    assert r.triggered and "[ClinVar expert panel]" in r.evidence


def test_ps3_expert_panel_only_with_curated_evidence(tmp_path, monkeypatch):
    monkeypatch.setattr(clinvar_sqlite, "query_expert_panel_criteria",
                        lambda *a, **k: {"PS3": "Moderate"})
    monkeypatch.setattr(clinvar_sqlite, "query_functional_pmids", lambda *a, **k: [], raising=False)
    monkeypatch.setattr(clinvar_sqlite, "query_functional_evidence", lambda *a, **k: 0, raising=False)
    ann = AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol="GAA",
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding")])
    r = PS3Evaluator(_cfg(tmp_path, False)).evaluate(V, ann)
    assert "[ClinVar expert panel]" not in (r.evidence or "")
    r = PS3Evaluator(_cfg(tmp_path, True)).evaluate(V, ann)
    assert r.triggered and r.strength == CriterionStrength.MODERATE
    assert "[ClinVar expert panel]" in r.evidence


def test_bs2_expert_panel_only_with_curated_evidence(tmp_path):
    ann = AnnotationData(
        gnomad=GnomADData(filter_pass=True, ac=0),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol="NOVCEP",
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding")])
    with patch("acmg_classifier.local_db.clinvar_sqlite.query_bs2_benign_evidence",
               return_value=(True, "Strong")):
        r = BS2Evaluator(_cfg(tmp_path, False)).evaluate(V, ann)
        assert not r.triggered and "no curated evidence supplied" in r.evidence
        assert BS2Evaluator(_cfg(tmp_path, True)).evaluate(V, ann).triggered


def test_supplement_still_applies_without_flag(tmp_path):
    sup = [SupplementEntry(variant_id=V.key, criterion=ACMGCriterion.PP1,
                           strength=CriterionStrength.MODERATE, evidence="curator")]
    r = PP1Evaluator(_cfg(tmp_path, False)).evaluate(V, AnnotationData(), sup)
    assert r.triggered and r.strength == CriterionStrength.MODERATE


def test_run_single_sets_flag_only_with_manual_evidence(monkeypatch):
    import acmg_classifier.pipeline.pipeline as pl
    seen = []
    monkeypatch.setattr(pl, "annotate_one", lambda *a, **k: (V, AnnotationData()))
    monkeypatch.setattr(pl, "classify_annotated",
                        lambda v, a, cfg, s=None: seen.append(cfg.curated_evidence) or MagicMock())
    monkeypatch.setattr("acmg_classifier.io.report_writer.print_report", lambda *a, **k: None)
    pl.run_single("chr1", 100, "C", "T", Config())
    sup = [SupplementEntry(variant_id=V.key, criterion=ACMGCriterion.PS3,
                           strength=CriterionStrength.STRONG, evidence="x")]
    pl.run_single("chr1", 100, "C", "T", Config(), supplement=sup)
    assert seen == [False, True]
