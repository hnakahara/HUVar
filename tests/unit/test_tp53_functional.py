"""TP53 VCEP functional rules (PS3/BS3) from systematic assay data
(Kato / Giacomelli / Kotler / Kawaguchi; resources/shared/tp53_functional.tsv)."""
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.tp53_functional import (
    SRC_TP53_FUNCTIONAL, TP53Assays, TP53Functional,
)
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion as C, Assembly, ConsequenceType, CriterionStrength as S,
)
from acmg_classifier.models.variant import VariantRecord

_TSV = ("# test\naa_change\tkato\tgiacomelli\tkotler\tkawaguchi\n"
        "R175H\tnonfunctional\tLOF\tLOF\t\n"
        "H214L\tnonfunctional\tLOF\tnoLOF\t\n"
        "A347D\tnonfunctional\tnoLOF\t\tLOF\n"
        "X1Y\tnonfunctional\tnoLOF\tnoLOF\tLOF\n"
        "R175C\tpartial\tnoLOF\tnoLOF\t\n"
        "R337H\tpartial\tnoLOF\t\tLOF\n"
        "P3Q\tpartial\tLOF\tLOF\t\n"
        "P295S\tfunctional\t\t\t\n"
        "P72R\tfunctional\tnoLOF\t\t\n"
        "M1I\t\tLOF\t\t\n"
        "M1K\t\tnoLOF\t\t\n")


@pytest.fixture
def tp53(tmp_path):
    p = tmp_path / "tp53_functional.tsv"
    p.write_text(_TSV, encoding="utf-8")
    return TP53Functional(p)


@pytest.mark.parametrize("aa,expect", [
    ("R175H", (C.PS3, S.STRONG)),
    ("H214L", (C.PS3, S.STRONG)),       # Kato counted in the majority: 2 of 3 LOF
    ("A347D", (C.PS3, S.STRONG)),       # Kato + Kawaguchi LOF vs Giacomelli
    ("X1Y", (C.PS3, S.SUPPORTING)),     # tie -> abnormal oligomerisation branch
    ("R175C", (C.BS3, S.SUPPORTING)),
    ("R337H", None),                    # partial, not all other assays no-LOF
    ("P3Q", (C.PS3, S.MODERATE)),
    ("P295S", (C.BS3, S.STRONG)),       # Kato functional alone
    ("P72R", (C.BS3, S.STRONG)),
    ("M1I", (C.PS3, S.SUPPORTING)),     # no Kato data
    ("M1K", (C.BS3, S.SUPPORTING)),
])
def test_flowchart(tp53, aa, expect):
    d = tp53.decide(tp53.lookup(aa))
    assert (d[:2] if d else None) == expect


def test_aa_change_parsing():
    assert TP53Functional.aa_change("R/H", 175) == "R175H"
    assert TP53Functional.aa_change(None, None, "ENSP1:p.Arg175His") == "R175H"
    assert TP53Functional.aa_change("R/*", 175) is None
    assert TP53Functional.aa_change(None, None, float("nan")) is None


def test_packaged_table_loads():
    t = TP53Functional(Path("does/not/exist/tp53_functional.tsv"))  # falls back to resources
    assert t and t.decide(t.lookup("R175H"))[:2] == (C.PS3, S.STRONG)


def _ann(gene="TP53", aa="R/H", pos=175):
    return AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_000546", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        amino_acids=aa, protein_position=pos)])


def _cfg(tmp_path):
    p = tmp_path / "tp53_functional.tsv"
    p.write_text(_TSV, encoding="utf-8")
    cfg = MagicMock()
    cfg.tp53_functional_tsv = p
    cfg.use_tp53_functional = True
    cfg.exclude_self_expert_panel = True
    cfg.clinvar_sqlite = tmp_path / "none.sqlite"
    return cfg


_V = VariantRecord(chrom="chr17", pos=7675088, ref="C", alt="T", assembly=Assembly.GRCH38)


def test_ps3_evaluator_uses_functional_data(tmp_path):
    from acmg_classifier.criteria.pathogenic.ps3 import PS3Evaluator
    r = PS3Evaluator(_cfg(tmp_path)).evaluate(_V, _ann())
    assert r.triggered and r.strength == S.STRONG and r.evidence.startswith(SRC_TP53_FUNCTIONAL)
    # data present but flowchart gives BS3 -> PS3 not met, no text-mining fallback
    r = PS3Evaluator(_cfg(tmp_path)).evaluate(_V, _ann(aa="P/S", pos=295))
    assert not r.triggered and SRC_TP53_FUNCTIONAL in r.evidence


def test_ps3_supplement_takes_priority(tmp_path):
    from acmg_classifier.criteria.pathogenic.ps3 import PS3Evaluator
    from acmg_classifier.models.supplement import SupplementEntry
    sup = [SupplementEntry(variant_id="chr17:7675088:C:T", criterion=C.PS3,
                           strength=S.SUPPORTING, evidence="curator")]
    r = PS3Evaluator(_cfg(tmp_path)).evaluate(_V, _ann(), sup)
    assert r.strength == S.SUPPORTING and "curator" in r.evidence


def test_bs3_from_manual_benign_evaluator(tmp_path):
    from acmg_classifier.criteria.benign.manual import ManualBenignEvaluator
    res = ManualBenignEvaluator(_cfg(tmp_path)).evaluate(_V, _ann(aa="P/S", pos=295))
    bs3 = next(r for r in res if r.criterion == C.BS3)
    assert bs3.triggered and bs3.strength == S.STRONG
    # other genes: BS3 stays manual-only
    res = ManualBenignEvaluator(_cfg(tmp_path)).evaluate(_V, _ann(gene="BRCA1", aa="P/S", pos=295))
    assert not next(r for r in res if r.criterion == C.BS3).triggered


def _func(c, s):
    return CriteriaResult.met(c, s, f"{SRC_TP53_FUNCTIONAL} test")


def test_caveat_splice_pp3_suppresses():
    from acmg_classifier.criteria.registry import _apply_tp53_functional_caveats
    res = [_func(C.BS3, S.STRONG),
           CriteriaResult.met(C.PP3, S.MODERATE, "SpliceAI max_delta=0.5 (Moderate) — missense with predicted splice impact")]
    _apply_tp53_functional_caveats(res, _ann())
    assert res[0].suppressed


@pytest.mark.parametrize("pvs1,expect_supp,expect_strength", [
    (S.VERY_STRONG, True, S.STRONG), (S.STRONG, False, S.MODERATE), (S.MODERATE, False, S.STRONG),
])
def test_caveat_pvs1(pvs1, expect_supp, expect_strength):
    from acmg_classifier.criteria.registry import _apply_tp53_functional_caveats
    res = [_func(C.PS3, S.STRONG), CriteriaResult.met(C.PVS1, pvs1, "splice")]
    _apply_tp53_functional_caveats(res, _ann())
    assert res[0].suppressed is expect_supp and res[0].strength == expect_strength


def test_caveat_leaves_curated_entries():
    from acmg_classifier.criteria.registry import _apply_tp53_functional_caveats
    res = [CriteriaResult.met(C.PS3, S.STRONG, "[manual] curator"),
           CriteriaResult.met(C.PVS1, S.VERY_STRONG, "splice")]
    _apply_tp53_functional_caveats(res, _ann())
    assert not res[0].suppressed


def test_off_by_default(tmp_path):
    """Opt-in (non-commercial data terms): without the flag the TP53 data are not used."""
    from acmg_classifier.config import Config
    from acmg_classifier.criteria.benign.manual import ManualBenignEvaluator
    assert Config(data_dir=Path("./data")).use_tp53_functional is False
    cfg = _cfg(tmp_path)
    cfg.use_tp53_functional = False
    res = ManualBenignEvaluator(cfg).evaluate(_V, _ann(aa="P/S", pos=295))
    assert not next(r for r in res if r.criterion == C.BS3).triggered


def test_cli_flag(tmp_path):
    from click.testing import CliRunner
    from acmg_classifier.cli import cli
    out = CliRunner().invoke(cli, ["classify", "--help"]).output
    assert "--with-tp53-functional" in out
