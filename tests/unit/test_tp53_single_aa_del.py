"""TP53 VCEP PP3/BP4 for single amino-acid in-frame deletions from the VCEP
"single aa BayesDel spreadsheet" (resources/shared/tp53_single_aa_del_bayesdel.tsv)."""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
from acmg_classifier.criteria.insilico_genes import TP53Codes
from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, SpliceScore
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, CriterionStrength, InSilicoTool,
)
from acmg_classifier.models.variant import VariantRecord

_BUNDLED = (Path(__file__).resolve().parents[2] / "resources" / "shared"
            / "tp53_single_aa_del_bayesdel.tsv")
_V = VariantRecord(chrom="17", pos=1, ref="ACTG", alt="A", assembly=Assembly.GRCH38)


def _cfg(tmp_path, bayesdel=True):
    cfg = MagicMock()
    cfg.insilico_tool = InSilicoTool.REVEL
    cfg.use_cadd = False
    cfg.use_bayesdel = bayesdel
    cfg.disease_prevalence_tsv = Path("/nonexistent_spec.tsv")
    cfg.tp53_codes_tsv = tmp_path / "missing_codes.tsv"   # deletion table: bundled copy
    return cfg


def _ann(hgvs_p, splice=None, csq=ConsequenceType.INFRAME_DELETION):
    return AnnotationData(
        consequences=[ConsequenceInfo(
            transcript_id="NM_000546.6", gene_id="ENSG", gene_symbol="TP53",
            consequence=csq, biotype="protein_coding",
            hgvs_c="NM_000546.6:c.1_3del", hgvs_p=f"NP_000537.3:{hgvs_p}")],
        splice=(SpliceScore(tool="spliceai", max_delta=splice, is_available=True)
                if splice is not None else None))


def test_bundled_table_loaded():
    codes = TP53Codes(Path("/nonexistent/tp53_pp3_bp4_codes.tsv"))
    assert codes and len(codes._del) > 300
    assert codes.deletion_bayesdel("p.Asp391del") > 0.16          # 0.2045
    assert codes.deletion_bayesdel("p.Ser392del") < 0.16          # 0.0810
    assert codes.deletion_bayesdel("NP_000537.3:p.Leu252del") is not None  # c.-keyed row
    assert _BUNDLED.exists()


def test_pp3_and_bp4(tmp_path):
    cfg = _cfg(tmp_path)
    r = PP3Evaluator(cfg).evaluate(_V, _ann("p.Asp391del"))
    assert r.triggered and r.strength == CriterionStrength.SUPPORTING and "0.205" in r.evidence
    assert not BP4Evaluator(cfg).evaluate(_V, _ann("p.Asp391del")).triggered

    assert not PP3Evaluator(cfg).evaluate(_V, _ann("p.Ser392del")).triggered
    r = BP4Evaluator(cfg).evaluate(_V, _ann("p.Ser392del", splice=0.05))
    assert r.triggered and r.strength == CriterionStrength.SUPPORTING
    # Predicted splice effect: BP4 withheld, PP3 by the splicing rule.
    assert not BP4Evaluator(cfg).evaluate(_V, _ann("p.Ser392del", splice=0.4)).triggered
    r = PP3Evaluator(cfg).evaluate(_V, _ann("p.Ser392del", splice=0.4))
    assert r.triggered and "splice" in r.evidence


def test_scope_and_gate(tmp_path):
    cfg = _cfg(tmp_path)
    # Multi-residue deletion: outside the single-aa rule.
    r = PP3Evaluator(cfg).evaluate(_V, _ann("p.Arg280_Arg282del"))
    assert not r.triggered and "single amino-acid" in r.evidence
    # BayesDel licence gate (--with-bayesdel under REVEL/AlphaMissense) applies.
    r = PP3Evaluator(_cfg(tmp_path, bayesdel=False)).evaluate(_V, _ann("p.Asp391del"))
    assert not r.triggered and r.evidence == "Consequence not applicable for PP3"
