"""PM4 per-gene VCEP applicability gate (disease_prevalence.tsv `pm4` column)."""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pm4 import PM4Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_TSV = (
    "gene_symbol\tpm4\n"
    "BRCA1\tnot_applicable\n"
    "MYH7\tapplicable\n"
)


def _cfg(tmp_path: Path) -> MagicMock:
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_TSV, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = p
    return cfg


def _ann(gene: str) -> AnnotationData:
    # No repeat record → the in-frame deletion is outside a repeat (PM4 eligible).
    return AnnotationData(
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.INFRAME_DELETION, biotype="protein_coding",
        )],
    )


def _snv() -> VariantRecord:
    return VariantRecord(chrom="chr1", pos=100, ref="ACT", alt="A", assembly=Assembly.GRCH38)


def test_pm4_withheld_for_not_applicable_gene(tmp_path):
    # In-frame deletion in BRCA1: VCEP declined PM4 → must NOT fire.
    r = PM4Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("BRCA1"))
    assert not r.triggered
    assert "not applicable" in r.evidence.lower()


def test_pm4_fires_for_applicable_gene(tmp_path):
    r = PM4Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("MYH7"))
    assert r.triggered


def test_pm4_fires_for_uncovered_gene(tmp_path):
    # A gene absent from the table is not gated (fires on the in-frame indel).
    r = PM4Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("GENE0"))
    assert r.triggered
