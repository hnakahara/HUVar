"""GALT VCEP PM2 exception list: listed known pathogenic variants never meet PM2."""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_V = VariantRecord(chrom="9", pos=1, ref="A", alt="G", assembly=Assembly.GRCH38)


def _ann(hgvs_c, hgvs_p):
    return AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_000155.4", gene_id="ENSG", gene_symbol="GALT",
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding", is_mane_select=True,
        hgvs_c=f"NM_000155.4:{hgvs_c}", hgvs_p=f"NP_000146.2:{hgvs_p}")])


def _ev():
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = Path("/nonexistent.tsv")
    return PM2Evaluator(cfg)


def test_listed_variant_not_met_even_if_absent():
    r = _ev().evaluate(_V, _ann("c.563A>G", "p.Gln188Arg"))
    assert not r.triggered and "PM2 exception list" in r.evidence and "c.563A>G" in r.evidence


def test_other_variant_unaffected():
    r = _ev().evaluate(_V, _ann("c.565C>T", "p.Arg189Cys"))   # gnomad None -> absent
    assert r.triggered
