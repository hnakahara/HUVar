"""BS2 for genes without a VCEP rule: inheritance from gene_inheritance.tsv,
dominant carriers counted only in LoF-constrained genes (LOEUF < 0.5)."""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.config import Config
from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.local_db.inheritance_db import load_inheritance_map

_DP = "gene_symbol\tbs2\tinheritance\tbs2_count\tbs2_female_only\tbs2_hom_only\nVCEPGENE\tapplicable\tAD\t\t\t\n"
_INH = "gene\tinheritance\nRECGENE\tAR\nDOMGENE\tAD\nMIXGENE\tAD/AR\nXGENE\tXL\nMTGENE\tMT\n"


def _cfg(tmp_path: Path, max_loeuf=0.5) -> MagicMock:
    load_inheritance_map.cache_clear()
    dp = tmp_path / "disease_prevalence.tsv"
    dp.write_text(_DP, encoding="utf-8")
    inh = tmp_path / "gene_inheritance.tsv"
    inh.write_text(_INH, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = dp
    cfg.gene_inheritance_tsv = inh
    cfg.bs2_min_homalt, cfg.bs2_min_hemi, cfg.bs2_min_het = 2, 2, 3
    cfg.bs2_ad_max_loeuf = max_loeuf
    return cfg


def _ann(gene: str, **gd) -> AnnotationData:
    return AnnotationData(
        gnomad=GnomADData(filter_pass=True, **gd),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding")],
    )


_V = VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


def test_default_threshold():
    assert Config(data_dir=Path("./data")).bs2_ad_max_loeuf == 0.5


def test_recessive_homozygotes(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path))
    assert ev.evaluate(_V, _ann("RECGENE", ac=50, nhomalt=2)).triggered
    # carriers alone never count for a recessive gene
    assert not ev.evaluate(_V, _ann("RECGENE", ac=50, nhomalt=1, loeuf=0.1)).triggered


def test_dominant_requires_lof_constraint(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path))
    assert ev.evaluate(_V, _ann("DOMGENE", ac=5, nhomalt=0, loeuf=0.3)).triggered
    r = ev.evaluate(_V, _ann("DOMGENE", ac=50, nhomalt=0, loeuf=0.8))
    assert not r.triggered and "LoF-constrained" in r.evidence
    assert not ev.evaluate(_V, _ann("DOMGENE", ac=50, nhomalt=0, loeuf=None)).triggered
    # boundary: LOEUF equal to the cutoff is not constrained
    assert not ev.evaluate(_V, _ann("DOMGENE", ac=50, nhomalt=0, loeuf=0.5)).triggered
    # constrained but too few carriers
    assert not ev.evaluate(_V, _ann("DOMGENE", ac=2, nhomalt=0, loeuf=0.3)).triggered


def test_dominant_homozygotes_do_not_count(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path))
    assert not ev.evaluate(_V, _ann("DOMGENE", ac=10, nhomalt=3, loeuf=0.9)).triggered


def test_mixed_ad_ar(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path))
    assert ev.evaluate(_V, _ann("MIXGENE", ac=10, nhomalt=2, loeuf=0.9)).triggered
    assert not ev.evaluate(_V, _ann("MIXGENE", ac=10, nhomalt=0, loeuf=0.9)).triggered
    assert ev.evaluate(_V, _ann("MIXGENE", ac=10, nhomalt=0, loeuf=0.2)).triggered


def test_x_linked_and_unknown(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path))
    assert ev.evaluate(_V, _ann("XGENE", ac=5, nhemi=2)).triggered
    assert ev.evaluate(_V, _ann("UNLISTED", ac=5, nhomalt=2)).triggered
    assert not ev.evaluate(_V, _ann("UNLISTED", ac=50, nhomalt=0, loeuf=0.1)).triggered
    assert ev.evaluate(_V, _ann("MTGENE", ac=5, nhomalt=2)).triggered


def test_vcep_gene_unchanged(tmp_path):
    # VCEP-applicable AD gene keeps the VCEP carrier route regardless of LOEUF
    ev = BS2Evaluator(_cfg(tmp_path))
    assert ev.evaluate(_V, _ann("VCEPGENE", ac=5, nhomalt=0, loeuf=0.9)).triggered


def test_custom_cutoff(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path, max_loeuf=1.0))
    assert ev.evaluate(_V, _ann("DOMGENE", ac=5, nhomalt=0, loeuf=0.8)).triggered
