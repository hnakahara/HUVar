"""BS2 per-gene VCEP observation counts + lowered global defaults."""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.config import Config
from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
from acmg_classifier.criteria.bs2_genes import BS2Applicability
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

# CDH1: cancer panel requiring >=10 unaffected individuals (high bar). LDLR: 3.
_TSV = (
    "gene_symbol\tbs2\tinheritance\tbs2_count\n"
    "CDH1\tapplicable\tAD\t10\n"
    "LDLR\tapplicable\tAD\t3\n"
    "GENE0\tapplicable\tAD\t\n"   # no per-gene count → global default
)


def _tsv(tmp_path: Path) -> Path:
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_TSV, encoding="utf-8")
    return p


def _cfg(tmp_path: Path) -> MagicMock:
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _tsv(tmp_path)
    cfg.bs2_min_homalt = 2
    cfg.bs2_min_hemi = 2
    cfg.bs2_min_het = 3
    return cfg


def _ann(gene: str, **gd) -> AnnotationData:
    return AnnotationData(
        gnomad=GnomADData(filter_pass=True, **gd),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        )],
    )


def _snv() -> VariantRecord:
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


def test_config_defaults_lowered():
    c = Config(data_dir=Path("./data"))
    assert (c.bs2_min_homalt, c.bs2_min_hemi, c.bs2_min_het) == (2, 2, 3)


def test_loader_reads_count(tmp_path):
    spec = BS2Applicability(_tsv(tmp_path))
    assert spec.count("CDH1") == 10
    assert spec.count("LDLR") == 3
    assert spec.count("GENE0") is None


class TestBS2CountOverride:
    def test_cancer_gene_below_vcep_count_not_met(self, tmp_path):
        # 6 healthy carriers would fire BS2 under the default het=3, but CDH1's
        # VCEP demands >=10 — must NOT fire (protects a pathogenic cancer variant).
        ev = BS2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("CDH1", ac=6, nhomalt=0))
        assert not r.triggered

    def test_cancer_gene_at_vcep_count_met(self, tmp_path):
        ev = BS2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("CDH1", ac=12, nhomalt=0))
        assert r.triggered

    def test_gene_without_count_uses_lowered_default(self, tmp_path):
        # GENE0 has no per-gene count → default het=3; 4 carriers fires.
        ev = BS2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("GENE0", ac=4, nhomalt=0))
        assert r.triggered
