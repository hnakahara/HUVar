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
    "gene_symbol\tbs2\tinheritance\tbs2_count\tbs2_female_only\n"
    "CDH1\tapplicable\tAD\t10\t\n"
    "LDLR\tapplicable\tAD\t3\t\n"
    "GENE0\tapplicable\tAD\t\t\n"   # no per-gene count → global default
    "TP53\tapplicable\tAD\t8\t1\n"  # counts only females (gnomAD AC_XX), bar >=8
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


class TestBS2FemaleOnly:
    def test_loader_reads_female_only(self, tmp_path):
        spec = BS2Applicability(_tsv(tmp_path))
        assert spec.female_only("TP53") is True
        assert spec.female_only("CDH1") is False

    def test_male_carriers_excluded(self, tmp_path):
        # 9 carriers across both sexes (the real TP53 R248W case) but only 5 are
        # female (AC_XX): below the >=8 female bar → BS2 must NOT fire.
        ev = BS2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("TP53", ac=9, nhomalt=0, ac_xx=5, nhomalt_xx=0))
        assert not r.triggered

    def test_enough_female_carriers_met(self, tmp_path):
        ev = BS2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("TP53", ac=20, nhomalt=0, ac_xx=9, nhomalt_xx=0))
        assert r.triggered
        assert "female" in r.evidence.lower()

    def test_female_homozygotes_subtracted(self, tmp_path):
        # AC_XX counts alleles; a homozygous female is 2 alleles but 1 carrier.
        # female carriers = AC_XX - nhomalt_XX = 9 - 1 = 8 → exactly meets >=8.
        ev = BS2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("TP53", ac=30, nhomalt=1, ac_xx=9, nhomalt_xx=1))
        assert r.triggered

    def test_missing_female_counts_withholds_bs2(self, tmp_path):
        # An older gnomAD DB lacks AC_XX (ac_xx is None) → cannot confirm the
        # female count, so BS2 is withheld rather than counting both sexes.
        ev = BS2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("TP53", ac=20, nhomalt=0))
        assert not r.triggered
        assert "AC_XX" in r.evidence
