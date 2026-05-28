"""Unit tests for variant left-alignment normalizer."""
from acmg_classifier.models.enums import Assembly
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.annotation.normalizer import left_align_and_trim


def test_snv_passthrough(tmp_path):
    v = VariantRecord(chrom="chr17", pos=100, ref="G", alt="A", assembly=Assembly.GRCH38)
    # SNVs pass through without touching the FASTA
    result = left_align_and_trim(v, tmp_path / "fake.fa")
    assert result.ref == "G"
    assert result.alt == "A"
    assert result.pos == 100
