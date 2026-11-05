"""Unit tests for variant left-alignment normalizer."""
from pathlib import Path
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


# --- normalisation against a reference (vt algorithm) --------------------------
import pysam  # noqa: E402
import pytest  # noqa: E402

from acmg_classifier.annotation.normalizer import normalize_alleles, normalize_record  # noqa: E402

#            1234567890123456789
_CHR1 = "CCTTGAAAACGTACGTTTTT" + "ACGT" * 20


@pytest.fixture()
def fasta(tmp_path):
    p = tmp_path / "ref.fa"
    p.write_text(">chr1\n" + _CHR1 + "\n")
    pysam.faidx(str(p))
    fa = pysam.FastaFile(str(p))
    yield fa
    fa.close()


def test_padded_deletion_from_a_split_multiallelic_record(fasta):
    # GATK "GAAAA>GAAA,G" split -> GAAAA>GAAA: one A of the poly-A deleted (pos 5 = G)
    assert normalize_alleles("chr1", 5, "GAAAA", "GAAA", fasta) == (5, "GA", "G")


def test_deletion_inside_a_repeat_is_left_aligned(fasta):
    # the same deletion written at the last A of the run
    assert normalize_alleles("chr1", 8, "AA", "A", fasta) == (5, "GA", "G")
    assert normalize_alleles("1", 8, "AA", "A", fasta) == (5, "GA", "G")   # bare contig


def test_insertion_in_a_repeat_is_left_aligned(fasta):
    assert normalize_alleles("chr1", 9, "A", "AA", fasta) == (5, "G", "GA")
    assert normalize_alleles("chr1", 16, "TT", "TTT", fasta) == (15, "G", "GT")


def test_snv_mnv_and_already_normal(fasta):
    assert normalize_alleles("chr1", 3, "T", "C", fasta) == (3, "T", "C")
    assert normalize_alleles("chr1", 5, "GA", "G", fasta) == (5, "GA", "G")
    # MNV with shared bases collapses to the differing base
    assert normalize_alleles("chr1", 10, "CGT", "CAT", fasta) == (11, "G", "A")


def test_without_reference_or_with_wrong_ref_only_trims(fasta):
    assert normalize_alleles("chr1", 8, "AAC", "AC", None) == (8, "AA", "A")
    # REF does not match the reference (other assembly): trimmed, not shifted
    assert normalize_alleles("chr1", 2, "GTTT", "GTT", fasta) == (2, "GT", "G")
    # no left base at the start of the contig: left as it was
    assert normalize_alleles("chr1", 1, "CC", "C", fasta) == (1, "CC", "C")


def test_normalize_record_keeps_the_input(fasta):
    v = VariantRecord(chrom="chr1", pos=5, ref="GAAAA", alt="GAAA", assembly=Assembly.GRCH38)
    n = normalize_record(v, fasta)
    assert (n.pos, n.ref, n.alt, n.key) == (5, "GA", "G", "chr1:5:GA:G")
    assert n.variant_id == n.key
    assert n.input_variant == "chr1:5:GAAAA:GAAA"
    same = VariantRecord(chrom="chr1", pos=5, ref="GA", alt="G", assembly=Assembly.GRCH38)
    assert normalize_record(same, fasta) is same


def test_read_vcf_normalises_split_multiallelic_records(tmp_path, fasta):
    from acmg_classifier.io.vcf_reader import read_vcf
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5\t.\tGAAAA\tGAAA,G\t50\tPASS\t.\n"
        "chr1\t3\t.\tT\tC\t50\tPASS\t.\n"
        "chr1\t9\t.\tA\t<DEL>\t50\tPASS\t.\n")
    recs = list(read_vcf(vcf, Assembly.GRCH38, fasta=Path(fasta.filename.decode())))
    keys = [(r.key, r.input_variant) for r in recs]
    assert keys == [("chr1:5:GA:G", "chr1:5:GAAAA:GAAA"),
                    ("chr1:5:GAAAA:G", None),
                    ("chr1:3:T:C", None),
                    ("chr1:9:A:<DEL>", None)]
    # without a FASTA the padding is still trimmed
    assert list(read_vcf(vcf, Assembly.GRCH38))[0].key == "chr1:5:GA:G"


def test_tsv_reports_the_input_representation():
    from acmg_classifier.io.tsv_writer import _HEADER, _build_row
    from acmg_classifier.models.classification import ClassificationResult
    from acmg_classifier.models.enums import Pathogenicity
    r = ClassificationResult(
        variant_id="chr1:5:GA:G", input_variant="chr1:5:GAAAA:GAAA", chrom="chr1", pos=5,
        ref="GA", alt="G", criteria_results=[], classification_2015=Pathogenicity.VUS,
        classification_2015_rules="", bayesian_score=0,
        classification_bayesian=Pathogenicity.VUS, warnings=["other"])
    row = dict(zip(_HEADER, _build_row(r)))
    assert row["input_variant"] == "chr1:5:GAAAA:GAAA"
    assert row["warnings"].startswith("NORMALIZED: VCF allele chr1:5:GAAAA:GAAA")
    assert row["warnings"].endswith("; other")
