"""VCEP "G to non-G at the last exon nucleotide" PVS1 rules (ATM / MLH1 / CDH1)."""
import pytest

from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType as C, CriterionStrength as S
from acmg_classifier.pvs1.last_nt import evaluate_last_nt


def _pc(gene, tx, hgvs_c, cons=C.MISSENSE):
    return ConsequenceInfo(
        transcript_id=tx, gene_id="ENSG", gene_symbol=gene, consequence=cons,
        biotype="protein_coding", is_mane_select=True, hgvs_c=f"{tx}:{hgvs_c}",
    )


def _fetch(seq):
    return lambda chrom, start, end: seq


class TestAtm:
    def test_frameshifting_exon_strong(self):
        # c.3284 ends a 131-nt exon (c.3154-3284): skipping shifts the frame
        s, ev = evaluate_last_nt(_pc("ATM", "NM_000051.4", "c.3284G>A"))
        assert s == S.STRONG and "one level down" in ev

    def test_inframe_heat_exon_moderate(self):
        # c.3576 ends an in-frame HEAT-repeat exon (list B donor: Strong) -> Moderate
        s, _ = evaluate_last_nt(_pc("ATM", "NM_000051.4", "c.3576G>A", C.SYNONYMOUS))
        assert s == S.MODERATE

    def test_gtrrgt_intron_excluded(self):
        pc = _pc("ATM", "NM_000051.4", "c.3284G>C")
        assert evaluate_last_nt(pc, "chr11", 100, "G", _fetch("GTAAGT")) is None

    def test_non_consensus_intron_applies(self):
        pc = _pc("ATM", "NM_000051.4", "c.3284G>C")
        s, ev = evaluate_last_nt(pc, "chr11", 100, "G", _fetch("GTATTT"))
        assert s == S.STRONG and "GTATTT" in ev

    def test_minus_strand_reverse_complement(self):
        # genomic C (transcript G); + strand bases ACTTAC upstream -> rc GTAAGT
        pc = _pc("ATM", "NM_000051.4", "c.3284G>C")
        assert evaluate_last_nt(pc, "chr11", 100, "C", _fetch("ACTTAC")) is None


def test_mlh1_strong():
    s, _ = evaluate_last_nt(_pc("MLH1", "NM_000249.4", "c.1667G>A"))
    assert s == S.STRONG


@pytest.mark.parametrize("h", ["c.1008G>A", "c.1008G>T", "c.1137G>A", "c.1137G>T"])
def test_cdh1_moderate_no_intron_check(h):
    s, _ = evaluate_last_nt(_pc("CDH1", "NM_004360.5", h), "chr16", 1, "G", _fetch("GTAAGT"))
    assert s == S.MODERATE


@pytest.mark.parametrize("gene,tx,h", [
    ("CDH1", "NM_004360.5", "c.1565C>T"),     # not a G
    ("ATM", "NM_000051.4", "c.3283G>A"),      # not the last exon base
    ("BRCA1", "NM_007294.4", "c.4096G>A"),    # gene without the rule
    ("MLH1", "NM_000249.4", "c.1667G>G"),
])
def test_not_applicable(gene, tx, h):
    assert evaluate_last_nt(_pc(gene, tx, h)) is None


def test_evaluator_hook_missense_gets_pvs1(tmp_path):
    """The PVS1 evaluator applies the rule to a VEP missense call; the reference
    FASTA is read through the evaluator (here a tiny FASTA on the + strand)."""
    from unittest.mock import MagicMock
    import pysam
    from acmg_classifier.criteria.pathogenic.pvs1 import PVS1Evaluator
    from acmg_classifier.models.annotation import AnnotationData, GnomADData
    from acmg_classifier.models.enums import Assembly
    from acmg_classifier.models.variant import VariantRecord

    fa = tmp_path / "g.fa"
    fa.write_text(">chr3\n" + "A" * 99 + "G" + "GTAAGT" + "A" * 20 + "\n")
    pysam.faidx(str(fa))
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text("gene_symbol\tpvs1\nMLH1\tapplicable\n", encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = p
    cfg.genome_fasta = fa
    ev = PVS1Evaluator(cfg)
    ann = AnnotationData(consequences=[_pc("MLH1", "NM_000249.4", "c.1667G>A")],
                         gnomad=GnomADData())
    v = VariantRecord(chrom="chr3", pos=100, ref="G", alt="A", assembly=Assembly.GRCH38)
    # GTAAGT follows -> consensus donor -> rule not met (falls through: not LoF)
    assert not ev.evaluate(v, ann).triggered
    fa.write_text(">chr3\n" + "A" * 99 + "G" + "GTGAGC" + "A" * 20 + "\n")
    pysam.faidx(str(fa))
    ev2 = PVS1Evaluator(cfg)
    r = ev2.evaluate(v, ann)
    assert r.triggered and r.strength == S.STRONG and "GTGAGC" in r.evidence
