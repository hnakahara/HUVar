"""PVS1 alternative-transcript rescue: only transcripts of the gene under
evaluation count. VEP also reports transcripts of neighbouring genes within
5 kb (e.g. RB1-DT next to RB1); those must not down-grade PVS1."""
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import ConsequenceType
from acmg_classifier.pvs1.transcript_evaluator import has_alternative_transcript_rescue


def _c(gene, ctype, *, mane=False, canonical=False, biotype="protein_coding", tx="NM_000321.3"):
    return ConsequenceInfo(
        transcript_id=tx, gene_id=f"G_{gene}", gene_symbol=gene, consequence=ctype,
        biotype=biotype, is_mane_select=mane, is_canonical=canonical, exon="1/27",
    )


def _ann(*cs):
    return AnnotationData(consequences=list(cs), gnomad=GnomADData())


def test_neighbouring_gene_does_not_rescue():
    ann = _ann(
        _c("RB1", ConsequenceType.FRAMESHIFT, mane=True, canonical=True),
        _c("RB1-DT", ConsequenceType.UPSTREAM, canonical=True, biotype="lncRNA", tx="NR_1"),
        _c("LPAR6", ConsequenceType.DOWNSTREAM, mane=True, canonical=True, tx="NM_2"),
    )
    assert has_alternative_transcript_rescue(ann) is False


def test_same_gene_coding_transcript_rescues():
    ann = _ann(
        _c("GENE1", ConsequenceType.STOP_GAINED, mane=True, canonical=True, tx="NM_10"),
        _c("GENE1", ConsequenceType.MISSENSE, canonical=True, tx="ENST10"),
    )
    assert has_alternative_transcript_rescue(ann) is True


def test_same_gene_noncoding_transcript_does_not_rescue():
    ann = _ann(
        _c("GENE1", ConsequenceType.STOP_GAINED, mane=True, canonical=True, tx="NM_10"),
        _c("GENE1", ConsequenceType.INTRON, canonical=True, biotype="retained_intron", tx="ENST11"),
    )
    assert has_alternative_transcript_rescue(ann) is False


def test_same_gene_minor_isoform_ignored():
    ann = _ann(
        _c("GENE1", ConsequenceType.STOP_GAINED, mane=True, canonical=True, tx="NM_10"),
        _c("GENE1", ConsequenceType.MISSENSE, tx="ENST12"),
    )
    assert has_alternative_transcript_rescue(ann) is False
