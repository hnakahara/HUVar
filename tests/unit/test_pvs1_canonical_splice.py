"""PVS1 for canonical +/-1,2 splice variants: splicing disruption is presumed
(ClinGen SVI PVS1 decision tree), independent of the splice-predictor score."""
from unittest.mock import MagicMock

from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, SpliceScore
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength as S
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pvs1.decision_tree import _splice_branch


def _pc(cq=ConsequenceType.SPLICE_ACCEPTOR):
    return ConsequenceInfo(transcript_id="NM_000321.3", gene_id="ENSG00000139687", gene_symbol="RB1",
                           consequence=cq, biotype="protein_coding", is_mane_select=True, intron="3/26")


def _run(score, tool="openspliceai", alt_rescue=False, cq=ConsequenceType.SPLICE_ACCEPTOR):
    pc = _pc(cq)
    ann = AnnotationData(consequences=[pc],
                         splice=SpliceScore(tool=tool, max_delta=score, is_available=True))
    v = VariantRecord(chrom="chr13", pos=48345079, ref="G", alt="A", assembly=Assembly.GRCH38)
    return _splice_branch(v, ann, MagicMock(), True, alt_rescue, pc)


def test_low_openspliceai_score_still_very_strong():
    # RB1 c.381-1G>A: OpenSpliceAI (80-nt model) 0.06 -> previously Moderate.
    s, ev = _run(0.06)
    assert s == S.VERY_STRONG
    assert "presumed" in ev and "OpenSpliceAI=0.060" in ev


def test_donor_and_high_score_very_strong():
    assert _run(0.98, tool="spliceai", cq=ConsequenceType.SPLICE_DONOR)[0] == S.VERY_STRONG


def test_rescue_transcript_gives_strong():
    assert _run(0.06, alt_rescue=True)[0] == S.STRONG


def test_no_predictor_very_strong():
    pc = _pc()
    ann = AnnotationData(consequences=[pc])
    v = VariantRecord(chrom="chr13", pos=48345079, ref="G", alt="A", assembly=Assembly.GRCH38)
    assert _splice_branch(v, ann, MagicMock(), True, False, pc)[0] == S.VERY_STRONG


def test_lof_mechanism_still_required():
    pc = _pc()
    ann = AnnotationData(consequences=[pc], splice=SpliceScore(tool="spliceai", max_delta=1.0, is_available=True))
    v = VariantRecord(chrom="chr13", pos=48345079, ref="G", alt="A", assembly=Assembly.GRCH38)
    assert _splice_branch(v, ann, MagicMock(), False, False, pc)[0] == S.NOT_MET
