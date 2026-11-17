"""PVS1 for canonical +/-1,2 splice variants: splicing disruption is presumed
(ClinGen SVI PVS1 decision tree), independent of the splice-predictor score."""
from unittest.mock import MagicMock

from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, SpliceScore
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength as S
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pvs1.decision_tree import _splice_branch


def _pc(cq=ConsequenceType.SPLICE_ACCEPTOR):
    return ConsequenceInfo(transcript_id="NM_000321.3", gene_id="ENSG00000139687", gene_symbol="RB1",
                           consequence=cq, biotype="protein_coding", is_mane_select=True, intron="2/26")


def _run(score, tool="openspliceai", alt_rescue=False, cq=ConsequenceType.SPLICE_ACCEPTOR,
         intron=None):
    pc = _pc(cq)
    if intron is not None:
        pc.intron = intron
    ann = AnnotationData(consequences=[pc],
                         splice=SpliceScore(tool=tool, max_delta=score, is_available=True))
    v = VariantRecord(chrom="chr13", pos=48345079, ref="G", alt="A", assembly=Assembly.GRCH38)
    return _splice_branch(v, ann, MagicMock(), True, alt_rescue, pc)


def test_low_openspliceai_score_still_very_strong():
    # RB1 intron-2 acceptor (skips the 116-nt exon 3: frameshift -> NMD):
    # OpenSpliceAI (80-nt model) 0.06 -> previously Moderate.
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


def test_intronic_deletion_beyond_canonical_not_pvs1():
    # LZTR1 c.651+10_651+46del: VEP splice_donor_variant (left-aligned allele in a
    # repeat), but the HGVS (3') change does not touch +1/+2 -> no PVS1.
    pc = ConsequenceInfo(transcript_id="NM_006767.4", gene_id="ENSG", gene_symbol="LZTR1",
                         consequence=ConsequenceType.SPLICE_DONOR, biotype="protein_coding",
                         is_mane_select=True, hgvs_c="NM_006767.4:c.651+10_651+46del")
    ann = AnnotationData(consequences=[pc])
    v = VariantRecord(chrom="chr22", pos=20989676, ref="GG", alt="G", assembly=Assembly.GRCH38)
    s, ev = _splice_branch(v, ann, MagicMock(), True, False, pc)
    assert s == S.NOT_MET and "PP3" in ev


def test_intronic_beyond_canonical_helper():
    from acmg_classifier.pvs1.decision_tree import intronic_beyond_canonical as f
    assert f("c.531+5_531+8del") and f("c.651+10_651+46del") and f("c.100-3C>G")
    assert not f("c.381-1G>A") and not f("c.539+2T>C") and not f("c.101-3_101-1del")
    assert not f("c.1135_1139+3del")          # crosses the exon/intron boundary
    assert not f("c.3756+905_3897+151del")    # removes a whole exon
    assert not f("c.1696del") and not f(None)  # exonic


def test_pp3_scores_noncanonical_splice_annotation(tmp_path):
    from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
    from acmg_classifier.config import Config
    pc = ConsequenceInfo(transcript_id="NM_x", gene_id="ENSG", gene_symbol="SLC6A8",
                         consequence=ConsequenceType.SPLICE_DONOR, biotype="protein_coding",
                         is_mane_select=True, hgvs_c="c.644+3_644+6del")
    ann = AnnotationData(consequences=[pc],
                         splice=SpliceScore(tool="spliceai", max_delta=0.99, is_available=True))
    v = VariantRecord(chrom="chrX", pos=1, ref="AG", alt="A", assembly=Assembly.GRCH38)
    r = PP3Evaluator(Config()).evaluate(v, ann)
    assert r.triggered, r.evidence


# --- SVI exon-skipping arm (in-frame / stop-codon exon / penultimate exon) ---

def test_inframe_exon_skip_follows_10pct_rule():
    # RB1 c.381-1G>A (intron 3 acceptor) skips exon 4: 120 nt, in-frame, 4.3% of
    # the protein, role unknown -> PVS1_Moderate.
    s, ev = _run(0.9, intron="3/26")
    assert s == S.MODERATE
    assert "in-frame skip of exon 4/27" in ev and "10% rule" in ev


def test_inframe_skip_with_rescue_transcript_not_above_moderate():
    assert _run(0.9, intron="3/26", alt_rescue=True)[0] == S.MODERATE


def test_skip_of_stop_codon_exon_uses_10pct_rule():
    # RB1 intron-26 acceptor skips exon 27 (last; holds the stop codon).
    s, ev = _run(0.9, intron="26/26")
    assert s == S.MODERATE and "holding the stop codon" in ev


def test_out_of_frame_penultimate_exon_skip_escapes_nmd():
    # RB1 intron-26 donor skips exon 26 (50 nt, out of frame): the frameshift
    # runs into the last exon -> NMD not expected -> 10% rule (Moderate).
    s, ev = _run(0.9, intron="26/26", cq=ConsequenceType.SPLICE_DONOR)
    assert s == S.MODERATE and "penultimate exon 26/27" in ev


def test_non_mane_transcript_keeps_very_strong():
    pc = _pc()
    pc.transcript_id = "NM_999999.1"
    pc.intron = "3/26"
    ann = AnnotationData(consequences=[pc])
    v = VariantRecord(chrom="chr13", pos=48345079, ref="G", alt="A", assembly=Assembly.GRCH38)
    assert _splice_branch(v, ann, MagicMock(), True, False, pc)[0] == S.VERY_STRONG
