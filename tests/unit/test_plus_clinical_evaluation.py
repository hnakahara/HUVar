"""Evaluating the criteria on a MANE Plus Clinical transcript (phase 2).

The criteria are evaluated on the Plus Clinical consequence (consequence type,
exon structure, isoform in-silico scores), while VCEP position rules and ClinVar
codon comparisons keep MANE Select numbering (``rule_consequence``) — the same
genomic residue. MECP2: MANE Select NM_001110792 (e1) numbers residues 12 higher
than the Plus Clinical NM_004992 (e2) in the shared exons.
"""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.annotation.orchestrator import AnnotationOrchestrator
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.pathogenic.pm1 import PM1Evaluator
from acmg_classifier.criteria.registry import CriteriaRegistry
from acmg_classifier.config import Config
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, ESM1bData, GnomADData, RevelData, TranscriptInSilico,
)
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pipeline.pipeline import (
    criteria_changes, plus_clinical_evaluations, transcript_note,
)

S = CriterionStrength
SEL, SEL_E = "NM_001110792.2", "ENST00000453960.7"
PC, PC_E = "NM_004992.4", "ENST00000303391.11"
_PM1_TSV = Path(__file__).resolve().parents[2] / "resources" / "shared" / "pm1_hotspots.tsv"
_V = VariantRecord(chrom="chrX", pos=154030000, ref="G", alt="A", assembly=Assembly.GRCH38)


def _c(tx, pos, csq=ConsequenceType.MISSENSE, sel=False, pc=False, aa="R/C"):
    return ConsequenceInfo(
        transcript_id=tx, gene_id="ENSG", gene_symbol="MECP2", consequence=csq,
        biotype="protein_coding", is_mane_select=sel, is_mane_plus_clinical=pc,
        protein_position=pos, amino_acids=aa if pos else None,
        hgvs_c=f"{tx}:c.{pos * 3 if pos else '1-10'}G>A",
    )


def _ann(sel_pos=110, pc_pos=98, sel_csq=ConsequenceType.MISSENSE):
    return AnnotationData(
        gnomad=GnomADData(),
        consequences=[_c(SEL, sel_pos, sel_csq, sel=True), _c(SEL_E, sel_pos, sel_csq, sel=True),
                      _c(PC, pc_pos, pc=True), _c(PC_E, pc_pos, pc=True)],
        esm1b=ESM1bData(llr=-5.0),
        transcript_insilico={
            "NM_001110792": TranscriptInSilico(esm1b=ESM1bData(llr=-5.0)),
            "NM_004992": TranscriptInSilico(esm1b=ESM1bData(llr=-12.0), revel=RevelData(score=0.9)),
        },
    )


def test_default_is_unchanged():
    ann = _ann()
    assert ann.primary_consequence.transcript_id == SEL
    assert ann.rule_consequence is ann.primary_consequence
    assert ann.rule_view() is ann
    assert transcript_note(ann) == []


def test_for_transcript_switches_primary_scores_and_keeps_rule_numbering():
    ann = _ann().for_transcript(PC)
    assert ann.primary_consequence.transcript_id == PC
    assert ann.rule_consequence.transcript_id == SEL
    assert ann.esm1b.llr == -12.0 and ann.revel.score == 0.9
    back = ann.for_transcript(None)
    assert back.primary_consequence.transcript_id == SEL and back.esm1b.llr == -5.0
    assert ann.rule_view().primary_consequence.transcript_id == SEL
    # Plus Clinical summary stays anchored on MANE Select
    assert [m.transcript_id for m in ann.mane_plus_clinical_summary()] == [PC]
    # version-less / Ensembl accession resolves too
    assert _ann().for_transcript("NM_004992").primary_consequence.transcript_id == PC


def test_unknown_transcript_falls_back_with_note():
    ann = _ann().for_transcript("NM_999999.1")
    assert ann.primary_consequence.transcript_id == SEL
    assert "no VEP consequence" in transcript_note(ann)[0]


def test_registry_gives_rule_view_to_rule_numbered_evaluators():
    seen = {}

    class _Rec(CriterionEvaluator):
        def __init__(self, crit, rule):
            self.crit, self.rule_numbered = crit, rule

        def evaluate(self, variant, annotation, supplement=None):
            seen[self.crit] = annotation.primary_consequence.transcript_id
            return CriteriaResult.not_met(self.crit, "x")

    reg = CriteriaRegistry(Config(data_dir=Path("./data")))
    reg._evaluators = [_Rec(ACMGCriterion.PM1, True), _Rec(ACMGCriterion.PP3, False)]
    reg.evaluate_all(_V, _ann().for_transcript(PC))
    assert seen == {ACMGCriterion.PM1: SEL, ACMGCriterion.PP3: PC}


def _pm1(ann):
    cfg = MagicMock()
    cfg.pm1_hotspots_tsv = _PM1_TSV
    return PM1Evaluator(cfg).evaluate(_V, ann)


def test_pm1_reads_mane_select_numbering():
    # MECP2 PM1 102-174 (MANE Select). e1 residue 110 == e2 residue 98.
    ann = _ann(sel_pos=110, pc_pos=98).for_transcript(PC)
    assert _pm1(ann.rule_view()).triggered          # what the registry passes
    assert not _pm1(ann).triggered                  # e2 numbering would miss it


def test_plus_clinical_only_region_cannot_meet_position_rules():
    ann = _ann(sel_pos=None, pc_pos=120, sel_csq=ConsequenceType.INTRON).for_transcript(PC)
    assert ann.primary_consequence.consequence == ConsequenceType.MISSENSE
    assert not _pm1(ann.rule_view()).triggered
    assert "cannot be met" in transcript_note(ann)[0]


def test_criteria_changes():
    base = [CriteriaResult.met(ACMGCriterion.PM1, S.MODERATE, "a"),
            CriteriaResult.not_met(ACMGCriterion.PP3, "b")]
    other = [CriteriaResult.not_met(ACMGCriterion.PM1, "a"),
             CriteriaResult.met(ACMGCriterion.PP3, S.STRONG, "b")]
    assert criteria_changes(base, other) == "PM1:Moderate->NotMet; PP3:NotMet->Strong"
    assert criteria_changes(base, base) == ""


def test_plus_clinical_evaluations_fill_classification():
    from acmg_classifier.models.enums import Pathogenicity as P
    ann = _ann()
    base = [CriteriaResult.not_met(ACMGCriterion.PP3, "b")]
    seen = []

    def evaluate(a):
        seen.append(a.primary_consequence.transcript_id)
        return ([CriteriaResult.met(ACMGCriterion.PP3, S.STRONG, "b")],
                P.LIKELY_PATHOGENIC, "PP3", 4, P.VUS)

    (info,) = plus_clinical_evaluations(ann, Config(data_dir=Path("./data")), base, evaluate)
    assert seen == [PC]
    assert info.classification_2015 == P.LIKELY_PATHOGENIC.value and info.bayesian_score == 4
    assert info.criteria_changes == "PP3:NotMet->Strong"
    off = Config(data_dir=Path("./data"), evaluate_plus_clinical=False)
    (info,) = plus_clinical_evaluations(ann, off, base, evaluate)
    assert info.classification_2015 is None and len(seen) == 1


def test_orchestrator_records_isoform_scores():
    orch = object.__new__(AnnotationOrchestrator)
    orch._esm1b_db = MagicMock()
    orch._esm1b_db.lookup.side_effect = lambda uid, pos, alt: ESM1bData(llr=float(-pos))
    orch._revel_reader = object()
    cons = _ann().consequences
    cons[1].uniprot_id, cons[3].uniprot_id = "P51608-2", "P51608"
    calls = []

    def fake_revel(reader, chrom, pos, ref, alt, aa=None, transcripts=None):
        calls.append((aa, set(transcripts or ())))
        return RevelData(score=0.5)

    out = orch._plus_clinical_insilico(_V, cons[0], cons, ESM1bData(llr=-1.0), None, fake_revel)
    assert out["NM_001110792"].esm1b.llr == -1.0
    assert out["NM_004992"].esm1b.llr == -98.0        # e2 residue 98, e2 isoform
    assert out["NM_004992"].revel.score == 0.5
    assert calls[0] == (("R", "C"), {"ENST00000303391"})
    # No Plus Clinical transcript -> nothing extra
    assert orch._plus_clinical_insilico(_V, cons[0], cons[:2], None, None, fake_revel) == {}
