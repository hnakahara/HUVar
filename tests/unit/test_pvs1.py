"""Unit tests for PVS1 decision tree."""
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pvs1.nmd_predictor import predicts_nmd, is_last_exon, is_penultimate_exon
from acmg_classifier.pvs1.transcript_evaluator import gene_has_lof_mechanism


def _consequence(ctype, exon="5/24", gene="BRCA1"):
    return ConsequenceInfo(
        transcript_id="NM_007294.4",
        gene_id="ENSG00000012048",
        gene_symbol=gene,
        consequence=ctype,
        biotype="protein_coding",
        is_mane_select=True,
        exon=exon,
    )


class TestPVS1ApplicabilityGate:
    """A VCEP that declines PVS1 (LoF not the disease mechanism — MYOC,
    RASopathy, cardiomyopathy, …) withholds it even for a null variant."""

    def _ev(self, tmp_path, gene_pvs1):
        from unittest.mock import MagicMock
        from acmg_classifier.criteria.pathogenic.pvs1 import PVS1Evaluator
        p = tmp_path / "disease_prevalence.tsv"
        p.write_text(
            "gene_symbol\tpvs1\n" + "\n".join(f"{g}\t{v}" for g, v in gene_pvs1) + "\n",
            encoding="utf-8",
        )
        cfg = MagicMock()
        cfg.disease_prevalence_tsv = p
        return PVS1Evaluator(cfg)

    def _ann(self, gene):
        c = _consequence(ConsequenceType.STOP_GAINED, exon="5/24", gene=gene)
        return AnnotationData(consequences=[c], gnomad=GnomADData(loeuf=0.10))

    def test_not_applicable_gene_withheld(self, tmp_path):
        from acmg_classifier.models.enums import ACMGCriterion
        ev = self._ev(tmp_path, [("MYOC", "not_applicable")])
        v = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)
        r = ev.evaluate(v, self._ann("MYOC"))
        assert not r.triggered
        assert r.criterion == ACMGCriterion.PVS1
        assert "not applicable" in r.evidence.lower()

    def test_applicable_gene_not_gated(self, tmp_path):
        # An applicable gene is NOT short-circuited (proceeds to the decision tree).
        ev = self._ev(tmp_path, [("MYOC", "not_applicable"), ("BRCA1", "applicable")])
        v = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)
        r = ev.evaluate(v, self._ann("BRCA1"))
        # BRCA1 LOEUF 0.10 + stop-gained → not blocked by the applicability gate
        # (the result then depends on the decision tree / ClinVar caps).
        assert "not applicable" not in r.evidence.lower()

    def test_committed_tsv_marks_gof_genes(self):
        import csv
        from pathlib import Path
        tsv = Path(__file__).resolve().parents[2] / "resources" / "shared" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        # Gain-of-function / dominant-negative genes: PVS1 not applicable.
        for gene in ("MYOC", "BRAF", "KRAS", "PTPN11", "MYH7", "TNNT2",
                     "PIK3CA", "VWF"):
            assert rows[gene]["pvs1"] == "not_applicable", gene
        # A haploinsufficiency gene keeps PVS1 applicable.
        assert rows["BRCA1"]["pvs1"] != "not_applicable"
        # ACTA1/RYR1: the Congenital Myopathies VCEP applies PVS1 to null variants
        # (LoF is a known mechanism) — must NOT be marked not_applicable.
        assert rows["ACTA1"]["pvs1"] == "applicable"
        assert rows["RYR1"]["pvs1"] == "applicable"


class TestPVS1VcepLofEstablished:
    """A gene whose VCEP explicitly applies PVS1 is treated as having an
    established LoF mechanism — the decision tree skips the ClinVar/LOEUF
    heuristic, so a null variant fires PVS1 even without those signals."""

    def _ev(self, tmp_path, gene, value):
        from unittest.mock import MagicMock
        from acmg_classifier.criteria.pathogenic.pvs1 import PVS1Evaluator
        p = tmp_path / "dp.tsv"
        p.write_text(f"gene_symbol\tpvs1\n{gene}\t{value}\n", encoding="utf-8")
        cfg = MagicMock()
        cfg.disease_prevalence_tsv = p
        return PVS1Evaluator(cfg)

    def _ann_no_lof_signal(self, gene):
        # No LOEUF and (with a MagicMock clinvar DB) no P/LP null count → the
        # heuristic alone would NOT establish LoF.
        c = _consequence(ConsequenceType.STOP_GAINED, exon="5/24", gene=gene)
        return AnnotationData(consequences=[c], gnomad=GnomADData(loeuf=None))

    def test_applicable_gene_fires_without_heuristic_signal(self, tmp_path):
        ev = self._ev(tmp_path, "ACTA1", "applicable")
        v = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)
        r = ev.evaluate(v, self._ann_no_lof_signal("ACTA1"))
        assert r.triggered and r.strength == CriterionStrength.VERY_STRONG

    def test_unmarked_gene_still_needs_heuristic(self, tmp_path):
        # A gene with no VCEP PVS1 marking and no LoF signal → not established.
        ev = self._ev(tmp_path, "ZZZ1", "")
        v = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)
        r = ev.evaluate(v, self._ann_no_lof_signal("ZZZ1"))
        assert not r.triggered


class TestNMDPredictor:
    def test_nmd_predicted_early_exon(self):
        c = _consequence(ConsequenceType.STOP_GAINED, exon="5/24")
        assert predicts_nmd(c) is True

    def test_nmd_not_last_exon(self):
        c = _consequence(ConsequenceType.STOP_GAINED, exon="24/24")
        assert predicts_nmd(c) is False

    def test_single_exon(self):
        c = _consequence(ConsequenceType.STOP_GAINED, exon="1/1")
        assert predicts_nmd(c) is False

    def test_is_last_exon(self):
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="12/12")
        assert is_last_exon(c) is True

    def test_is_penultimate_exon(self):
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="11/12")
        assert is_penultimate_exon(c) is True


class TestGeneLoFMechanism:
    def test_lof_intolerant_low_loeuf(self):
        assert gene_has_lof_mechanism(None, gnomad_loeuf=0.10) is True

    def test_lof_tolerant_high_loeuf(self):
        assert gene_has_lof_mechanism(None, gnomad_loeuf=0.80) is False

    def test_no_signals_defaults_to_false(self):
        # No ClinVar P/LP null variants AND no LOEUF → LoF NOT established.
        assert gene_has_lof_mechanism(None, gnomad_loeuf=None) is False


class TestPVS1DecisionTree:
    def setup_method(self):
        from unittest.mock import MagicMock
        self.cfg = MagicMock()

    def test_frameshift_nmd_no_rescue_very_strong(self):
        from unittest.mock import patch
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="5/24")
        gd = GnomADData(loeuf=0.10)
        ann = AnnotationData(consequences=[c], gnomad=gd)
        v = VariantRecord(chrom="chr17", pos=100, ref="G", alt="GA", assembly=Assembly.GRCH38)
        # Bypass the ClinGen SVI strength cap by simulating a gene with
        # established LoF mechanism (>= _MIN_PLP_NULL_FOR_FULL_PVS1 P/LP nulls).
        with patch(
            "acmg_classifier.local_db.clinvar_sqlite.query_pathogenic_null_count",
            return_value=5,
        ), patch(
            "acmg_classifier.local_db.clinvar_sqlite.query_pathogenic_missense_count",
            return_value=0,
        ):
            strength, evidence = evaluate_pvs1(v, ann, self.cfg)
        assert strength == CriterionStrength.VERY_STRONG

    def test_start_loss_moderate(self):
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        c = _consequence(ConsequenceType.START_LOST, exon="1/24")
        ann = AnnotationData(consequences=[c])
        v = VariantRecord(chrom="chr17", pos=100, ref="G", alt="A", assembly=Assembly.GRCH38)
        strength, evidence = evaluate_pvs1(v, ann, self.cfg)
        assert strength == CriterionStrength.MODERATE

    def test_last_exon_frameshift_no_domain_not_met(self):
        # SVI: a last-exon truncation escapes NMD; without evidence that a
        # critical region is removed (no functional domain), PVS1 is N/A. The
        # old behaviour over-applied Moderate (an eRepo PVS1 false-positive
        # source on APC/MYOC-style last-exon truncations).
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="24/24")
        gd = GnomADData(loeuf=0.10)
        ann = AnnotationData(consequences=[c], gnomad=gd)
        v = VariantRecord(chrom="chr17", pos=100, ref="GA", alt="G", assembly=Assembly.GRCH38)
        strength, evidence = evaluate_pvs1(v, ann, self.cfg)
        assert strength == CriterionStrength.NOT_MET

    def test_last_exon_frameshift_with_domain_strong(self):
        # A functional domain in the truncated last-exon region → still PVS1_Strong.
        from unittest.mock import patch
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="24/24")
        c = c.model_copy(update={"domains": ["Pfam:PF00001"]})
        gd = GnomADData(loeuf=0.10)
        ann = AnnotationData(consequences=[c], gnomad=gd)
        v = VariantRecord(chrom="chr17", pos=100, ref="GA", alt="G", assembly=Assembly.GRCH38)
        # Established LoF mechanism (>= P/LP nulls) so the SVI cap does not apply.
        with patch(
            "acmg_classifier.local_db.clinvar_sqlite.query_pathogenic_null_count",
            return_value=5,
        ), patch(
            "acmg_classifier.local_db.clinvar_sqlite.query_pathogenic_missense_count",
            return_value=0,
        ):
            strength, evidence = evaluate_pvs1(v, ann, self.cfg)
        assert strength == CriterionStrength.STRONG
