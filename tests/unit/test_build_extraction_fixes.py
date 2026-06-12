"""Regression tests for cspec-extraction bug fixes in build_disease_thresholds.py.

Covers:
  * #9  _af_basis must flag "males" only for an in-males FREQUENCY rule, not for
        a "hemizygotes" COUNT clause (OTC/SLC6A8) or prevalence note (ABCD1).
  * #10 BP3 applicability must not be negated by a "not applicable" sub-clause
        when a positive "can be applied to <region>" clause is present (VHL),
        and "AA14-AA48" residue ranges must be parsed.
"""
import importlib.util
from pathlib import Path

_BUILD = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BUILD)
b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b)


def _rs(label: str, desc: str) -> dict:
    """A minimal rule set with one applicable criterion code carrying *desc*."""
    return {
        "criteriaCodes": [{
            "label": label,
            "evidenceStrengths": [
                {"label": "Strong", "applicability": "Applicable", "description": desc},
            ],
        }],
    }


class TestAfBasisMales:
    def test_in_males_frequency_is_males(self):
        rs = _rs("BS1", "Allele frequency in males is greater than 0.0001.")
        assert b._af_basis(rs) == "males"

    def test_hemizygote_count_is_not_males(self):
        # OTC-style: "(male) hemizygotes" is a COUNT rule, not a frequency basis.
        rs = _rs("BA1", "AF above 1.0% Grpmax FAF OR >=10 (female) homozygotes "
                        "or (male) hemizygotes in gnomAD.")
        assert b._af_basis(rs) == ""

    def test_prevalence_in_hemizygotes_is_not_males(self):
        # ABCD1-style: "Total Grpmax FAF" with a prevalence note mentioning
        # hemizygotes — overall FAF, not males.
        rs = _rs("BA1", "Use a Total Grpmax FAF cutoff of >=0.00017 "
                        "(prevalence in hemizygotes is 1 in 5000).")
        assert b._af_basis(rs) == ""


class TestBP3VHL:
    def test_positive_clause_keeps_applicable(self):
        rs = _rs("BP3", "BP3 can be applied to the GXEEX repeat (AA14-AA48). "
                        "Otherwise the rest of the gene has no repeats and BP3 "
                        "is not applicable.")
        assert b._bp3_applicability(rs) == "applicable"

    def test_aa_prefixed_range_parsed(self):
        assert b._bp_residue_ranges("repeat motif (AA14-AA48)") == "14-48"

    def test_plain_not_applicable_still_negated(self):
        rs = _rs("BP3", "BP3 is not applicable for this gene.")
        assert b._bp3_applicability(rs) == "not_applicable"


class TestBS1Strength:
    def test_very_strong_tier_when_no_strong(self):
        # GN023-style: BS1 applicable only at Very Strong + Supporting → the
        # chosen tier (first applicable) is Very Strong, not the Strong default.
        rs = {"criteriaCodes": [{
            "label": "BS1",
            "evidenceStrengths": [
                {"label": "Very Strong", "applicability": "Applicable",
                 "description": "MAF of >=0.003 (0.3%) for autosomal recessive."},
                {"label": "Supporting", "applicability": "Applicable",
                 "description": "MAF of >=0.0007 (0.07%) for autosomal recessive."},
            ],
        }]}
        assert b._bs1_strength(rs) == "VeryStrong"

    def test_strong_tier_default(self):
        rs = {"criteriaCodes": [{
            "label": "BS1",
            "evidenceStrengths": [
                {"label": "Strong", "applicability": "Applicable",
                 "description": "AF >= 0.001"},
            ],
        }]}
        assert b._bs1_strength(rs) == "Strong"

    def test_evaluator_emits_tier_strength(self, tmp_path):
        from acmg_classifier.criteria.allele_frequency import DiseaseThresholds
        from acmg_classifier.models.enums import CriterionStrength
        p = tmp_path / "disease_prevalence.tsv"
        p.write_text(
            "gene_symbol\tbs1_threshold\tbs1_strength\n"
            "MYO15A\t0.003\tVeryStrong\n",
            encoding="utf-8",
        )
        gt = DiseaseThresholds(p).get("MYO15A")
        assert gt.bs1_strength == CriterionStrength.VERY_STRONG


class TestBS1ParenPercent:
    """#: a parenthesised, operator-prefixed percentage cutoff
    "(>0.0185%)" must parse to its proportion (PIK3-pathway VCEP GN018 BS1)."""

    def test_paren_gt_percent(self):
        thr, _ = b._threshold_from_desc(
            "Allele frequency (>0.0185%). An allele frequency (>0.0185%) "
            "was approved. (Supplemental Table 3)."
        )
        assert thr == 0.000185

    def test_plain_gt_percent(self):
        thr, _ = b._threshold_from_desc("Allele frequency > 0.0185%.")
        assert thr == 0.000185


class TestEmptyPlaceholderSpec:
    """An empty Pilot/In-Prep spec (0 criteriaCodes) must not shadow a populated
    spec's BS1/BA1, even when it is more gene-specific. Verified on the committed
    TSV: these genes each have an empty single-gene pilot spec AND a populated
    grouped Released spec; the grouped thresholds must survive."""

    def _rows(self):
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            return {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}

    def test_pik3_pathway_bs1_ba1_survive(self):
        rows = self._rows()
        for gene in ("AKT3", "MTOR", "PIK3CA", "PIK3R2"):
            assert rows[gene]["bs1_threshold"] == "0.000185", gene
            assert rows[gene]["ba1_threshold"] == "0.000926", gene

    def test_hearing_loss_panel_survives(self):
        rows = self._rows()
        for gene in ("CDH23", "COCH", "GJB2", "KCNQ4", "MYO6", "MYO7A",
                     "SLC26A4", "TECTA", "USH2A", "MYO15A", "OTOF"):
            assert rows[gene]["bs1_threshold"] == "0.003", gene
            assert rows[gene]["ba1_threshold"] == "0.005", gene

    def test_mito_panel_survives(self):
        rows = self._rows()
        assert rows["POLG"]["bs1_threshold"] == "0.005"
        assert rows["ETHE1"]["bs1_threshold"] == "0.0002"
        assert rows["PDHA1"]["bs1_threshold"] == "0.000092"
        assert rows["SLC19A3"]["bs1_threshold"] == "0.0005"


class TestBS1VariantExclusion:
    """MYOC BS1 "Does not apply to p.Gln368Ter" — a recurrent disease allele
    whose population frequency must not be read as benign. Curated into the
    bs1_exclude column (absent from the JSON-LD BS1 description)."""

    def test_committed_tsv_has_exclusion(self):
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        assert rows["MYOC"]["bs1_exclude"] == "p.Gln368Ter"

    def test_resolve_row_loads_exclude(self):
        from acmg_classifier.criteria.allele_frequency import _resolve_row
        gt = _resolve_row({"gene_symbol": "MYOC", "bs1_threshold": "0.001",
                           "bs1_exclude": "p.Gln368Ter"})
        assert gt.bs1_exclude == "p.Gln368Ter"

    def _ev(self, tmp_path):
        from unittest.mock import MagicMock
        from acmg_classifier.criteria.benign.bs1 import BS1Evaluator
        p = tmp_path / "disease_prevalence.tsv"
        p.write_text(
            "gene_symbol\tbs1_threshold\tbs1_exclude\n"
            "MYOC\t0.001\tp.Gln368Ter\n",
            encoding="utf-8",
        )
        cfg = MagicMock()
        cfg.disease_prevalence_tsv = p
        return BS1Evaluator(cfg)

    def _ann(self, hgvs_p):
        from acmg_classifier.models.annotation import (
            AnnotationData, GnomADData, ConsequenceInfo,
        )
        from acmg_classifier.models.enums import ConsequenceType
        return AnnotationData(
            gnomad=GnomADData(faf95_popmax=0.02, filter_pass=True),  # >> 0.001
            consequences=[ConsequenceInfo(
                transcript_id="NM_x", gene_id="ENSG", gene_symbol="MYOC",
                consequence=ConsequenceType.STOP_GAINED, biotype="protein_coding",
                hgvs_p=hgvs_p,
            )],
        )

    def test_excluded_variant_not_met(self, tmp_path):
        from acmg_classifier.models.variant import VariantRecord
        from acmg_classifier.models.enums import Assembly
        ev = self._ev(tmp_path)
        v = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)
        r = ev.evaluate(v, self._ann("NP_000252.1:p.Gln368Ter"))
        assert not r.triggered
        assert "excludes BS1" in r.evidence

    def test_star_notation_also_excluded(self, tmp_path):
        from acmg_classifier.models.variant import VariantRecord
        from acmg_classifier.models.enums import Assembly
        ev = self._ev(tmp_path)
        v = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)
        r = ev.evaluate(v, self._ann("NP_000252.1:p.Gln368*"))
        assert not r.triggered

    def test_other_variant_still_fires(self, tmp_path):
        from acmg_classifier.models.variant import VariantRecord
        from acmg_classifier.models.enums import Assembly
        ev = self._ev(tmp_path)
        v = VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)
        # A different high-frequency MYOC variant is not excluded → BS1 fires.
        r = ev.evaluate(v, self._ann("NP_000252.1:p.Gln368His"))
        assert r.triggered


class TestRPGRThresholdCorrection:
    """RPGR's cspec carries typo'd AF thresholds (BS1 8.3e-5; a legacy '5%' BA1
    boilerplate → 0.05). The VCEP's published classifications use BS1 > 5e-6 and
    BA1 > 5e-5 (BA1 = 10 x BS1), applied via the curated override."""

    def _rows(self):
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            return {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}

    def test_corrected_thresholds(self):
        r = self._rows()["RPGR"]
        assert r["bs1_threshold"] == "0.000005"
        assert r["ba1_threshold"] == "0.00005"
        # BA1 = 10 x BS1, and the male AF basis is preserved.
        assert r["af_basis"] == "males"

    def test_ba1_above_bs1(self):
        # Sanity: the correction restores the normal BA1 >= BS1 ordering.
        r = self._rows()["RPGR"]
        assert float(r["ba1_threshold"]) > float(r["bs1_threshold"])


def test_bp1_truncating_includes_lof_classes():
    """#13: RASopathy GoF BP1 truncating target covers splice / start-loss /
    whole-gene deletion, not only nonsense/frameshift."""
    from acmg_classifier.criteria.benign.bp1 import _TRUNCATING
    from acmg_classifier.models.enums import ConsequenceType
    assert ConsequenceType.SPLICE_ACCEPTOR in _TRUNCATING
    assert ConsequenceType.SPLICE_DONOR in _TRUNCATING
    assert ConsequenceType.START_LOST in _TRUNCATING
    assert ConsequenceType.TRANSCRIPT_ABLATION in _TRUNCATING
