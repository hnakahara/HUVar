"""PP2 VCEP-applicability: cspec extraction + evaluator gating.

The ``pp2`` column of ``disease_prevalence.tsv`` carries each VCEP's explicit
PP2 decision; it overrides the statistical heuristic to curb over-assignment.
"""
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pp2_genes import PP2Applicability
from acmg_classifier.criteria.pathogenic.pp2 import PP2Evaluator
from acmg_classifier.criteria.registry import CriteriaRegistry
from acmg_classifier.models.annotation import (
    AnnotationData, GnomADData, ConsequenceInfo,
)
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    Assembly, ACMGCriterion, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.models.variant import VariantRecord


# --- load the build script as a module (it lives under scripts/, not the pkg) -
_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


def _pp2_code(applicability, desc):
    return {"label": "PP2", "evidenceStrengths": [
        {"label": "Supporting", "applicability": applicability, "description": desc},
    ]}


class TestPP2Extraction:
    def test_gene_specific_exclusion(self):
        rs = {
            "genes": [{"label": "MTOR"}, {"label": "PIK3R2"}],
            "criteriaCodes": [_pp2_code(
                "Applicable",
                "Award PP2 if z-score > 3.09 (applicable to MTOR but not PIK3R2).",
            )],
        }
        assert bdt._pp2_applicability(rs) == {
            "MTOR": "applicable", "PIK3R2": "not_applicable",
        }

    def test_blanket_negation(self):
        rs = {
            "genes": [{"label": "KCNQ1"}],
            "criteriaCodes": [_pp2_code(
                "Applicable",
                "Not applicable due to presence of benign variation (z-score 1.83).",
            )],
        }
        assert bdt._pp2_applicability(rs) == {"KCNQ1": "not_applicable"}

    def test_declined_when_no_applicable_strength(self):
        rs = {
            "genes": [{"label": "FOO"}],
            "criteriaCodes": [_pp2_code("Not Applicable for this VCEP", "")],
        }
        assert bdt._pp2_applicability(rs) == {"FOO": "not_applicable"}

    def test_plain_applicable(self):
        rs = {
            "genes": [{"label": "PTEN"}],
            "criteriaCodes": [_pp2_code(
                "Applicable", "Missense is a common mechanism; low benign rate.",
            )],
        }
        assert bdt._pp2_applicability(rs) == {"PTEN": "applicable"}

    def test_no_pp2_code_returns_empty(self):
        rs = {"genes": [{"label": "X"}], "criteriaCodes": [{"label": "BA1"}]}
        assert bdt._pp2_applicability(rs) == {}


class TestPP2SpecificityResolution:
    """The most gene-specific spec's PP2 decision wins (single-gene VCEP over a
    grouped panel); on a specificity tie the conservative not_applicable wins."""

    def test_single_gene_not_applicable_beats_grouped_applicable(self):
        # RASopathy: grouped GN004 (12 genes) "applicable" vs single-gene
        # GN039 NRAS "not_applicable" — the gene-specific decision must win.
        grouped = (12, "applicable", "")
        specific = (1, "not_applicable", "")
        assert bdt._pp2_more_specific(specific, grouped)
        assert not bdt._pp2_more_specific(grouped, specific)

    def test_single_gene_applicable_beats_grouped(self):
        grouped = (12, "not_applicable", "")
        specific = (1, "applicable", "PM2,PP3")
        assert bdt._pp2_more_specific(specific, grouped)

    def test_tie_prefers_not_applicable(self):
        # Two single-gene specs disagree (e.g. ACTA1 GN147 vs GN169) — keep the
        # conservative not_applicable.
        na = (1, "not_applicable", "")
        ap = (1, "applicable", "")
        assert bdt._pp2_more_specific(na, ap)
        assert not bdt._pp2_more_specific(ap, na)

    def test_tie_same_decision_is_not_replaced(self):
        assert not bdt._pp2_more_specific((1, "applicable", ""), (1, "applicable", ""))


class TestPP2Requires:
    def test_co_requirement_extracted(self):
        rs = {
            "genes": [{"label": "BMPR2"}],
            "criteriaCodes": [_pp2_code("Applicable", "PM2_supporting and PP3 must be met.")],
        }
        assert bdt._pp2_requires(rs) == "PM2,PP3"

    def test_unconditional_has_no_requirement(self):
        rs = {
            "genes": [{"label": "PTEN"}],
            "criteriaCodes": [_pp2_code("Applicable", "Missense is a common mechanism.")],
        }
        assert bdt._pp2_requires(rs) == ""


def _consequence(gene, ctype=ConsequenceType.MISSENSE):
    return ConsequenceInfo(
        transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
        consequence=ctype, biotype="protein_coding", is_mane_select=True,
        protein_position=100, amino_acid_change="V100I", codon_position=100,
    )


def _cfg(tmp_path):
    tsv = tmp_path / "disease_prevalence.tsv"
    tsv.write_text(
        "gene_symbol\tpp2\n"
        "PTEN\tapplicable\n"
        "KCNQ1\tnot_applicable\n",
        encoding="utf-8",
    )
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = tmp_path / "absent.sqlite"   # forces heuristic to fail
    return cfg


class TestPP2Loader:
    def test_loads_pp2_column(self, tmp_path):
        cfg = _cfg(tmp_path)
        a = PP2Applicability(cfg.disease_prevalence_tsv)
        assert a.get("PTEN") == "applicable"
        assert a.get("KCNQ1") == "not_applicable"
        assert a.get("UNSEEN") == ""

    def test_missing_file_is_empty(self, tmp_path):
        a = PP2Applicability(tmp_path / "nope.tsv")
        assert a.get("PTEN") == ""

    def test_requires_parsed(self, tmp_path):
        tsv = tmp_path / "dp.tsv"
        tsv.write_text(
            "gene_symbol\tpp2\tpp2_requires\n"
            "BMPR2\tapplicable\tPM2,PP3\n"
            "PTEN\tapplicable\t\n",
            encoding="utf-8",
        )
        a = PP2Applicability(tsv)
        assert a.requires("BMPR2") == [ACMGCriterion.PM2, ACMGCriterion.PP3]
        assert a.requires("PTEN") == []
        assert a.requires("UNSEEN") == []


def _bmpr2_tsv(tmp_path):
    tsv = tmp_path / "disease_prevalence.tsv"
    tsv.write_text(
        "gene_symbol\tpp2\tpp2_requires\n"
        "BMPR2\tapplicable\tPM2,PP3\n"
        "PTEN\tapplicable\t\n",
        encoding="utf-8",
    )
    return tsv


def _registry(tmp_path):
    # Bypass full __init__ (which builds every evaluator) — we only exercise the
    # PP2 co-requirement post-hoc pass, which needs just self._pp2.
    reg = object.__new__(CriteriaRegistry)
    reg._pp2 = PP2Applicability(_bmpr2_tsv(tmp_path))
    return reg


def _ann(gene):
    return AnnotationData(gnomad=GnomADData(), consequences=[_consequence(gene)])


class TestRegistryPP2CoRequirement:
    def test_pp2_kept_when_requirements_met(self, tmp_path):
        reg = _registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PP2, CriterionStrength.SUPPORTING),
            CriteriaResult.met(ACMGCriterion.PM2),
            CriteriaResult.met(ACMGCriterion.PP3),
        ]
        reg._apply_pp2_co_requirements(results, _ann("BMPR2"))
        pp2 = next(r for r in results if r.criterion == ACMGCriterion.PP2)
        assert not pp2.suppressed

    def test_pp2_suppressed_when_pp3_missing(self, tmp_path):
        reg = _registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PP2, CriterionStrength.SUPPORTING),
            CriteriaResult.met(ACMGCriterion.PM2),
            CriteriaResult.not_met(ACMGCriterion.PP3, "no in-silico support"),
        ]
        reg._apply_pp2_co_requirements(results, _ann("BMPR2"))
        pp2 = next(r for r in results if r.criterion == ACMGCriterion.PP2)
        assert pp2.suppressed
        assert "requires PM2+PP3" in pp2.evidence

    def test_pp2_suppressed_when_both_missing(self, tmp_path):
        reg = _registry(tmp_path)
        results = [CriteriaResult.met(ACMGCriterion.PP2, CriterionStrength.SUPPORTING)]
        reg._apply_pp2_co_requirements(results, _ann("BMPR2"))
        assert results[0].suppressed

    def test_gene_without_requirement_untouched(self, tmp_path):
        reg = _registry(tmp_path)
        results = [CriteriaResult.met(ACMGCriterion.PP2, CriterionStrength.SUPPORTING)]
        reg._apply_pp2_co_requirements(results, _ann("PTEN"))
        assert not results[0].suppressed


class TestPP2EvaluatorGating:
    def _snv(self):
        return VariantRecord(chrom="chr10", pos=1, ref="G", alt="A", assembly=Assembly.GRCH38)

    def test_vcep_applicable_meets(self, tmp_path):
        ev = PP2Evaluator(_cfg(tmp_path))
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_consequence("PTEN")])
        r = ev.evaluate(self._snv(), ann)
        assert r.triggered
        assert r.strength == CriterionStrength.SUPPORTING
        assert "VCEP designates PP2 applicable" in r.evidence

    def test_vcep_not_applicable_blocks(self, tmp_path):
        ev = PP2Evaluator(_cfg(tmp_path))
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_consequence("KCNQ1")])
        r = ev.evaluate(self._snv(), ann)
        assert not r.triggered
        assert "not applicable" in r.evidence

    def test_unknown_gene_falls_back_to_heuristic(self, tmp_path):
        # No VCEP row + no ClinVar DB → heuristic cannot qualify → not met,
        # but crucially it reached the heuristic (not the VCEP gate).
        ev = PP2Evaluator(_cfg(tmp_path))
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_consequence("NOVCEP")])
        r = ev.evaluate(self._snv(), ann)
        assert not r.triggered
        assert "VCEP" not in r.evidence

    def test_non_missense_not_met(self, tmp_path):
        ev = PP2Evaluator(_cfg(tmp_path))
        ann = AnnotationData(
            gnomad=GnomADData(),
            consequences=[_consequence("PTEN", ConsequenceType.SYNONYMOUS)],
        )
        r = ev.evaluate(self._snv(), ann)
        assert not r.triggered
        assert "Not a missense" in r.evidence

    def test_supplement_overrides_not_applicable(self, tmp_path):
        ev = PP2Evaluator(_cfg(tmp_path))
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_consequence("KCNQ1")])
        supp = [SupplementEntry(
            variant_id="chr10:1:G:A", criterion=ACMGCriterion.PP2,
            strength=CriterionStrength.SUPPORTING, evidence="expert override",
        )]
        r = ev.evaluate(self._snv(), ann, supp)
        assert r.triggered
