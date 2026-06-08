"""BS2 VCEP-applicability: cspec extraction for genes where BS2 cannot be
derived from gnomAD population counts.

Our BS2 evaluator is gnomAD-based, so a VCEP must resolve to ``not_applicable``
when it bars population data or scores BS2 purely on clinical phenotype/points
with no gnomAD-countable rule (RPE65, the RASopathy point specs, the Fanconi
"points per proband" cancer specs, ...).
"""
import importlib.util
from pathlib import Path

# --- load the build script as a module (it lives under scripts/, not the pkg) -
_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


def _bs2_code(desc, applicability="Applicable"):
    return {"label": "BS2", "evidenceStrengths": [
        {"label": "Strong", "applicability": applicability, "description": desc},
    ]}


def _rs(gene, *codes):
    return {"genes": [{"label": gene}], "criteriaCodes": list(codes)}


class TestBS2NotApplicable:
    def test_bars_population_data(self):
        # RASopathy GN004 wording.
        rs = _rs("NRAS", _bs2_code(
            "general population data should not be used for this criterion."))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_gnomad_not_considered(self):
        # RPE65 GN120 wording.
        rs = _rs("RPE65", _bs2_code(
            "Variant is present in >= 3 homozygotes without any features of the "
            "phenotype. Presence in databases such as gnomAD are not considered."))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_not_applicable_due_to(self):
        # KCNQ1 GN112 wording.
        rs = _rs("KCNQ1", _bs2_code("Not applicable due to incomplete penetrance."))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_pure_point_scoring(self):
        # Single-gene RASopathy specs (NRAS/GN039 etc.) inherit point-based,
        # phenotype-driven scoring with no gnomAD-countable rule.
        rs = _rs("NRAS", _bs2_code("\\-4 Points."))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_points_per_proband_phenotype(self):
        # BRCA1/2, PALB2: Fanconi-Anemia "points per proband", no gnomAD count.
        rs = _rs("BRCA1", _bs2_code(
            "Applied in absence of features of recessive disease (Fanconi Anemia "
            "phenotype). Approach to assign points per proband. BS2 = >= 4 points"))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_no_applicable_strength_declined(self):
        rs = _rs("FOO", _bs2_code("", applicability="Not Applicable"))
        assert bdt._bs2_applicability(rs) == "not_applicable"


class TestBS2StaysApplicable:
    def test_points_with_homozygous_gnomad_path(self):
        # APC: points OR a homozygous-state count — gnomAD-derivable, keep on.
        rs = _rs("APC", _bs2_code(
            ">= 10 points for healthy individuals OR >= 2 times in homozygous state."))
        assert bdt._bs2_applicability(rs) == "applicable"

    def test_plain_count_based(self):
        rs = _rs("LDLR", _bs2_code(
            "Variant is identified in >= 3 heterozygous or >= 1 homozygous "
            "well-phenotyped, untreated, normolipidemic adults (unrelated)."))
        assert bdt._bs2_applicability(rs) == "applicable"

    def test_hemizygotes_in_gnomad(self):
        # PDHA1: ">= 16 hemizygotes in gnomAD" — explicitly gnomAD-countable.
        rs = _rs("PDHA1", _bs2_code(
            "Observed in at least two healthy male adults ... AND/OR >= 16 "
            "hemizygotes in gnomAD"))
        assert bdt._bs2_applicability(rs) == "applicable"

    def test_no_bs2_code(self):
        rs = _rs("FOO")
        assert bdt._bs2_applicability(rs) == ""


class TestResolvedTSV:
    """Spot-check the committed TSV reflects the extraction."""

    def test_committed_tsv_marks_target_genes_not_applicable(self):
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        for gene in ("RPE65", "NRAS", "BRAF", "KCNQ1", "BRCA1", "BRCA2", "PALB2"):
            assert rows[gene]["bs2"] == "not_applicable", gene
