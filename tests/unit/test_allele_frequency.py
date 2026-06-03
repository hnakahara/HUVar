"""Unit tests for Whiffin/Ware disease-specific AF thresholds."""
import math
from pathlib import Path

from acmg_classifier.criteria.allele_frequency import (
    DiseaseThresholds,
    GeneThresholds,
    compute_max_credible_af,
)


class TestComputeMaxCredibleAF:
    def test_dominant_formula(self):
        # G = 0.0005 * 0.5 * 0.5 / 1.0 = 1.25e-4 ; AD → G/2 = 6.25e-5
        af = compute_max_credible_af(0.0005, 0.5, 0.5, 1.0, recessive=False)
        assert af == 6.25e-05

    def test_recessive_formula(self):
        # Same G = 1.25e-4 ; AR → sqrt(G)
        af = compute_max_credible_af(0.0005, 0.5, 0.5, 1.0, recessive=True)
        assert math.isclose(af, math.sqrt(1.25e-4))

    def test_reduced_penetrance_raises_threshold(self):
        full = compute_max_credible_af(0.0005, 1.0, 1.0, 1.0, recessive=False)
        half = compute_max_credible_af(0.0005, 1.0, 1.0, 0.5, recessive=False)
        assert half == 2 * full  # lower penetrance → more carriers tolerated

    def test_missing_param_returns_none(self):
        assert compute_max_credible_af(None, 1.0, 1.0, 1.0, recessive=False) is None
        assert compute_max_credible_af(0.0005, 1.0, 1.0, None, recessive=False) is None

    def test_out_of_range_returns_none(self):
        assert compute_max_credible_af(0.0005, 1.0, 1.0, 0.0, recessive=False) is None  # penetrance 0
        assert compute_max_credible_af(0.0005, 1.5, 1.0, 1.0, recessive=False) is None  # het > 1
        assert compute_max_credible_af(1.0, 1.0, 1.0, 1.0, recessive=False) is None     # prevalence not <1


def _write_tsv(tmp_path: Path, rows: str) -> Path:
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(rows, encoding="utf-8")
    return p


class TestDiseaseThresholds:
    def test_computes_from_params_with_floor_and_cap(self, tmp_path: Path):
        # maxAF (AD) = 6.25e-5 → bs1 = max(6.25e-5, 5e-4) = 5e-4 (floor),
        #                        ba1 = min(0.05, 10*6.25e-5 = 6.25e-4) = 6.25e-4
        tsv = _write_tsv(
            tmp_path,
            "gene_symbol\tinheritance\tprevalence\tallelic_het\tgenetic_het\tpenetrance\n"
            "GENEA\tAD\t0.0005\t0.5\t0.5\t1.0\n",
        )
        dt = DiseaseThresholds(tsv)
        t = dt.get("GENEA")
        assert t.bs1 == 0.0005           # floored
        assert t.ba1 == 6.25e-4

    def test_direct_override_takes_precedence(self, tmp_path: Path):
        tsv = _write_tsv(
            tmp_path,
            "gene_symbol\tbs1_threshold\tba1_threshold\n"
            "GENEB\t0.01\t0.02\n",
        )
        t = DiseaseThresholds(tsv).get("GENEB")
        assert t.bs1 == 0.01
        assert t.ba1 == 0.02

    def test_recessive_gene_uses_sqrt(self, tmp_path: Path):
        tsv = _write_tsv(
            tmp_path,
            "gene_symbol\tinheritance\tprevalence\tallelic_het\tgenetic_het\tpenetrance\n"
            "GENEC\tAR\t0.0005\t1.0\t1.0\t1.0\n",
        )
        t = DiseaseThresholds(tsv).get("GENEC")
        # maxAF = sqrt(0.0005) ≈ 0.02236 → bs1 above floor, ba1 capped at 0.05
        assert math.isclose(t.bs1, math.sqrt(0.0005))
        assert t.ba1 == 0.05            # 10*0.02236 > 0.05 → capped

    def test_ba1_cap_applied(self, tmp_path: Path):
        tsv = _write_tsv(
            tmp_path,
            "gene_symbol\tinheritance\tprevalence\tallelic_het\tgenetic_het\tpenetrance\n"
            "GENED\tAD\t0.01\t1.0\t1.0\t1.0\n",
        )
        # maxAF = 0.005 → 10x = 0.05 → cap exactly 0.05
        t = DiseaseThresholds(tsv).get("GENED")
        assert t.ba1 == 0.05

    def test_unknown_gene_uses_flat_defaults(self, tmp_path: Path):
        tsv = _write_tsv(tmp_path, "gene_symbol\tprevalence\n")
        t = DiseaseThresholds(tsv).get("NOPE")
        assert t == GeneThresholds(ba1=0.05, bs1=0.005)

    def test_missing_file_yields_defaults(self, tmp_path: Path):
        dt = DiseaseThresholds(tmp_path / "absent.tsv")
        assert dt.get("ANY") == GeneThresholds(ba1=0.05, bs1=0.005)

    def test_incomplete_params_fall_back_to_default(self, tmp_path: Path):
        # prevalence present but penetrance missing → cannot compute → defaults
        tsv = _write_tsv(
            tmp_path,
            "gene_symbol\tinheritance\tprevalence\n"
            "GENEE\tAD\t0.0005\n",
        )
        t = DiseaseThresholds(tsv).get("GENEE")
        assert t == GeneThresholds(ba1=0.05, bs1=0.005)

    def test_af_basis_males_read_from_column(self, tmp_path: Path):
        tsv = _write_tsv(
            tmp_path,
            "gene_symbol\tbs1_threshold\tba1_threshold\taf_basis\n"
            "RPGR\t0.000083\t0.05\tmales\n"
            "MECP2\t0.0000083\t0.000083\t\n",
        )
        dt = DiseaseThresholds(tsv)
        assert dt.get("RPGR").af_basis == "males"
        assert dt.get("MECP2").af_basis == ""        # blank → overall population
        assert dt.get("UNKNOWN").af_basis == ""      # default thresholds

    def test_het_defaults_to_one_when_omitted(self, tmp_path: Path):
        # Only prevalence + penetrance given; het defaults to 1.0.
        # maxAF (AD) = 0.001*1*1/1 /2 = 5e-4 → bs1 = max(5e-4,5e-4)=5e-4
        tsv = _write_tsv(
            tmp_path,
            "gene_symbol\tinheritance\tprevalence\tpenetrance\n"
            "GENEF\tAD\t0.001\t1.0\n",
        )
        t = DiseaseThresholds(tsv).get("GENEF")
        assert math.isclose(t.bs1, 5e-4)
        assert math.isclose(t.ba1, min(0.05, 10 * 5e-4))
