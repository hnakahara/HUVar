"""GrpMax FAF extraction across gnomAD versions (v4 direct vs v2.1.1 per-pop max)."""
import math

import pytest

from acmg_classifier.setup.gnomad_builder import (
    _faf95_from_counts,
    _faf95_noncancer_popmax,
    _faf95_popmax,
    _nhemi,
)


class _FakeINFO:
    def __init__(self, d):
        self._d = d

    def get(self, key):
        return self._d.get(key)


class _FakeVariant:
    def __init__(self, info):
        self.INFO = _FakeINFO(info)


def test_v4_uses_direct_joint_field():
    v = _FakeVariant({"fafmax_faf95_max_joint": 0.0012, "faf95_nfe": 0.5})
    # The direct joint field wins over any per-pop value.
    assert _faf95_popmax(v) == 0.0012


def test_v4_falls_back_to_exome_fafmax():
    v = _FakeVariant({"fafmax_faf95_max": 0.0009})
    assert _faf95_popmax(v) == 0.0009


def test_v2_computes_max_over_continental_pops():
    # v2.1.1 has no popmax-FAF field — take the max over faf95_<pop>.
    v = _FakeVariant({
        "faf95_afr": 1.0e-5,
        "faf95_amr": 8.0e-5,
        "faf95_eas": 0.0,
        "faf95_nfe": 3.0e-5,
        "faf95_sas": 2.0e-5,
    })
    assert _faf95_popmax(v) == 8.0e-5


def test_v2_zero_is_kept_not_dropped():
    # A genuine FAF95 of 0.0 (wide CI on a sparse variant) must be returned as
    # 0.0, not treated as "missing" — that is exactly the over-firing the
    # per-pop computation prevents.
    v = _FakeVariant({"faf95_afr": 0.0, "faf95_nfe": 0.0})
    assert _faf95_popmax(v) == 0.0


def test_no_faf_fields_returns_none():
    v = _FakeVariant({"AF": 0.01})
    assert _faf95_popmax(v) is None


class TestFaf95FromCounts:
    """Poisson 95% CI filtering AF recomputed from a raw (AC, AN) pair — needed
    for the non-cancer subset, whose faf95 gnomAD does NOT precompute."""

    def test_singleton_matches_poisson_closed_form(self):
        # For AC=1 the Poisson 95% lower bound has a closed form: the regularised
        # lower incomplete gamma P(1, x) = 1 - e^-x, so x = -ln(0.95). This is the
        # textbook gnomAD/Whiffin result and anchors the method to the published
        # definition (a full cross-check against the VCF's precomputed faf95_<pop>
        # can be run on the server after the rebuild).
        pytest.importorskip("scipy")
        expected = -math.log(0.95) / 100_000
        assert _faf95_from_counts(1, 100_000) == pytest.approx(expected, rel=1e-6)

    def test_zero_ac_is_zero(self):
        # No observed alt alleles → the lower confidence bound is exactly 0.
        assert _faf95_from_counts(0, 50_000) == 0.0

    def test_missing_or_zero_an_is_none(self):
        # No callable denominator → None (no FAF), never a divide-by-zero.
        assert _faf95_from_counts(5, 0) is None
        assert _faf95_from_counts(5, None) is None

    def test_missing_ac_with_valid_an_is_zero(self):
        # AN present but AC absent = no alt alleles observed in that group → the
        # lower confidence bound is 0.0 (harmless in the popmax max()).
        assert _faf95_from_counts(None, 10_000) == 0.0

    def test_faf_below_point_estimate(self):
        # The 95% lower bound must sit below the point AF = AC/AN.
        pytest.importorskip("scipy")
        faf = _faf95_from_counts(50, 100_000)
        assert faf is not None and 0.0 < faf < 50 / 100_000


class TestFaf95NoncancerPopmax:
    """Non-cancer popmax FAF = max over the continental groups of the per-group
    Poisson FAF computed from AC_non_cancer_<pop>/AN_non_cancer_<pop>."""

    def test_takes_max_over_continental_groups(self):
        pytest.importorskip("scipy")
        # nfe carries the highest AC/AN ratio → it sets the popmax FAF. asj/fin are
        # founder groups and are intentionally NOT in the popmax basis.
        v = _FakeVariant({
            "AC_non_cancer_afr": 1, "AN_non_cancer_afr": 100_000,
            "AC_non_cancer_nfe": 40, "AN_non_cancer_nfe": 100_000,
            "AC_non_cancer_asj": 500, "AN_non_cancer_asj": 10_000,  # excluded
        })
        expected = _faf95_from_counts(40, 100_000)
        assert _faf95_noncancer_popmax(v) == pytest.approx(expected)

    def test_no_noncancer_counts_returns_none(self):
        # A site with no non-cancer per-group counts has no non-cancer FAF.
        v = _FakeVariant({"AF": 0.01, "faf95_nfe": 0.5})
        assert _faf95_noncancer_popmax(v) is None


class TestNhemi:
    """Hemizygote count is derived from AC_XY on chrX/chrY (gnomAD has no
    nhemi field); on autosomes it must stay None."""

    def test_x_uses_ac_xy_joint(self):
        v = _FakeVariant({"AC_joint_XY": 7, "AC_XY": 99})
        assert _nhemi(v, "X") == 7

    def test_x_falls_back_to_ac_xy(self):
        v = _FakeVariant({"AC_XY": 5})
        assert _nhemi(v, "X") == 5

    def test_y_uses_ac_xy(self):
        v = _FakeVariant({"AC_XY": 3})
        assert _nhemi(v, "Y") == 3

    def test_v2_male_field(self):
        # gnomAD v2.1.1 names the male allele count AC_male, not AC_XY.
        v = _FakeVariant({"AC_male": 4})
        assert _nhemi(v, "X") == 4

    def test_autosome_returns_none(self):
        # AC_XY on an autosome is a diploid male allele count, NOT hemizygotes.
        v = _FakeVariant({"AC_XY": 12})
        assert _nhemi(v, "1") is None

    def test_genuine_nhemi_field_preferred(self):
        v = _FakeVariant({"nhemi_joint": 9, "AC_XY": 2})
        assert _nhemi(v, "X") == 9

    def test_missing_returns_none(self):
        v = _FakeVariant({"AC": 100})
        assert _nhemi(v, "X") is None
