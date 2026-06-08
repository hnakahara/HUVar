"""GrpMax FAF extraction across gnomAD versions (v4 direct vs v2.1.1 per-pop max)."""
from acmg_classifier.setup.gnomad_builder import _faf95_popmax


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
