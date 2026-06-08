"""gnomAD multi-row merge (GRCh37 exomes+genomes "either meets" semantics)."""
from acmg_classifier.local_db.gnomad_db import _merge_rows, _pass_filter


# Row layout: (af, an, ac, nhomalt, nhemi, popmax_af, popmax_pop, faf95_popmax,
#              af_xy, ac_xx, nhomalt_xx, filters)
def _row(af, popmax_af, faf95, nhomalt=0, pop="afr", filters="PASS", af_xy=None,
         ac_xx=None, nhomalt_xx=None):
    return (af, 100000, int(af * 100000), nhomalt, 0, popmax_af, pop, faf95,
            af_xy, ac_xx, nhomalt_xx, filters)


class TestPassFilter:
    def test_pass_variants(self):
        assert _pass_filter(None) and _pass_filter("") and _pass_filter("PASS") and _pass_filter(".")

    def test_fail(self):
        assert not _pass_filter("AC0") and not _pass_filter("InbreedingCoeff")


class TestMergeRows:
    def test_single_row_passthrough(self):
        gd = _merge_rows([_row(0.001, 0.002, 0.0015)])
        assert gd.af == 0.001 and gd.popmax_af == 0.002 and gd.filter_pass

    def test_per_field_max_either_meets(self):
        # exomes: higher FAF; genomes: higher homozygotes -> each field maxed.
        exo = _row(0.0003, 0.00032, 0.00022, nhomalt=2, pop="afr")
        gen = _row(0.0001, 0.00015, 0.00012, nhomalt=9, pop="nfe")
        gd = _merge_rows([exo, gen])
        assert gd.popmax_af == 0.00032        # higher AF wins (BA1/BS1)
        assert gd.faf95_popmax == 0.00022
        assert gd.nhomalt == 9                 # higher homozygote count wins (BS2)
        assert gd.popmax_pop == "afr"          # pop of the higher-AF row
        assert gd.filter_pass

    def test_filtered_row_excluded_from_merge(self):
        # The high-AF row is FILTERED; only the PASS row contributes the frequency.
        filtered = _row(0.02, 0.03, 0.025, filters="AC0")
        passed = _row(0.0001, 0.00012, 0.0001, filters="PASS")
        gd = _merge_rows([filtered, passed])
        assert gd.popmax_af == 0.00012 and gd.filter_pass

    def test_all_filtered_reports_filter_failed(self):
        gd = _merge_rows([_row(0.01, 0.02, 0.015, filters="AC0")])
        assert not gd.filter_pass

    def test_female_counts_merged_by_max(self):
        # ac_xx / nhomalt_xx (female counts for BS2) are merged by per-field MAX
        # like the other counts; None contributes nothing.
        a = _row(0.0003, 0.0003, 0.0002, ac_xx=5, nhomalt_xx=0)
        b = _row(0.0001, 0.0001, 0.0001, ac_xx=None, nhomalt_xx=1)
        gd = _merge_rows([a, b])
        assert gd.ac_xx == 5 and gd.nhomalt_xx == 1
