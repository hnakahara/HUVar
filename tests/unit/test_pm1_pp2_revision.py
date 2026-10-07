"""PM1 / PP2 heuristic revisions (Genome Medicine revision, Reviewer 3).

PM1 fallback: count only P/LP *missense* changes, never the query variant itself,
with a configurable window / minimum. PP2 fallback: gnomAD common missense
variants (precomputed per gene) count on the benign side; thresholds come from
Config.
"""
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.local_db.clinvar_sqlite import query_hotspot_cluster, query_pp2_eligible

_SCHEMA = """
CREATE TABLE variants (
    variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER
);
"""


def _db(tmp_path: Path, rows) -> Path:
    p = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, hgvs_p, "
        "amino_acid_change, codon_position, clinical_significance, review_status, "
        "star_rating) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit()
    con.close()
    return p


def _r(vid, pos, hgvs_p, codon, sig="Pathogenic", gene="G", stars=1):
    return (vid, "1", pos, "A", "G", gene, f"NM_1:{hgvs_p}", hgvs_p[2:], codon, sig,
            "criteria provided", stars)


class TestPM1Heuristic:
    def test_three_missense_is_hotspot(self, tmp_path):
        db = _db(tmp_path, [_r("1", 10, "p.Arg100His", 100), _r("2", 11, "p.Gly105Asp", 105),
                            _r("3", 12, "p.Leu110Pro", 110)])
        ok, ev = query_hotspot_cluster(db, "G", 102)
        assert ok and "3 other P/LP missense" in ev

    def test_truncating_and_synonymous_not_counted(self, tmp_path):
        db = _db(tmp_path, [_r("1", 10, "p.Arg100His", 100), _r("2", 11, "p.Gly105Ter", 105),
                            _r("3", 12, "p.Leu110ProfsTer4", 110), _r("4", 13, "p.Leu111=", 111)])
        ok, _ = query_hotspot_cluster(db, "G", 102)
        assert not ok

    def test_query_variant_excluded(self, tmp_path):
        db = _db(tmp_path, [_r("1", 10, "p.Arg100His", 100), _r("2", 11, "p.Gly105Asp", 105),
                            _r("3", 12, "p.Leu110Pro", 110)])
        ok, _ = query_hotspot_cluster(db, "G", 110, chrom="chr1", pos=12, ref="A", alt="G")
        assert not ok  # only two *other* variants remain

    def test_benign_missense_blocks(self, tmp_path):
        db = _db(tmp_path, [_r("1", 10, "p.Arg100His", 100), _r("2", 11, "p.Gly105Asp", 105),
                            _r("3", 12, "p.Leu110Pro", 110),
                            _r("4", 13, "p.Ala108Thr", 108, sig="Likely benign")])
        ok, ev = query_hotspot_cluster(db, "G", 102)
        assert not ok and "Benign missense" in ev

    def test_configurable_window_and_minimum(self, tmp_path):
        db = _db(tmp_path, [_r("1", 10, "p.Arg100His", 100), _r("2", 11, "p.Gly105Asp", 105),
                            _r("3", 12, "p.Leu130Pro", 130)])
        assert not query_hotspot_cluster(db, "G", 102, window=10)[0]
        assert query_hotspot_cluster(db, "G", 102, window=10, min_path_variants=2)[0]


def _pp2_rows(gene, n_path):
    return [_r(f"P{i}", 100 + i, f"p.Val{100 + i}Ile", 100 + i, gene=gene) for i in range(n_path)]


class TestPP2CommonMissense:
    def test_gnomad_common_missense_counts_as_benign(self, tmp_path):
        db = _db(tmp_path, _pp2_rows("GC", 20))
        assert query_pp2_eligible(db, "GC")[0]
        ok, ev = query_pp2_eligible(db, "GC", common_missense=5)
        assert not ok and "gnomAD common 5" in ev

    def test_thresholds_are_parameters(self, tmp_path):
        db = _db(tmp_path, _pp2_rows("GT", 6))
        assert not query_pp2_eligible(db, "GT")[0]
        assert query_pp2_eligible(db, "GT", min_path=5)[0]


class TestPP2Evaluator:
    def test_evaluator_reads_gene_stats(self, tmp_path, monkeypatch):
        from acmg_classifier.criteria.pathogenic import pp2 as pp2_mod
        from acmg_classifier.models.enums import ConsequenceType

        db = _db(tmp_path, _pp2_rows("GE", 20))
        stats = tmp_path / "pp2_gene_stats.tsv"
        stats.write_text("gene\tcommon_missense\nGE\t4\n", encoding="utf-8")
        pp2_mod.load_common_missense.cache_clear()

        cfg = MagicMock()
        cfg.clinvar_sqlite = db
        cfg.pp2_gene_stats_tsv = stats
        monkeypatch.setattr(pp2_mod, "PP2Applicability", lambda _p: MagicMock(get=lambda g: None))
        ev = pp2_mod.PP2Evaluator(cfg)
        ann = MagicMock()
        ann.primary_consequence.consequence = ConsequenceType.MISSENSE
        ann.primary_consequence.gene_symbol = "GE"
        ann.gnomad = None
        res = ev.evaluate(MagicMock(), ann, None)
        assert not res.triggered
        assert "gnomAD common 4" in res.evidence

        stats.unlink()
        pp2_mod.load_common_missense.cache_clear()
        assert ev.evaluate(MagicMock(), ann, None).triggered


class TestVcepPredictorNote:
    def _spec(self):
        from acmg_classifier.criteria.insilico_genes import InSilicoGeneSpec
        return InSilicoGeneSpec(None)

    def _cfg(self, tool, bayesdel=False, cadd=False):
        from acmg_classifier.models.enums import InSilicoTool
        cfg = MagicMock()
        cfg.insilico_tool = InSilicoTool(tool)
        cfg.use_bayesdel = bayesdel
        cfg.use_cadd = cadd
        return cfg

    def test_bayesdel_gene_under_esm1b(self):
        from acmg_classifier.criteria.insilico_genes import vcep_predictor_note
        note = vcep_predictor_note("BRCA1", self._cfg("esm1b"), self._spec(), None, "PP3")
        assert note and "BayesDel" in note and "--with-bayesdel" in note

    def test_no_note_when_vcep_predictor_active(self):
        from acmg_classifier.criteria.insilico_genes import vcep_predictor_note
        assert vcep_predictor_note("BRCA1", self._cfg("revel", bayesdel=True),
                                   self._spec(), None, "PP3") is None

    def test_revel_cutoff_gene_under_esm1b(self):
        from acmg_classifier.criteria.insilico_genes import vcep_predictor_note
        from acmg_classifier.criteria.revel_genes import RevelRule
        from acmg_classifier.models.enums import CriterionStrength
        rs = MagicMock()
        rs.get.return_value = RevelRule(pp3={CriterionStrength.SUPPORTING: 0.7}, bp4={})
        note = vcep_predictor_note("MYH7", self._cfg("esm1b"), self._spec(), rs, "PP3")
        assert note and "REVEL cutoffs" in note
        assert vcep_predictor_note("MYH7", self._cfg("esm1b"), self._spec(), rs, "BP4") is None
        assert vcep_predictor_note("MYH7", self._cfg("revel"), self._spec(), rs, "PP3") is None
