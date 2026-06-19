"""PM1 manually-curated hotspots (build_pm1_hotspots._CURATED) and committed TSV.

Specs whose residue/region lists the free-text miner cannot reliably extract
(OTC's 21-residue table, GALT's range written as endpoints, the KCNQ4/PDHA1
multi-gene-panel hotspots, KCNQ1's pore helix, JAK3's JH2 residues) are pinned
via a curation override. These tests assert the override is applied in the build
and lands correctly in the committed TSV / loader.
"""
import csv
import importlib.util
from pathlib import Path

from acmg_classifier.criteria.pm1_hotspots import PM1Hotspots
from acmg_classifier.models.enums import CriterionStrength

_ROOT = Path(__file__).resolve().parents[2]
_BUILD = _ROOT / "scripts" / "build_pm1_hotspots.py"
_spec = importlib.util.spec_from_file_location("build_pm1_hotspots", _BUILD)
bp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bp)

_TSV = _ROOT / "resources" / "shared" / "pm1_hotspots.tsv"


def _committed():
    with _TSV.open(encoding="utf-8") as f:
        return {(r["gene_symbol"], r["strength"]): (r["regions"], r["residues"])
                for r in csv.DictReader(f, delimiter="\t")}


class TestCuratedConstant:
    def test_galt_is_range_not_endpoints(self):
        regions, residues = bp._CURATED[("GALT", "Moderate")]
        assert regions == [(171, 188)] and residues == []

    def test_otc_full_residue_set(self):
        _, residues = bp._CURATED[("OTC", "Moderate")]
        assert len(residues) == 21
        assert {90, 268, 330} <= set(residues)   # endpoints of the table

    def test_pdha1_residue_count(self):
        _, residues = bp._CURATED[("PDHA1", "Moderate")]
        assert len(residues) == 59
        assert {88, 118, 316} <= set(residues)

    def test_jak3_and_kcnq(self):
        assert bp._CURATED[("JAK3", "Moderate")] == ([], [651, 759])
        assert bp._CURATED[("KCNQ4", "Moderate")] == ([(271, 292)], [])
        assert bp._CURATED[("KCNQ1", "Moderate")] == ([(300, 320)], [])


class TestCommittedTSV:
    def test_curated_rows_present(self):
        com = _committed()
        assert com[("GALT", "Moderate")] == ("171-188", "")
        assert com[("KCNQ4", "Moderate")] == ("271-292", "")
        assert com[("KCNQ1", "Moderate")] == ("300-320", "")
        assert com[("JAK3", "Moderate")] == ("", "651,759")
        # OTC: full set replaced the single mined residue 268.
        assert com[("OTC", "Moderate")][1].startswith("90,91,92,93,117")
        assert "268" in com[("OTC", "Moderate")][1].split(",")


class TestLoaderOnCommitted:
    def _h(self):
        return PM1Hotspots(_TSV)

    def test_galt_range_interior_awarded(self):
        h = self._h()
        # Interior residue 180 (previously unawarded when stored as endpoints).
        assert h.lookup("GALT", 180) == CriterionStrength.MODERATE
        assert h.lookup("GALT", 171) == CriterionStrength.MODERATE
        assert h.lookup("GALT", 200) is None

    def test_otc_extra_residue_awarded(self):
        h = self._h()
        assert h.lookup("OTC", 92) == CriterionStrength.MODERATE   # newly added
        assert h.lookup("OTC", 268) == CriterionStrength.MODERATE  # originally only one
        assert h.lookup("OTC", 500) is None

    def test_new_genes_lookup(self):
        h = self._h()
        assert h.lookup("KCNQ4", 280) == CriterionStrength.MODERATE
        assert h.lookup("KCNQ1", 310) == CriterionStrength.MODERATE
        assert h.lookup("JAK3", 651) == CriterionStrength.MODERATE
        assert h.lookup("PDHA1", 292) == CriterionStrength.MODERATE
        assert h.lookup("KCNQ4", 300) is None
