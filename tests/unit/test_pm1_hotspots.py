"""PM1 VCEP hotspot regions: cspec mining, loader, and evaluator gating."""
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pm1_hotspots import PM1Hotspots
from acmg_classifier.criteria.pathogenic.pm1 import PM1Evaluator
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_BUILD = Path(__file__).resolve().parents[2] / "scripts" / "build_pm1_hotspots.py"
_spec = importlib.util.spec_from_file_location("build_pm1_hotspots", _BUILD)
bp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bp)


class TestParseRegions:
    def test_codon_range(self):
        ranges, res = bp.parse_regions("Applicable to missense variants in codons 167-931.")
        assert (167, 931) in ranges

    def test_thousands_separator_and_multiple(self):
        ranges, _ = bp.parse_regions("Residues 1-552 (N-terminal) and 2,101-2,458 (central).")
        assert (1, 552) in ranges and (2101, 2458) in ranges

    def test_aa_residues_three_and_one_letter(self):
        _, res = bp.parse_regions("Residues R107, K110, Arg158, Tyr268 are critical.")
        assert {107, 110, 158, 268} <= set(res)

    def test_codon_list_with_markdown_transcript(self):
        # Markdown-escaped transcript must not break the codon-list parse.
        _, res = bp.parse_regions(
            "Missense within codons using transcript NM\\_00546.4: 175, 245, 248, 249, 273, 282."
        )
        assert set(res) == {175, 245, 248, 249, 273, 282}

    def test_occurrence_range_skipped(self):
        ranges, _ = bp.parse_regions("seen in cancerhotspots.org with 2-9 somatic occurrences")
        assert ranges == []

    def test_enumeration_noise_dropped(self):
        # "Exons 1-3 and exons 1-4" / "Cys2-Cys3" must not become residues/ranges.
        ranges, res = bp.parse_regions("Exons 1-3 and exons 1-4; critical Gly between Cys2-Cys3.")
        assert ranges == [] and res == []

    def test_bracketed_aa_ranges_preserved(self):
        # RASopathy-style "[AA 10-17]" notation: brackets carry the real
        # hotspot ranges and must NOT be stripped as citation noise.
        ranges, _ = bp.parse_regions(
            "Applicable only to domains in the supplementary table "
            "(P-loop \\[AA 10-17\\], SW1 \\[AA 25-40\\], SW2 \\[AA 57-64\\], "
            "SAK \\[AA 145-156\\]). Not applicable to specific amino acid residues (see PM5)."
        )
        assert (10, 17) in ranges and (25, 40) in ranges
        assert (57, 64) in ranges and (145, 156) in ranges

    def test_aa_prefixed_single_residues(self):
        # PTPN11-style isolated "AA 247" residues alongside bracketed ranges.
        _, res = bp.parse_regions(
            "Interacting residues \\[AA 7-9, AA 247, AA 251, AA 256\\]."
        )
        assert {247, 251, 256} <= set(res)

    def test_numeric_citation_brackets_still_stripped(self):
        # Pure-number citation brackets ("[12, 13]", "[1-3]") remain noise.
        ranges, res = bp.parse_regions("Critical domain per references \\[12, 13\\] \\[1-3\\].")
        assert ranges == [] and res == []


def _tsv(tmp_path):
    p = tmp_path / "pm1_hotspots.tsv"
    p.write_text(
        "gene_symbol\tstrength\tregions\tresidues\n"
        "MYH7\tModerate\t167-931\t\n"
        "RUNX1\tStrong\t\t107,110,134\n"
        "RUNX1\tSupporting\t89-204\t\n"
        "GAA\tModerate\t\t282,404\n"
        "ABCA4\tnot_applicable\t\t\n",
        encoding="utf-8",
    )
    return p


class TestPM1HotspotsLoader:
    def test_range_lookup(self, tmp_path):
        h = PM1Hotspots(_tsv(tmp_path))
        assert h.lookup("MYH7", 500) == CriterionStrength.MODERATE
        assert h.lookup("MYH7", 1000) is None

    def test_residue_lookup(self, tmp_path):
        h = PM1Hotspots(_tsv(tmp_path))
        assert h.lookup("RUNX1", 107) == CriterionStrength.STRONG   # residue (Strong)
        assert h.lookup("RUNX1", 150) == CriterionStrength.SUPPORTING  # 89-204 range

    def test_strongest_wins(self, tmp_path):
        h = PM1Hotspots(_tsv(tmp_path))
        # 110 is a Strong residue AND inside the Supporting 89-204 range.
        assert h.lookup("RUNX1", 110) == CriterionStrength.STRONG

    def test_not_applicable_and_has_gene(self, tmp_path):
        h = PM1Hotspots(_tsv(tmp_path))
        assert h.is_not_applicable("ABCA4") is True
        assert h.has_gene("MYH7") is True
        assert h.has_gene("ABCA4") is False     # not_applicable, no positive rows
        assert h.has_gene("UNSEEN") is False

    def test_missing_file(self, tmp_path):
        h = PM1Hotspots(tmp_path / "nope.tsv")
        assert h.lookup("MYH7", 500) is None
        assert h.is_not_applicable("ABCA4") is False


def _cfg(tmp_path, clinvar=None):
    cfg = MagicMock()
    cfg.pm1_hotspots_tsv = _tsv(tmp_path)
    cfg.clinvar_sqlite = clinvar or (tmp_path / "absent.sqlite")
    return cfg


def _ann(gene, pos, ctype=ConsequenceType.MISSENSE):
    return AnnotationData(
        gnomad=GnomADData(),
        consequences=[ConsequenceInfo(
            transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene, consequence=ctype,
            biotype="protein_coding", is_mane_select=True, protein_position=pos,
        )],
    )


def _snv():
    return VariantRecord(chrom="chr1", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38)


class TestPM1Evaluator:
    def test_curated_range_fires_at_strength(self, tmp_path):
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("MYH7", 500))
        assert r.triggered and r.strength == CriterionStrength.MODERATE
        assert "VCEP PM1 hotspot" in r.evidence

    def test_curated_miss_withheld_no_fallback(self, tmp_path):
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("MYH7", 1000))
        assert not r.triggered and "not in a VCEP PM1 hotspot" in r.evidence

    def test_not_applicable_blocks(self, tmp_path):
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("ABCA4", 100))
        assert not r.triggered and "not applicable" in r.evidence

    def test_residue_strong(self, tmp_path):
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("RUNX1", 134))
        assert r.triggered and r.strength == CriterionStrength.STRONG

    def test_non_missense_not_met(self, tmp_path):
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("MYH7", 500, ConsequenceType.SYNONYMOUS)
        )
        assert not r.triggered and "Not a missense/in-frame" in r.evidence

    def test_inframe_indel_withheld_for_missense_only_gene(self, tmp_path):
        # MYH7's VCEP scopes PM1 to missense → an in-frame indel earns nothing.
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("MYH7", 500, ConsequenceType.INFRAME_DELETION)
        )
        assert not r.triggered and "missense variants only" in r.evidence

    def test_inframe_indel_eligible_for_indel_gene(self, tmp_path):
        # GAA's VCEP extends PM1 to in-frame deletions → fires at a hotspot residue.
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("GAA", 282, ConsequenceType.INFRAME_DELETION)
        )
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_runx1_inframe_indel_withheld(self, tmp_path):
        # RUNX1 PM1 is missense-only: an in-frame indel at an RHD residue (134)
        # that WOULD be a Strong missense hotspot must not earn PM1.
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("RUNX1", 134, ConsequenceType.INFRAME_DELETION)
        )
        assert not r.triggered and "missense variants only" in r.evidence

    def test_uncurated_gene_heuristic_off_by_default(self, tmp_path):
        r = PM1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("NOVCEP", 100))
        assert not r.triggered and "heuristic disabled" in r.evidence

    def test_uncurated_gene_uses_fallback_when_enabled(self, tmp_path):
        # Opt-in: no curated rows + absent ClinVar DB -> heuristic runs and finds
        # nothing, but crucially it reached the fallback.
        cfg = _cfg(tmp_path)
        cfg.pm1_heuristic = True
        r = PM1Evaluator(cfg).evaluate(_snv(), _ann("NOVCEP", 100))
        assert not r.triggered and "hotspot cluster" in r.evidence


class TestCuratedMultiGeneSpecs:
    """PM1 definitions from multi-gene specs that the miner skips (Genome
    Medicine revision): Hearing Loss, Brain Malformations, SCID (IL2RG),
    Hemoglobinopathy."""

    def _hot(self):
        tsv = Path(__file__).resolve().parents[2] / "resources" / "shared" / "pm1_hotspots.tsv"
        return PM1Hotspots(tsv)

    def test_hearing_loss_only_kcnq4(self):
        h = self._hot()
        for g in ("CDH23", "COCH", "GJB2", "MYO15A", "MYO6", "MYO7A", "OTOF",
                  "SLC26A4", "TECTA", "USH2A"):
            assert h.is_not_applicable(g), g
        assert h.lookup("KCNQ4", 280) is not None
        assert h.lookup("KCNQ4", 300) is None

    def test_brain_malformation_domains(self):
        h = self._hot()
        from acmg_classifier.models.enums import CriterionStrength
        assert h.lookup("MTOR", 1500) == CriterionStrength.SUPPORTING
        assert h.lookup("MTOR", 2000) is None
        assert h.lookup("PIK3CA", 1047) is not None
        assert h.lookup("AKT3", 17) is not None
        assert h.lookup("PIK3R2", 400) is not None

    def test_il2rg_residues_and_polar_tm(self):
        from acmg_classifier.models.enums import CriterionStrength
        h = self._hot()
        assert h.lookup("IL2RG", 62) == CriterionStrength.STRONG
        assert h.lookup("IL2RG", 239) == CriterionStrength.STRONG
        # transmembrane 263-283: only charged/polar residues introduced
        assert h.lookup("IL2RG", 270, "D") == CriterionStrength.STRONG
        assert h.lookup("IL2RG", 270, "V") is None
        assert h.lookup("IL2RG", 270) is None

    def test_hemoglobin(self):
        h = self._hot()
        assert h.lookup("HBA2", 104) is not None and h.lookup("HBA2", 105) is None
        assert h.is_not_applicable("HBB")


def test_alt_restricted_row(tmp_path):
    from acmg_classifier.models.enums import CriterionStrength
    p = tmp_path / "pm1_hotspots.tsv"
    p.write_text("gene_symbol\tstrength\tregions\tresidues\talt_aa\n"
                 "G1\tModerate\t10-20\t\t\nG1\tStrong\t30-40\t\tDE\n", encoding="utf-8")
    h = PM1Hotspots(p)
    assert h.lookup("G1", 15) == CriterionStrength.MODERATE
    assert h.lookup("G1", 35, "E") == CriterionStrength.STRONG
    assert h.lookup("G1", 35, "A") is None


class TestPM1CoRequirements:
    """IL2RG: PM1 only together with PM2 and without BA1/BS1/BS2."""

    def _results(self, *codes):
        from acmg_classifier.models.criteria import CriteriaResult
        from acmg_classifier.models.enums import ACMGCriterion
        out = [CriteriaResult.met(ACMGCriterion.PM1, strength=CriterionStrength.STRONG, evidence="x")]
        for c in codes:
            out.append(CriteriaResult.met(ACMGCriterion(c), evidence="y"))
        return out

    def _run(self, gene, *codes):
        from acmg_classifier.criteria.registry import _apply_pm1_co_requirements
        res = self._results(*codes)
        _apply_pm1_co_requirements(res, _ann(gene, 62))
        return res[0]

    def test_il2rg_with_pm2_kept(self):
        assert not self._run("IL2RG", "PM2").suppressed

    def test_il2rg_without_pm2_suppressed(self):
        r = self._run("IL2RG")
        assert r.suppressed and "requires PM2" in r.evidence

    def test_il2rg_with_bs1_suppressed(self):
        r = self._run("IL2RG", "PM2", "BS1")
        assert r.suppressed and "not with BS1" in r.evidence

    def test_other_gene_unaffected(self):
        assert not self._run("MYH7").suppressed
