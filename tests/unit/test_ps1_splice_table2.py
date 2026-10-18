"""PS1 for splicing per ClinGen SVI Walker et al. 2023 Table 2, the ±1,2
dinucleotide comparator, exonic (synonymous) variants, and the VCEP comparator
rules added in the Genome Medicine revision (LDLR/DICER1 Pathogenic-only)."""
import sqlite3
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator, SPLICE_TAG, canonical_offset
from acmg_classifier.criteria.registry import apply_ps1_splice_table2
from acmg_classifier.criteria.vcep_comparators import filter_comparators, load_rules
from acmg_classifier.local_db.clinvar_sqlite import query_splice_dinucleotide_partner
from acmg_classifier.models.annotation import AnnotationData, ClinVarRecord, ConsequenceInfo, GnomADData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, Assembly, ConsequenceType, CriterionStrength as S
from acmg_classifier.models.variant import VariantRecord

_SCHEMA = """CREATE TABLE variants (variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT, codon_position INTEGER,
clinical_significance TEXT, review_status TEXT, star_rating INTEGER);"""


def _db(tmp_path, rows):
    p = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, hgvs_c, hgvs_p,"
        " clinical_significance, review_status, star_rating) VALUES (?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit(); con.close()
    return p


def _row(vid, pos, ref, alt, hgvs_c, sig, hgvs_p=None, gene="HNF1A", stars=2):
    return (vid, "12", pos, ref, alt, gene, hgvs_c, hgvs_p, sig, "x", stars)


@pytest.fixture(autouse=True)
def _no_vcep_rules(monkeypatch, request):
    if request.cls is not None and request.cls.__name__ == "TestVcepComparatorRules":
        return
    import acmg_classifier.criteria.vcep_comparators as vc
    monkeypatch.setattr(vc, "requires_vcep_comparator", lambda *a, **k: None)


def _cfg(db, tmp_path):
    cfg = MagicMock()
    cfg.clinvar_sqlite = db
    cfg.ps1_min_stars = 2
    tsv = tmp_path / "dp.tsv"
    tsv.write_text("gene_symbol\tps1\tps1_splice\nHNF1A\tapplicable\tcanonical\n", encoding="utf-8")
    cfg.disease_prevalence_tsv = tsv
    return cfg


def _ann(cq, hgvs_c, hgvs_p=None):
    return AnnotationData(gnomad=GnomADData(), consequences=[ConsequenceInfo(
        transcript_id="NM_000545.8", gene_id="G", gene_symbol="HNF1A", consequence=cq,
        biotype="protein_coding", is_mane_select=True, hgvs_c=f"NM_000545.8:{hgvs_c}", hgvs_p=hgvs_p)])


def _v(pos, ref="G", alt="T"):
    return VariantRecord(chrom="chr12", pos=pos, ref=ref, alt=alt, assembly=Assembly.GRCH38)


class TestCanonicalOffset:
    @pytest.mark.parametrize("c,exp", [("c.526+1G>T", ("+", 1)), ("c.527-2A>G", ("-", 2)),
                                       ("c.526+12G>T", None), ("c.526+5G>A", None), ("c.100G>A", None)])
    def test_offsets(self, c, exp):
        pc = ConsequenceInfo(transcript_id="t", gene_id="g", gene_symbol="X",
                             consequence=ConsequenceType.INTRON, biotype="pc", hgvs_c=c)
        assert canonical_offset(pc) == exp


class TestDinucleotidePartner:
    def test_partner_found_either_strand(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", 101, "T", "C", "NM_1:c.526+2T>C", "Pathogenic"),   # + strand partner
            _row("2", 99, "A", "G", "NM_1:c.526+2T>C", "Pathogenic"),    # - strand partner
            _row("3", 101, "A", "G", "NM_1:c.526+3A>G", "Pathogenic"),   # not the dinucleotide
            _row("4", 99, "A", "G", "NM_1:c.526+2T>C", "Benign"),
        ])
        hits = query_splice_dinucleotide_partner(db, "HNF1A", "chr12", 100, "+", 1, min_stars=2)
        assert sorted(h.variation_id for h in hits) == ["1", "2"]

    def test_acceptor(self, tmp_path):
        db = _db(tmp_path, [_row("1", 99, "A", "C", "c.527-1G>T", "Likely pathogenic")])
        hits = query_splice_dinucleotide_partner(db, "HNF1A", "12", 100, "-", 2, min_stars=2)
        assert [h.variation_id for h in hits] == ["1"]


class TestEvaluatorTags:
    def test_canonical_partner_comparator(self, tmp_path):
        db = _db(tmp_path, [_row("9", 101, "T", "C", "c.526+2T>C", "Pathogenic")])
        r = PS1Evaluator(_cfg(db, tmp_path)).evaluate(_v(100), _ann(ConsequenceType.SPLICE_DONOR, "c.526+1G>T"))
        assert r.triggered and f"{SPLICE_TAG}canonical:P]" in r.evidence
        assert "same ±1,2 dinucleotide" in r.evidence

    def test_noncanonical_lp(self, tmp_path):
        db = _db(tmp_path, [_row("9", 100, "G", "C", "c.526+5G>C", "Likely pathogenic")])
        r = PS1Evaluator(_cfg(db, tmp_path)).evaluate(_v(100, "G", "A"), _ann(ConsequenceType.SPLICE_REGION, "c.526+5G>A"))
        assert r.strength == S.MODERATE and f"{SPLICE_TAG}noncanonical:LP]" in r.evidence

    def test_synonymous_uses_splice_comparator_only(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", 100, "G", "A", "c.525G>A", "Pathogenic", hgvs_p="p.Lys175="),
            _row("2", 100, "G", "C", "c.525G>C", "Pathogenic", hgvs_p="p.Lys175Asn"),
        ])
        r = PS1Evaluator(_cfg(db, tmp_path)).evaluate(
            _v(100, "G", "T"), _ann(ConsequenceType.SYNONYMOUS, "c.525G>T", "p.Lys175="))
        assert r.triggered and "1" in r.evidence and ", 2" not in r.evidence

    def test_synonymous_missense_comparator_only_not_met(self, tmp_path):
        db = _db(tmp_path, [_row("2", 100, "G", "C", "c.525G>C", "Pathogenic", hgvs_p="p.Lys175Asn")])
        r = PS1Evaluator(_cfg(db, tmp_path)).evaluate(
            _v(100, "G", "T"), _ann(ConsequenceType.SYNONYMOUS, "c.525G>T", "p.Lys175="))
        assert not r.triggered


def _ps1(site, comp, strength):
    return CriteriaResult.met(ACMGCriterion.PS1, strength, f"splice {SPLICE_TAG}{site}:{comp}]")


class TestWalkerTable2:
    def _run(self, ps1, pvs1=None, pp3=False):
        res = [ps1]
        if pvs1 is not None:
            res.append(CriteriaResult.met(ACMGCriterion.PVS1, pvs1, "pvs1"))
        if pp3:
            res.append(CriteriaResult.met(ACMGCriterion.PP3, S.SUPPORTING, "pp3"))
        apply_ps1_splice_table2(res)
        return ps1

    def test_canonical_pvs1_vs_p_supporting(self):
        r = self._run(_ps1("canonical", "P", S.STRONG), S.VERY_STRONG)
        assert r.strength == S.SUPPORTING and not r.suppressed

    def test_canonical_pvs1_vs_lp_na(self):
        assert self._run(_ps1("canonical", "LP", S.MODERATE), S.VERY_STRONG).suppressed

    def test_canonical_reduced_pvs1_p_strong(self):
        r = self._run(_ps1("canonical", "P", S.STRONG), S.MODERATE)
        assert r.strength == S.STRONG and not r.suppressed

    def test_canonical_reduced_pvs1_lp_na(self):
        assert self._run(_ps1("canonical", "LP", S.MODERATE), S.STRONG).suppressed

    def test_canonical_without_pvs1_na(self):
        assert self._run(_ps1("canonical", "P", S.STRONG)).suppressed

    def test_noncanonical_requires_pp3(self):
        assert self._run(_ps1("noncanonical", "P", S.STRONG)).suppressed
        r = self._run(_ps1("noncanonical", "LP", S.MODERATE), pp3=True)
        assert r.strength == S.MODERATE and not r.suppressed

    def test_untagged_ps1_untouched(self):
        ps1 = CriteriaResult.met(ACMGCriterion.PS1, S.STRONG, "ClinVar same AA")
        assert not self._run(ps1).suppressed


class TestVcepComparatorRules:
    def test_new_rules_loaded(self):
        rules = load_rules()
        for key in [("LDLR", "PM5"), ("LDLR", "PS1"), ("DICER1", "PS1"), ("RPE65", "PM5"),
                    ("GUCY2D", "PS1"), ("AIPL1", "PM5"), ("RPGR", "PS1"), ("RS1", "PM5"),
                    ("KCNQ1", "PS1"), ("GP9", "PM5")]:
            assert key in rules, key
        assert rules[("LDLR", "PM5")].pathogenic_only and rules[("DICER1", "PS1")].pathogenic_only
        assert "PS1" in rules[("KCNQ1", "PS1")].exclude_codes
        assert ("GATM", "PM5") not in rules  # curator may classify the comparator

    def test_ldlr_pathogenic_only(self):
        lp = ClinVarRecord(variation_id="1", clinical_significance="Likely pathogenic", review_status="x", star_rating=3,
                           chrom="19", pos=1, ref="A", alt="G")
        p = ClinVarRecord(variation_id="2", clinical_significance="Pathogenic", review_status="x", star_rating=3,
                          chrom="19", pos=2, ref="A", alt="G")
        low = ClinVarRecord(variation_id="3", clinical_significance="Pathogenic", review_status="x", star_rating=2,
                            chrom="19", pos=3, ref="A", alt="G")
        kept, note = filter_comparators("LDLR", "PM5", [lp, p, low], "GRCh38")
        assert [h.variation_id for h in kept] == ["2"] and note


def test_missense_only_ps1_genes():
    from pathlib import Path
    from acmg_classifier.criteria.ps1_genes import PS1Spec
    spec = PS1Spec(Path(__file__).resolve().parents[2] / "resources" / "shared" / "disease_prevalence.tsv")
    for g in ("LDLR", "FBN1", "PIK3CD", "KCNQ1"):
        assert spec.splice_mode(g) == "", g
    assert spec.splice_mode("HNF1A") == "canonical"


class TestVcepOwnSpliceRule:
    """Genes whose VCEP has its own splice PS1 rule (e.g. ABCD1: P -> Strong,
    LP canonical splice comparator at the same position -> Moderate) are not
    re-weighted by Walker Table 2 and do not use the dinucleotide comparator."""

    def _cfg(self, db, tmp_path):
        cfg = _cfg(db, tmp_path)
        cfg.disease_prevalence_tsv.write_text(
            "gene_symbol\tps1\tps1_splice\nABCD1\tapplicable\tcanonical\n", encoding="utf-8")
        return cfg

    def _ann(self, c):
        return AnnotationData(gnomad=GnomADData(), consequences=[ConsequenceInfo(
            transcript_id="NM_000033.4", gene_id="G", gene_symbol="ABCD1",
            consequence=ConsequenceType.SPLICE_DONOR, biotype="protein_coding",
            is_mane_select=True, hgvs_c=f"NM_000033.4:{c}")])

    def test_abcd1_lp_same_position_moderate_kept_with_pvs1(self, tmp_path):
        db = _db(tmp_path, [("1", "12", 100, "T", "A", "ABCD1", "c.900+2T>A", None,
                             "Likely pathogenic", "x", 3)])
        ps1 = PS1Evaluator(self._cfg(db, tmp_path)).evaluate(
            _v(100, "T", "C"), self._ann("c.900+2T>C"))
        assert ps1.triggered and ps1.strength == S.MODERATE
        assert f"{SPLICE_TAG}vcep:LP]" in ps1.evidence
        res = [ps1, CriteriaResult.met(ACMGCriterion.PVS1, S.STRONG, "pvs1")]
        apply_ps1_splice_table2(res)
        assert not ps1.suppressed and ps1.strength == S.MODERATE

    def test_abcd1_no_dinucleotide_partner(self, tmp_path):
        db = _db(tmp_path, [("1", "12", 101, "T", "A", "ABCD1", "c.900+2T>A", None,
                             "Pathogenic", "x", 3)])
        ps1 = PS1Evaluator(self._cfg(db, tmp_path)).evaluate(
            _v(100, "G", "C"), self._ann("c.900+1G>C"))
        assert not ps1.triggered


class TestBmpr2SpliceMotif:
    """BMPR2 VCEP: same donor (last 3 exonic bases, +3..+6) / acceptor (first exonic
    base, -20..-3) motif as a P (PS1_Moderate) or LP (PS1_Supporting) variant."""

    def _cfg(self, db, tmp_path):
        cfg = _cfg(db, tmp_path)
        cfg.disease_prevalence_tsv.write_text(
            "gene_symbol\tps1\tps1_splice\nBMPR2\tapplicable\tnoncanonical\n", encoding="utf-8")
        return cfg

    def _ann(self, c, cq=ConsequenceType.SPLICE_REGION):
        return AnnotationData(gnomad=GnomADData(), consequences=[ConsequenceInfo(
            transcript_id="NM_001204.7", gene_id="G", gene_symbol="BMPR2", consequence=cq,
            biotype="protein_coding", is_mane_select=True, hgvs_c=f"NM_001204.7:{c}")])

    def _row(self, vid, pos, c, sig, p=None):
        return (vid, "12", pos, "A", "G", "BMPR2", f"NM_001204.7:{c}", p, sig, "x", 2)

    def test_same_donor_motif_p_moderate(self, tmp_path):
        db = _db(tmp_path, [self._row("1", 102, "c.418+5G>A", "Pathogenic")])
        r = PS1Evaluator(self._cfg(db, tmp_path)).evaluate(_v(100, "A", "T"), self._ann("c.418+3A>T"))
        assert r.triggered and r.strength == S.MODERATE and f"{SPLICE_TAG}vcep:P]" in r.evidence
        res = [r]  # no PP3: the BMPR2 rule does not require it
        apply_ps1_splice_table2(res)
        assert not r.suppressed

    def test_same_position_lp_supporting(self, tmp_path):
        db = _db(tmp_path, [self._row("1", 100, "c.1276+3A>G", "Likely pathogenic")])
        r = PS1Evaluator(self._cfg(db, tmp_path)).evaluate(_v(100, "A", "T"), self._ann("c.1276+3A>T"))
        assert r.triggered and r.strength == S.SUPPORTING

    def test_acceptor_motif(self, tmp_path):
        db = _db(tmp_path, [self._row("1", 98, "c.968-5A>G", "Pathogenic")])
        r = PS1Evaluator(self._cfg(db, tmp_path)).evaluate(_v(100, "C", "G"), self._ann("c.968-3C>G"))
        assert r.triggered and r.strength == S.MODERATE

    def test_canonical_comparator_not_in_motif(self, tmp_path):
        db = _db(tmp_path, [self._row("1", 103, "c.418+1G>A", "Pathogenic")])
        r = PS1Evaluator(self._cfg(db, tmp_path)).evaluate(_v(100, "A", "T"), self._ann("c.418+3A>T"))
        assert not r.triggered

    def test_other_motif_not_used(self, tmp_path):
        db = _db(tmp_path, [self._row("1", 110, "c.419-5A>G", "Pathogenic")])
        r = PS1Evaluator(self._cfg(db, tmp_path)).evaluate(_v(100, "A", "T"), self._ann("c.418+3A>T"))
        assert not r.triggered

    def test_deep_intronic_not_in_motif(self, tmp_path):
        db = _db(tmp_path, [self._row("1", 101, "c.419-25A>G", "Pathogenic")])
        r = PS1Evaluator(self._cfg(db, tmp_path)).evaluate(_v(100, "A", "T"), self._ann("c.419-24A>T", ConsequenceType.INTRON))
        assert not r.triggered


def test_bmpr2_motif_exonic_positions():
    from acmg_classifier.criteria.ps1_splice_motif import motif_of
    from acmg_classifier.pvs1.last_nt import _exon_bounds
    b = _exon_bounds("NM_001204.7")
    if not b:
        pytest.skip("MANE exon structure not available")
    end = b[0][1]
    assert motif_of(f"c.{end}G>A", "NM_001204.7") == ("D", end)
    assert motif_of(f"c.{end - 2}G>A", "NM_001204.7") == ("D", end)
    assert motif_of(f"c.{end - 3}G>A", "NM_001204.7") is None
    assert motif_of(f"c.{end + 1}G>A", "NM_001204.7") == ("A", end + 1)
    assert motif_of(f"c.{end}+2T>C", "NM_001204.7") is None
