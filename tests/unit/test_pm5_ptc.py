"""PM5 for protein-truncating variants (PM5_PTC; criteria/pm5_ptc.py).

Group A (Supporting): RUNX1 (downstream of c.98, SpliceAI <= 0.20), CDH1 (NMD
predicted or PTC upstream of p.Glu836Ter), ATM (PTC upstream of p.Arg3047), PALB2 (PTC upstream of p.Tyr1183).
Group B: BRCA1/BRCA2 per-exon strength from ENIGMA Specifications Table 4,
only together with PVS1.
"""
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
from acmg_classifier.criteria.pm5_genes import PM5Spec
from acmg_classifier.criteria.pm5_ptc import brca_exon_strength, evaluate_pm5_ptc
from acmg_classifier.criteria.registry import CriteriaRegistry
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, GnomADData, SpliceScore,
)
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion, Assembly, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.variant import VariantRecord

S = CriterionStrength
_TX = {
    "RUNX1": "NM_001754.5", "CDH1": "NM_004360.5", "ATM": "NM_000051.4",
    "PALB2": "NM_024675.4", "BRCA1": "NM_007294.4", "BRCA2": "NM_000059.4",
    "TP53": "NM_000546.6",
}


def _pc(gene, hgvs_p, hgvs_c=None, cons=ConsequenceType.STOP_GAINED, tx=None, pos=None):
    if pos is None:
        import re
        m = re.search(r"(\d+)", hgvs_p)
        pos = int(m.group(1)) if m else None
    return ConsequenceInfo(
        transcript_id=tx or _TX[gene], gene_id="ENSG", gene_symbol=gene,
        consequence=cons, biotype="protein_coding", is_mane_select=True,
        protein_position=pos, hgvs_p=f"NP_0:{hgvs_p}",
        hgvs_c=(f"{tx or _TX[gene]}:{hgvs_c}" if hgvs_c else None),
    )


def _ann(pc, splice=None):
    return AnnotationData(gnomad=GnomADData(), consequences=[pc], splice=splice)


def _eval(pc, splice=None):
    return evaluate_pm5_ptc(pc, _ann(pc, splice))


class TestScope:
    def test_gene_without_rule_returns_none(self):
        assert _eval(_pc("TP53", "p.Arg196Ter", "c.586C>T")) is None

    def test_missense_returns_none(self):
        assert _eval(_pc("ATM", "p.Arg100Cys", "c.298C>T",
                         cons=ConsequenceType.MISSENSE)) is None

    def test_splice_not_automated(self):
        assert _eval(_pc("PALB2", "p.?", "c.48+1G>A",
                         cons=ConsequenceType.SPLICE_DONOR, pos=None)) is None

    def test_other_transcript_not_met(self):
        s, ev = _eval(_pc("ATM", "p.Arg100Ter", "c.298C>T", tx="NM_999999.1"))
        assert s == S.NOT_MET and "NM_000051" in ev


class TestGroupA:
    def test_atm_upstream(self):
        assert _eval(_pc("ATM", "p.Arg100Ter", "c.298C>T"))[0] == S.SUPPORTING

    def test_atm_at_or_after_r3047(self):
        assert _eval(_pc("ATM", "p.Arg3047Ter", "c.9139C>T"))[0] == S.NOT_MET

    def test_atm_frameshift_uses_new_stop(self):
        # p.Leu3040fsTer10 -> PTC codon 3049 (not upstream of R3047)
        pc = _pc("ATM", "p.Leu3040ProfsTer10", "c.9118del", cons=ConsequenceType.FRAMESHIFT)
        assert _eval(pc)[0] == S.NOT_MET
        pc = _pc("ATM", "p.Leu3000ProfsTer10", "c.8998del", cons=ConsequenceType.FRAMESHIFT)
        assert _eval(pc)[0] == S.SUPPORTING

    def test_frameshift_without_stop_not_met(self):
        pc = _pc("ATM", "p.Leu3000ProfsTer?", "c.8998del", cons=ConsequenceType.FRAMESHIFT)
        assert _eval(pc)[0] == S.NOT_MET

    def test_palb2_cutoff(self):
        assert _eval(_pc("PALB2", "p.Tyr1182Ter", "c.3546C>G"))[0] == S.SUPPORTING
        assert _eval(_pc("PALB2", "p.Tyr1183Ter", "c.3549C>G"))[0] == S.NOT_MET

    def test_runx1_downstream_of_c98(self):
        assert _eval(_pc("RUNX1", "p.Arg80Ter", "c.238C>T"))[0] == S.SUPPORTING

    def test_runx1_upstream_of_c98_not_met(self):
        assert _eval(_pc("RUNX1", "p.Gln20Ter", "c.58C>T"))[0] == S.NOT_MET
        pc = _pc("RUNX1", "p.Ser33fsTer5", "c.98del", cons=ConsequenceType.FRAMESHIFT)
        assert _eval(pc)[0] == S.NOT_MET

    def test_runx1_splice_impact_not_met(self):
        pc = _pc("RUNX1", "p.Arg80Ter", "c.238C>T")
        sp = SpliceScore(tool="spliceai", max_delta=0.45, is_available=True)
        assert _eval(pc, sp)[0] == S.NOT_MET
        sp = SpliceScore(tool="spliceai", max_delta=0.10, is_available=True)
        assert _eval(pc, sp)[0] == S.SUPPORTING

    def test_cdh1_nmd(self):
        # Early nonsense -> NMD predicted -> Supporting
        assert _eval(_pc("CDH1", "p.Arg100Ter", "c.298C>T"))[0] == S.SUPPORTING

    def test_cdh1_nmd_escape_not_met(self):
        # Last exon (16) of NM_004360 -> NMD escape, downstream of p.Glu836Ter
        assert _eval(_pc("CDH1", "p.Glu870Ter", "c.2608G>T"))[0] == S.NOT_MET
        assert _eval(_pc("CDH1", "p.Trp865Ter", "c.2594G>A"))[0] == S.NOT_MET

    def test_cdh1_nmd_escape_upstream_of_last_pathogenic(self):
        # NMD escape, but PTC upstream of p.Glu836Ter -> Supporting
        s, ev = _eval(_pc("CDH1", "p.Tyr835Ter", "c.2505T>G"))
        assert s == S.SUPPORTING and "p.Glu836Ter" in ev
        pc = _pc("CDH1", "p.Leu831AlafsTer4", "c.2490dup", cons=ConsequenceType.FRAMESHIFT)
        assert _eval(pc)[0] == S.SUPPORTING          # PTC 834

    def test_cdh1_reference_variant_and_downstream_ptc(self):
        # p.Glu836Ter itself is not upstream of itself
        assert _eval(_pc("CDH1", "p.Glu836Ter", "c.2506G>T"))[0] == S.NOT_MET
        # frameshift at codon 836 whose new stop is 846 -> not met
        pc = _pc("CDH1", "p.Glu836ValfsTer11", "c.2505_2506dup", cons=ConsequenceType.FRAMESHIFT)
        assert _eval(pc)[0] == S.NOT_MET


class TestBrcaTable4:
    @pytest.mark.parametrize("gene,codon,exon,strength", [
        ("BRCA1", 1, 2, S.STRONG),            # c.1
        ("BRCA1", 150, 7, S.SUPPORTING),      # c.448 (exon 7 c.442-547)
        ("BRCA1", 190, 8, None),              # c.568 (exon 8 c.548-593)
        ("BRCA1", 210, 9, None),              # c.628 (exon 9 c.594-670)
        ("BRCA1", 1000, 10, S.STRONG),
        ("BRCA1", 1854, 23, S.STRONG),        # last exon, < p.I1855
        ("BRCA1", 1855, 23, None),            # last exon, > p.L1854
        ("BRCA2", 120, 4, S.MODERATE),        # c.358 (exon 4 c.317-425)
        ("BRCA2", 165, 6, None),              # c.493 (exon 6 c.476-516)
        ("BRCA2", 2300, 12, None),            # c.6898 (exon 12 c.6842-6937)
        ("BRCA2", 2900, 21, S.SUPPORTING),    # c.8698 (exon 21 c.8633-8754)
        ("BRCA2", 3200, 26, S.SUPPORTING),    # c.9598 (exon 26 c.9502-9648)
        ("BRCA2", 3309, 27, S.STRONG),
        ("BRCA2", 3310, 27, None),
    ])
    def test_exon_strength(self, gene, codon, exon, strength):
        assert brca_exon_strength(gene, codon) == (exon, strength)

    def test_exon_boundary(self):
        # BRCA1 exon 10 ends at c.4096: codon 1365 starts at c.4093 (exon 10),
        # codon 1366 at c.4096 (exon 10), codon 1367 at c.4099 (exon 11).
        assert brca_exon_strength("BRCA1", 1366)[0] == 10
        assert brca_exon_strength("BRCA1", 1367)[0] == 11

    def test_frameshift_uses_ptc_exon(self):
        # Variant in exon 7 (codon 150) but new stop in exon 10 -> Strong
        pc = _pc("BRCA1", "p.Glu150GlyfsTer80", "c.448del", cons=ConsequenceType.FRAMESHIFT)
        s, ev = _eval(pc)
        assert s == S.STRONG and "exon 10" in ev

    def test_evaluator_moderate(self, tmp_path):
        cfg = MagicMock()
        cfg.disease_prevalence_tsv = tmp_path / "none.tsv"
        cfg.pm5_ptc = True
        pc = _pc("BRCA2", "p.Gln120Ter", "c.358C>T")
        r = PM5Evaluator(cfg).evaluate(
            VariantRecord(chrom="chr13", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38),
            _ann(pc),
        )
        assert r.triggered and r.strength == S.MODERATE

    def test_evaluator_disabled(self, tmp_path):
        cfg = MagicMock()
        cfg.disease_prevalence_tsv = tmp_path / "none.tsv"
        cfg.pm5_ptc = False
        pc = _pc("BRCA2", "p.Gln120Ter", "c.358C>T")
        r = PM5Evaluator(cfg).evaluate(
            VariantRecord(chrom="chr13", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38),
            _ann(pc),
        )
        assert not r.triggered


class TestBrcaRequiresPvs1:
    def _reg(self, tmp_path):
        reg = object.__new__(CriteriaRegistry)
        reg._pm5 = PM5Spec(tmp_path / "none.tsv")
        return reg

    def _pm5(self):
        return CriteriaResult.met(ACMGCriterion.PM5, S.STRONG, "PM5_PTC (BRCA1, ENIGMA Table 4): x")

    def test_suppressed_without_pvs1(self, tmp_path):
        results = [self._pm5(), CriteriaResult.not_met(ACMGCriterion.PVS1)]
        self._reg(tmp_path)._apply_pm5_ptc_requires_pvs1(
            results, _ann(_pc("BRCA1", "p.Arg100Ter", "c.298C>T")))
        assert results[0].suppressed and "requires PVS1" in results[0].evidence

    def test_kept_with_pvs1(self, tmp_path):
        results = [self._pm5(), CriteriaResult.met(ACMGCriterion.PVS1, S.STRONG)]
        self._reg(tmp_path)._apply_pm5_ptc_requires_pvs1(
            results, _ann(_pc("BRCA1", "p.Arg100Ter", "c.298C>T")))
        assert not results[0].suppressed

    def test_manual_pm5_untouched(self, tmp_path):
        results = [CriteriaResult.met(ACMGCriterion.PM5, S.STRONG, "[manual] eRepo"),
                   CriteriaResult.not_met(ACMGCriterion.PVS1)]
        self._reg(tmp_path)._apply_pm5_ptc_requires_pvs1(
            results, _ann(_pc("BRCA1", "p.Arg100Ter", "c.298C>T")))
        assert not results[0].suppressed

    def test_group_a_not_affected(self, tmp_path):
        results = [CriteriaResult.met(ACMGCriterion.PM5, S.SUPPORTING, "PM5_PTC (ATM ...)"),
                   CriteriaResult.not_met(ACMGCriterion.PVS1)]
        self._reg(tmp_path)._apply_pm5_ptc_requires_pvs1(
            results, _ann(_pc("ATM", "p.Arg100Ter", "c.298C>T")))
        assert not results[0].suppressed
