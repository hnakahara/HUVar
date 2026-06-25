"""Per-gene auxiliary BayesDel/CADD/2-of-3 PP3/BP4 rules (Phase 2)."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.insilico_genes import (
    InSilicoGeneSpec, InSilicoScores, combo_active,
)
from acmg_classifier.models.annotation import (
    AnnotationData, CADDData, ConsequenceInfo, RevelData,
)
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, CriterionStrength, InSilicoTool,
)
from acmg_classifier.models.variant import VariantRecord

SPEC = InSilicoGeneSpec()


def _sc(consequence=ConsequenceType.MISSENSE, pp=None, revel=None, am=None,
        cadd=None, bayesdel=None, splice=None, hgvs_c=None, hgvs_p=None) -> InSilicoScores:
    return InSilicoScores(
        consequence=consequence, protein_position=pp, revel=revel,
        alphamissense=am, cadd=cadd, bayesdel=bayesdel, splice_max=splice,
        hgvs_c=hgvs_c, hgvs_p=hgvs_p,
    )


def _cfg(tool=InSilicoTool.REVEL, use_bayesdel=False, use_cadd=False):
    cfg = MagicMock()
    cfg.insilico_tool = tool
    cfg.use_bayesdel = use_bayesdel
    cfg.use_cadd = use_cadd
    return cfg


# --- gating -----------------------------------------------------------------

class TestComboActive:
    def test_inactive_under_esm1b_even_when_flag_on(self):
        # Licence gate: never under ESM1B (the commercial-safe path).
        assert not combo_active(SPEC, "CTLA4", _cfg(InSilicoTool.ESM1B, use_cadd=True))
        assert not combo_active(SPEC, "BRCA1", _cfg(InSilicoTool.ESM1B, use_bayesdel=True))

    def test_requires_matching_aux_flag(self):
        # CADD gene needs use_cadd; BayesDel gene needs use_bayesdel.
        assert combo_active(SPEC, "CTLA4", _cfg(use_cadd=True))
        assert not combo_active(SPEC, "CTLA4", _cfg(use_cadd=False))
        assert combo_active(SPEC, "BRCA1", _cfg(use_bayesdel=True))
        assert not combo_active(SPEC, "BRCA1", _cfg(use_bayesdel=False))

    def test_active_under_alphamissense(self):
        assert combo_active(SPEC, "BMPR2", _cfg(InSilicoTool.ALPHAMISSENSE, use_cadd=True))

    def test_uncovered_gene(self):
        assert not combo_active(SPEC, "MYH7", _cfg(use_cadd=True, use_bayesdel=True))
        assert SPEC.requires("MYH7") is None


# --- two-tool agreement (CTLA4 / PIK3CD / PIK3R1) ---------------------------

class TestTwoToolAgreement:
    def test_ctla4_pp3_requires_both_tools(self):
        # REVEL alone (≥0.75) must NOT meet PP3 — authoritative not-met.
        strength, _ = SPEC.pp3("CTLA4", _sc(revel=0.90, cadd=10.0))
        assert strength is None
        # Both agree → Supporting.
        strength, _ = SPEC.pp3("CTLA4", _sc(revel=0.90, cadd=25.0))
        assert strength == CriterionStrength.SUPPORTING

    def test_ctla4_bp4_strict_inequality(self):
        # CTLA4 BP4 uses strict < (REVEL<0.25 ∧ CADD<20).
        assert SPEC.bp4("CTLA4", _sc(revel=0.25, cadd=10.0))[0] is None
        assert SPEC.bp4("CTLA4", _sc(revel=0.10, cadd=10.0))[0] == CriterionStrength.SUPPORTING

    def test_pik3r1_thresholds(self):
        assert SPEC.pp3("PIK3R1", _sc(revel=0.70, cadd=26.0))[0] == CriterionStrength.SUPPORTING
        assert SPEC.pp3("PIK3R1", _sc(revel=0.70, cadd=25.0))[0] is None  # CADD below 26.0

    def test_bp4_blocked_by_splice_impact(self):
        # Predicted splice impact blocks the benign call.
        assert SPEC.bp4("PIK3CD", _sc(revel=0.10, cadd=10.0, splice=0.30))[0] is None
        assert SPEC.bp4("PIK3CD", _sc(revel=0.10, cadd=10.0))[0] == CriterionStrength.SUPPORTING


# --- BMPR2 2-of-3 -----------------------------------------------------------

class TestBMPR2:
    def test_pp3_two_of_three(self):
        # CADD + REVEL pass (AM missing) → 2/3 met.
        assert SPEC.pp3("BMPR2", _sc(cadd=26.0, revel=0.70))[0] == CriterionStrength.SUPPORTING
        # Only one passes → not met.
        assert SPEC.pp3("BMPR2", _sc(cadd=26.0, revel=0.50))[0] is None

    def test_bp4_synonymous_cadd_only(self):
        sc = _sc(consequence=ConsequenceType.SYNONYMOUS, cadd=20.0)
        assert SPEC.bp4("BMPR2", sc)[0] == CriterionStrength.SUPPORTING
        sc = _sc(consequence=ConsequenceType.SYNONYMOUS, cadd=23.0)
        assert SPEC.bp4("BMPR2", sc)[0] is None


# --- BRCA1/2 BayesDel, domain-gated ----------------------------------------

class TestBRCA:
    def test_pp3_in_domain(self):
        # BRCA1 RING domain (aa 2-101), BayesDel ≥ 0.28.
        assert SPEC.pp3("BRCA1", _sc(pp=50, bayesdel=0.30))[0] == CriterionStrength.SUPPORTING
        # Below cutoff → not met.
        assert SPEC.pp3("BRCA1", _sc(pp=50, bayesdel=0.20))[0] is None

    def test_pp3_outside_domain_not_met(self):
        assert SPEC.pp3("BRCA1", _sc(pp=500, bayesdel=0.90))[0] is None

    def test_brca2_bp4_in_domain(self):
        # BRCA2 DNA-binding domain (aa 2481-3186), BayesDel ≤ 0.18, no splice.
        assert SPEC.bp4("BRCA2", _sc(pp=2500, bayesdel=0.10))[0] == CriterionStrength.SUPPORTING
        assert SPEC.bp4("BRCA2", _sc(pp=2500, bayesdel=0.10, splice=0.30))[0] is None


# --- TP53 (VCEP precomputed code table) -------------------------------------

def _tp53_spec(tmp_path):
    from acmg_classifier.criteria.insilico_genes import InSilicoGeneSpec, TP53Codes
    tsv = tmp_path / "tp53.tsv"
    tsv.write_text(
        "hgvs_c\thgvs_p\tagvgd\tbayesdel\tcode\n"
        "c.7G>A\tp.Glu3Lys\tClass C65\t0.4500\tPP3\n"
        "c.523C>T\tp.Arg175Cys\tClass C65\t0.5500\tPP3_moderate\n"
        "c.4G>C\tp.Glu2Gln\tClass C0\t0.0858\tBP4\n"
        "c.9G>T\tp.Glu3Asp\tClass C0\t0.0136\tBP4_moderate\n"
        "c.7G>C\tp.Glu3Gln\tClass C15\t0.2000\tNo evidence\n"
        # An ambiguous protein change (two codes) → resolvable by hgvs_c only.
        "c.100A>T\tp.Xaa1Yyy\tClass C65\t0.4000\tPP3\n"
        "c.101A>G\tp.Xaa1Yyy\tClass C0\t0.1000\tBP4\n",
        encoding="utf-8",
    )
    return InSilicoGeneSpec(TP53Codes(tsv))


class TestTP53:
    def test_no_table_degrades_to_unavailable(self):
        # The default spec has no TP53 table → PP3/BP4 authoritatively not met.
        strength, note = SPEC.pp3("TP53", _sc(hgvs_c="NM_000546.6:c.7G>A"))
        assert strength is None and "unavailable" in note
        assert SPEC.bp4("TP53", _sc(hgvs_c="NM_000546.6:c.4G>C"))[0] is None

    def test_pp3_codes(self, tmp_path):
        spec = _tp53_spec(tmp_path)
        strength, note = spec.pp3("TP53", _sc(hgvs_c="NM_000546.6:c.7G>A"))
        assert strength == CriterionStrength.SUPPORTING
        # Evidence explains WHY: code + the Align-GVGD class + BayesDel behind it.
        assert "PP3" in note and "Align-GVGD=Class C65" in note and "BayesDel=0.4500" in note
        assert spec.pp3("TP53", _sc(hgvs_c="c.523C>T"))[0] == CriterionStrength.MODERATE
        # A BP4-coded variant is not PP3.
        assert spec.pp3("TP53", _sc(hgvs_c="c.4G>C"))[0] is None

    def test_bp4_codes(self, tmp_path):
        spec = _tp53_spec(tmp_path)
        assert spec.bp4("TP53", _sc(hgvs_c="c.4G>C"))[0] == CriterionStrength.SUPPORTING
        assert spec.bp4("TP53", _sc(hgvs_c="c.9G>T"))[0] == CriterionStrength.MODERATE
        assert spec.bp4("TP53", _sc(hgvs_c="c.7G>C"))[0] is None  # "No evidence"

    def test_protein_fallback_and_ambiguity(self, tmp_path):
        spec = _tp53_spec(tmp_path)
        # hgvs_c missing → fall back to unambiguous protein change.
        assert spec.pp3("TP53", _sc(hgvs_p="NP_000537.3:p.Glu3Lys"))[0] == CriterionStrength.SUPPORTING
        # Ambiguous protein change must NOT resolve via hgvs_p alone.
        assert spec.pp3("TP53", _sc(hgvs_p="p.Xaa1Yyy"))[0] is None
        # …but the exact hgvs_c still resolves it.
        assert spec.pp3("TP53", _sc(hgvs_c="c.100A>T"))[0] == CriterionStrength.SUPPORTING


# --- ABCA4 (REVEL for missense via default path; CADD for synonymous/indel) --

class TestABCA4:
    def test_missense_falls_through(self):
        # Missense is handled by the default per-gene REVEL path → None.
        assert SPEC.pp3("ABCA4", _sc(consequence=ConsequenceType.MISSENSE, revel=0.9)) is None
        assert SPEC.bp4("ABCA4", _sc(consequence=ConsequenceType.MISSENSE, revel=0.1)) is None

    def test_synonymous_cadd_tiers(self):
        syn = ConsequenceType.SYNONYMOUS
        assert SPEC.pp3("ABCA4", _sc(consequence=syn, cadd=29.0))[0] == CriterionStrength.MODERATE
        assert SPEC.pp3("ABCA4", _sc(consequence=syn, cadd=26.0))[0] == CriterionStrength.SUPPORTING
        assert SPEC.pp3("ABCA4", _sc(consequence=syn, cadd=20.0))[0] is None
        assert SPEC.bp4("ABCA4", _sc(consequence=syn, cadd=15.0))[0] == CriterionStrength.MODERATE
        assert SPEC.bp4("ABCA4", _sc(consequence=syn, cadd=19.0))[0] == CriterionStrength.SUPPORTING
        assert SPEC.bp4("ABCA4", _sc(consequence=syn, cadd=25.0))[0] is None


# --- end-to-end wiring through the real PP3 / BP4 evaluators -----------------

def _eval_cfg(tool=InSilicoTool.REVEL, use_cadd=True, use_bayesdel=False):
    cfg = MagicMock()
    cfg.insilico_tool = tool
    cfg.use_cadd = use_cadd
    cfg.use_bayesdel = use_bayesdel
    # Non-existent path → RevelSpec degrades to "no per-gene rule" cleanly.
    cfg.disease_prevalence_tsv = Path("/nonexistent_spec.tsv")
    return cfg


def _ctla4_ann(revel, cadd):
    return AnnotationData(
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol="CTLA4",
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        )],
        revel=RevelData(score=revel),
        cadd=CADDData(phred=cadd),
    )


class TestEvaluatorWiring:
    def test_pp3_authoritative_blocks_revel_only_firing(self):
        from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
        ev = PP3Evaluator(_eval_cfg())
        snv = VariantRecord(chrom="2", pos=1, ref="A", alt="G", assembly=Assembly.GRCH38)
        # REVEL alone is high (would meet the default REVEL PP3) but CADD is low →
        # the CTLA4 2-tool rule authoritatively withholds PP3.
        assert not ev.evaluate(snv, _ctla4_ann(0.90, 10.0)).triggered
        # Both tools agree → PP3 Supporting.
        r = ev.evaluate(snv, _ctla4_ann(0.90, 25.0))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_pp3_inactive_without_cadd_flag_uses_default(self):
        from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
        ev = PP3Evaluator(_eval_cfg(use_cadd=False))
        snv = VariantRecord(chrom="2", pos=1, ref="A", alt="G", assembly=Assembly.GRCH38)
        # Without --with-cadd the combo is inactive → default genome-wide REVEL
        # path fires on the high REVEL alone.
        assert ev.evaluate(snv, _ctla4_ann(0.90, 10.0)).triggered

    def test_bp4_authoritative(self):
        from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
        ev = BP4Evaluator(_eval_cfg())
        snv = VariantRecord(chrom="2", pos=1, ref="A", alt="G", assembly=Assembly.GRCH38)
        # Benign-range REVEL but CADD not below 20 → CTLA4 BP4 withheld.
        assert not ev.evaluate(snv, _ctla4_ann(0.10, 25.0)).triggered
        r = ev.evaluate(snv, _ctla4_ann(0.10, 10.0))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_tp53_codes_through_evaluator_and_gate(self, tmp_path):
        from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
        tsv = tmp_path / "tp53.tsv"
        tsv.write_text(
            "hgvs_c\thgvs_p\tagvgd\tbayesdel\tcode\n"
            "c.523C>T\tp.Arg175Cys\tClass C65\t0.5500\tPP3_moderate\n",
            encoding="utf-8")
        ann = AnnotationData(consequences=[ConsequenceInfo(
            transcript_id="NM_000546.6", gene_id="ENSG", gene_symbol="TP53",
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
            hgvs_c="NM_000546.6:c.523C>T", hgvs_p="NP_000537.3:p.Arg175Cys",
        )])
        snv = VariantRecord(chrom="17", pos=1, ref="G", alt="A", assembly=Assembly.GRCH38)

        cfg = _eval_cfg(use_cadd=False, use_bayesdel=True)
        cfg.tp53_codes_tsv = tsv
        r = PP3Evaluator(cfg).evaluate(snv, ann)
        assert r.triggered and r.strength == CriterionStrength.MODERATE

        # ESM1B → licence gate off; TP53 codes are NOT consulted (and no ESM1b
        # score present → not met), never the precomputed PP3.
        cfg_esm = _eval_cfg(tool=InSilicoTool.ESM1B, use_cadd=False, use_bayesdel=True)
        cfg_esm.tp53_codes_tsv = tsv
        assert not PP3Evaluator(cfg_esm).evaluate(snv, ann).triggered
