"""PM2 per-gene VCEP thresholds / strength / FAF basis (disease_prevalence.tsv)."""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
from acmg_classifier.criteria.pm2_genes import PM2Spec
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_HEADER = "gene_symbol\tpm2_threshold\tpm2_strength\tpm2_basis\tpm2_subpop\tpm2_zygosity\n"
_ROWS = (
    "LDLR\t0.0002\tModerate\t\t\t\n"    # Moderate strength, raw-AF threshold
    "GCK\t0.000003\t\tfaf\t\t\n"         # FAF-basis, Supporting
    "KRAS\t0\t\t\t\t\n"                  # must be ABSENT
    "MYO15A\t0.00007\t\t\t\t\n"          # plain per-gene threshold
    "RUNX1\t0.00005\t\tfaf\tpoint\t\n"   # FAF-basis + all-subpopulations point rule
    "MYH7\t0.00004\t\tfaf\tci95\t\n"     # FAF-basis + HCM upper-95%-CI rule
    "SLC6A8\t0.00002\t\tfaf\t\thomhemi:0\n"  # 0 homo- or hemizygotes
    "ADA\t0.0001742\t\tfaf\t\thom:0\n"       # no homozygotes
)


def _spec_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_HEADER + _ROWS, encoding="utf-8")
    return p


def _cfg(tmp_path: Path) -> MagicMock:
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _spec_tsv(tmp_path)
    cfg.gene_inheritance_tsv = tmp_path / "missing_inheritance.tsv"  # → dominant default
    return cfg


def _consequence(gene: str) -> ConsequenceInfo:
    return ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
    )


def _ann(gene: str, **gnomad) -> AnnotationData:
    return AnnotationData(
        gnomad=GnomADData(filter_pass=True, **gnomad),
        consequences=[_consequence(gene)],
    )


def _snv() -> VariantRecord:
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


class TestUpperAf95:
    """The Poisson upper-95%-CI helper must reproduce the Cardiomyopathy/HCM
    VCEP's published AC/AN equivalence table for the 0.00004 cutoff."""

    def test_vcep_equivalence_table(self):
        from acmg_classifier.criteria.pathogenic.pm2 import _upper_af_95
        # AC≤1 in AN≥120,000 / ≤2 in ≥160,000 / ≤3 in ≥195,000 / ≤4 in ≥230,000
        # → upper bound ≈ 0.00004 (just below).
        for ac, an in [(1, 120000), (2, 160000), (3, 195000), (4, 230000)]:
            u = _upper_af_95(ac, an)
            assert 3.8e-5 <= u <= 4.0e-5, (ac, an, u)

    def test_higher_count_exceeds(self):
        from acmg_classifier.criteria.pathogenic.pm2 import _upper_af_95
        # One more allele at the same AN pushes the upper bound over 0.00004.
        assert _upper_af_95(5, 230000) > 4e-5

    def test_none_when_missing(self):
        from acmg_classifier.criteria.pathogenic.pm2 import _upper_af_95
        assert _upper_af_95(None, 100000) is None
        assert _upper_af_95(2, None) is None
        assert _upper_af_95(2, 0) is None


class TestPM2Spec:
    def test_moderate_gene(self, tmp_path):
        spec = PM2Spec(_spec_tsv(tmp_path))
        assert spec.get("LDLR").strength == CriterionStrength.MODERATE
        assert spec.get("LDLR").threshold == 0.0002
        assert spec.get("LDLR").use_faf is False

    def test_faf_gene(self, tmp_path):
        spec = PM2Spec(_spec_tsv(tmp_path))
        assert spec.get("GCK").use_faf is True
        assert spec.get("GCK").strength == CriterionStrength.SUPPORTING

    def test_unknown_gene_none(self, tmp_path):
        assert PM2Spec(_spec_tsv(tmp_path)).get("TP53") is None

    def test_subpop_mode_loaded(self, tmp_path):
        spec = PM2Spec(_spec_tsv(tmp_path))
        assert spec.get("RUNX1").subpop_mode == "point"
        assert spec.get("MYH7").subpop_mode == "ci95"
        assert spec.get("GCK").subpop_mode == ""

    def test_zygosity_loaded(self, tmp_path):
        spec = PM2Spec(_spec_tsv(tmp_path))
        assert (spec.get("SLC6A8").zyg_scope, spec.get("SLC6A8").zyg_max) == ("homhemi", 0)
        assert (spec.get("ADA").zyg_scope, spec.get("ADA").zyg_max) == ("hom", 0)
        assert spec.get("GCK").zyg_scope == ""


class TestPM2PerGeneEvaluator:
    def test_moderate_strength_emitted(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # popmax 1e-4 < 2e-4 threshold → met at Moderate.
        r = ev.evaluate(_snv(), _ann("LDLR", ac=5, af=1e-4, popmax_af=1e-4, faf95_popmax=1e-4))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_faf_basis_uses_faf_not_raw(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # Raw popmax (1e-4) is ABOVE the 3e-6 cutoff, but FAF95 (1e-6) is below.
        # FAF-basis gene must fire PM2 on the FAF value.
        r = ev.evaluate(_snv(), _ann("GCK", ac=20, af=1e-4, popmax_af=1e-4, faf95_popmax=1e-6))
        assert r.triggered and "FAF95" in r.evidence

    def test_absent_required_gene_present_not_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # threshold 0 → a present (AC>0, AF>0) variant must NOT get PM2.
        r = ev.evaluate(_snv(), _ann("KRAS", ac=10, af=5e-5, popmax_af=5e-5, faf95_popmax=5e-5))
        assert not r.triggered

    def test_absent_required_gene_absent_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("KRAS", ac=0, af=0.0, popmax_af=0.0, faf95_popmax=0.0))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_per_gene_threshold_above_not_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # MYO15A cutoff 7e-5; popmax 1e-4 is above → not met (would also exceed).
        r = ev.evaluate(_snv(), _ann("MYO15A", ac=50, af=1e-4, popmax_af=1e-4, faf95_popmax=1e-4))
        assert not r.triggered
        assert "MYO15A VCEP" in r.evidence

    def test_runx1_faf_available_is_authoritative(self, tmp_path):
        # cSpec: "evaluate PM2 using the GrpMax FAF when available." FAF95=2.1e-5
        # < 5e-5 → PM2 fires even though the GrpMax POINT AF (1.2e-4) exceeds the
        # threshold — the point/all-subpopulation rule is only the FAF-absent
        # fallback.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("RUNX1", ac=2, af=8.56e-6,
                                     popmax_af=1.206e-4, faf95_popmax=2.09e-5))
        assert r.triggered and "FAF95" in r.evidence

    def test_runx1_faf_absent_requires_all_subpop(self, tmp_path):
        # GrpMax FAF unavailable → require all subpopulations meet the threshold;
        # the POINT AF (1.2e-4) exceeds 5e-5 → PM2 must NOT fire.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("RUNX1", ac=2, af=8.56e-6,
                                     popmax_af=1.206e-4, faf95_popmax=None))
        assert not r.triggered

    def test_runx1_subpop_allows_truly_rare(self, tmp_path):
        # Both FAF95 and GrpMax point AF below threshold → PM2 fires normally.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("RUNX1", ac=1, af=1e-6,
                                     popmax_af=1e-5, faf95_popmax=8e-6))
        assert r.triggered

    def test_runx1_subpop_absent_still_met(self, tmp_path):
        # An absent variant (AC=0) still gets PM2 — the subpop check doesn't block it.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("RUNX1", ac=0, af=0.0,
                                     popmax_af=0.0, faf95_popmax=0.0))
        assert r.triggered

    def test_hcm_ci95_upper_exceeds_blocks(self, tmp_path):
        # MYH7 (ci95): FAF95 below 4e-5 would fire, but the GrpMax 95%-CI upper
        # bound exceeds it. AC=4 in AN=100000 → upper 9.154/100000 = 9.2e-5 ≥
        # 4e-5 → PM2 must NOT fire.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("MYH7", ac=4, af=2e-5, popmax_af=4e-5,
                                     faf95_popmax=1e-5, ac_grpmax=4, an_grpmax=100000))
        assert not r.triggered
        assert "95%CI upper" in r.evidence

    def test_hcm_ci95_upper_below_threshold_met(self, tmp_path):
        # AC=4 in AN=230000 → upper 9.154/230000 = 3.98e-5 < 4e-5 → PM2 fires
        # (matches the VCEP AC/AN equivalence table).
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("MYH7", ac=4, af=1.7e-5, popmax_af=1.7e-5,
                                     faf95_popmax=1e-5, ac_grpmax=4, an_grpmax=230000))
        assert r.triggered

    def test_hcm_ci95_falls_back_to_point_when_grpmax_absent(self, tmp_path):
        # Old DB without GrpMax AC/AN → fall back to the POINT-AF proxy.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("MYH7", ac=10, af=4e-5, popmax_af=8e-5,
                                     faf95_popmax=1e-5))
        assert not r.triggered
        assert "CI unavailable" in r.evidence

    def test_slc6a8_homhemi_blocks(self, tmp_path):
        # The reported FP: rare (FAF95 < 2e-5) but a homo/hemizygote is present
        # → SLC6A8 "0 homo- or hemizygotes" withholds PM2.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("SLC6A8", ac=3, af=5e-6, popmax_af=1e-5,
                                     faf95_popmax=8e-6, nhemi=1))
        assert not r.triggered
        assert "hemizygote" in r.evidence.lower()

    def test_slc6a8_no_homhemi_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("SLC6A8", ac=2, af=5e-6, popmax_af=1e-5,
                                     faf95_popmax=8e-6, nhemi=0, nhomalt=0))
        assert r.triggered

    def test_ada_hom_only_blocks(self, tmp_path):
        # ADA "no homozygotes": a homozygote present → withhold PM2 even though rare.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("ADA", ac=10, af=5e-5, popmax_af=1e-4,
                                     faf95_popmax=5e-5, nhomalt=1))
        assert not r.triggered
        assert "homozygote" in r.evidence.lower()

    def test_ada_hemi_does_not_block_hom_only(self, tmp_path):
        # ADA is hom-only; a hemizygote (irrelevant for autosomal) must not block.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("ADA", ac=5, af=5e-5, popmax_af=1e-4,
                                     faf95_popmax=5e-5, nhomalt=0, nhemi=3))
        assert r.triggered

    def test_zygosity_blocks_even_when_faf_zero(self, tmp_path):
        # FAF95≈0 (AC=0 branch) but a homozygote recorded → still withheld.
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("SLC6A8", ac=2, af=1e-6, popmax_af=0.0,
                                     faf95_popmax=0.0, nhomalt=1))
        assert not r.triggered
