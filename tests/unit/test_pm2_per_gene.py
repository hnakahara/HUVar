"""PM2 per-gene VCEP thresholds / strength / FAF basis (disease_prevalence.tsv)."""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
from acmg_classifier.criteria.pm2_genes import PM2Spec
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_HEADER = "gene_symbol\tpm2_threshold\tpm2_strength\tpm2_basis\n"
_ROWS = (
    "LDLR\t0.0002\tModerate\t\n"      # Moderate strength, raw-AF threshold
    "GCK\t0.000003\t\tfaf\n"           # FAF-basis, Supporting
    "KRAS\t0\t\t\n"                    # must be ABSENT
    "MYO15A\t0.00007\t\t\n"            # plain per-gene threshold
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
