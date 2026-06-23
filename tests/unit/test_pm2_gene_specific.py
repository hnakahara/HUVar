"""PM2 gene-specific cSpec wording (F8/F9, RYR1, ATM, PTEN, RUNX1).

Each rule is transcribed verbatim from the gene's cSpec and applied only to that
gene:
  F8/F9  — "absent in males" (hemizygous count).
  RYR1   — AD "absent (1 allele allowed)".
  ATM    — "n=1 in a single subpopulation -> PM2 applies".
  PTEN   — "<1e-5 global; subpopulation with >=2 alleles must be <2e-5".
  RUNX1  — GrpMax FAF when available, else all-subpopulation (point) threshold.
"""
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_HEADER = "gene_symbol\tpm2_threshold\tpm2_strength\tpm2_basis\tpm2_subpop\n"
_ROWS = (
    "F8\t0\t\t\t\n"
    "RYR1\t0\t\t\t\n"
    "ATM\t0.00001\t\t\t\n"
    "PTEN\t0.00001\t\t\t\n"
    "RUNX1\t0.00005\t\tfaf\tpoint\n"
)


def _cfg(tmp_path):
    p = tmp_path / "dp.tsv"
    p.write_text(_HEADER + _ROWS, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = p
    cfg.gene_inheritance_tsv = tmp_path / "missing.tsv"
    return cfg


def _ann(gene, **gd):
    return AnnotationData(
        gnomad=GnomADData(filter_pass=True, **gd),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        )],
    )


def _snv():
    return VariantRecord(chrom="chrX", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


class TestMaleAbsent:
    def test_f8_absent_in_males_fires(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # Present in females (ac=5) but no hemizygotes → absent in males → PM2.
        r = ev.evaluate(_snv(), _ann("F8", ac=5, nhemi=0, popmax_af=1e-5, af=1e-5))
        assert r.triggered and "absent in males" in r.evidence

    def test_f8_present_in_males_not_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("F8", ac=5, nhemi=2, popmax_af=1e-5, af=1e-5))
        assert not r.triggered


class TestAllowOneAllele:
    def test_ryr1_one_allele_allowed(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("RYR1", ac=1, popmax_af=5e-6, af=5e-6))
        assert r.triggered and "1 allele allowed" in r.evidence

    def test_ryr1_two_alleles_not_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("RYR1", ac=2, popmax_af=5e-6, af=5e-6))
        assert not r.triggered


class TestATMSingleAllele:
    def test_atm_single_allele_overrides_inflated_point(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # popmax point (2e-5) exceeds 1e-5, but a single subpop allele → PM2.
        r = ev.evaluate(_snv(), _ann("ATM", ac=1, ac_grpmax=1, popmax_af=2e-5, af=6e-7))
        assert r.triggered and "single allele" in r.evidence

    def test_atm_multiple_alleles_not_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("ATM", ac=3, ac_grpmax=3, popmax_af=2e-5, af=2e-5))
        assert not r.triggered


class TestPTENSubpop:
    def test_pten_single_allele_judged_on_global(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # Inflated popmax point (2e-5) but global AF (5e-6) < 1e-5 and single
        # allele in the subpop → PM2.
        r = ev.evaluate(_snv(), _ann("PTEN", ac=1, ac_grpmax=1, popmax_af=2e-5, af=5e-6))
        assert r.triggered and "global AF" in r.evidence

    def test_pten_multi_allele_subpop_above_cutoff_not_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("PTEN", ac=4, ac_grpmax=4, popmax_af=3e-5, af=8e-6))
        assert not r.triggered

    def test_pten_multi_allele_subpop_below_cutoff_met(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        r = ev.evaluate(_snv(), _ann("PTEN", ac=2, ac_grpmax=2, popmax_af=1.5e-5, af=8e-6))
        assert r.triggered and "subpop AF" in r.evidence


class TestRUNX1FAFPriority:
    def test_faf_available_overrides_point_block(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # GrpMax FAF (1e-5) < 5e-5 → PM2 applies even though the POINT AF
        # (1e-4) exceeds the threshold (point block is only the FAF-absent fallback).
        r = ev.evaluate(_snv(), _ann("RUNX1", ac=3, faf95_popmax=1e-5, popmax_af=1e-4))
        assert r.triggered and "FAF95" in r.evidence

    def test_faf_unavailable_requires_all_subpop(self, tmp_path):
        ev = PM2Evaluator(_cfg(tmp_path))
        # No GrpMax FAF → fall back to the point AF (1e-4) which exceeds 5e-5.
        r = ev.evaluate(_snv(), _ann("RUNX1", ac=3, faf95_popmax=None, popmax_af=1e-4))
        assert not r.triggered
