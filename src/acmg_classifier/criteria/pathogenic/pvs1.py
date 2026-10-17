"""PVS1 — null variant in gene where LoF is disease mechanism (ClinGen 2019 decision tree)."""
from __future__ import annotations
from pathlib import Path
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Gate for the PVS1 decision tree: only these VEP consequences are candidates
# for "predicted null variant" per ACMG 2015. INFRAME_INSERTION/DELETION are
# handled by PM4, not PVS1. SPLICE_REGION (±1-2 outside donor/acceptor core)
# is intentionally excluded — its impact is uncertain and is covered by PP3
# splice predictors instead.
_LOF_CONSEQUENCES = {
    ConsequenceType.FRAMESHIFT,
    ConsequenceType.STOP_GAINED,
    ConsequenceType.SPLICE_ACCEPTOR,
    ConsequenceType.SPLICE_DONOR,
    ConsequenceType.START_LOST,
    ConsequenceType.TRANSCRIPT_ABLATION,
}


class PVS1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        from acmg_classifier.criteria.pvs1_genes import PVS1Applicability
        self._spec = PVS1Applicability(cfg.disease_prevalence_tsv)
        # Optional exon-aware splice-strength overrides; empty when the TSV is
        # absent, leaving the flat per-gene VCEP splice defaults untouched.
        from acmg_classifier.pvs1.vcep_pvs1_exons import SpliceExonOverrides
        self._splice_overrides = SpliceExonOverrides(
            getattr(cfg, "vcep_pvs1_splice_exons_tsv", None)
        )

    def _fetch_seq(self, chrom: str, start: int, end: int) -> str | None:
        """Reference sequence (1-based inclusive, + strand) from the genome
        FASTA, or None when it is unavailable."""
        fa = getattr(self, "_fasta", False)
        if fa is False:
            fa = None
            try:
                path = Path(self._cfg.genome_fasta)
                if path.exists():
                    import pysam  # type: ignore
                    fa = pysam.FastaFile(str(path))
            except Exception:  # noqa: BLE001 — optional resource
                fa = None
            self._fasta = fa
        if fa is None:
            return None
        try:
            name = chrom if chrom in fa.references else (
                chrom[3:] if chrom.startswith("chr") else f"chr{chrom}")
            return fa.fetch(name, start - 1, end).upper()
        except Exception:  # noqa: BLE001
            return None

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # Cheap pre-filter before invoking the full decision tree: if the
        # primary transcript's consequence is not a candidate LoF type, PVS1
        # cannot fire and we skip the expensive NMD / last-exon / cryptic-
        # splice analysis in pvs1/decision_tree.py.
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.PVS1, "No primary consequence")

        # APC has a gene-specific PVS1 tree (InSiGHT / Tayoun): a codon-range gate
        # for truncating variants and an explicit allele-specific strength table
        # for canonical-splice / "G to non-G last nucleotide" changes. The latter
        # can fire on exonic changes VEP calls missense/synonymous, so the APC
        # handler runs BEFORE the generic LoF-consequence gate. It returns None
        # for variants its special rules don't cover (then fall back to generic).
        if pc.gene_symbol == "APC":
            from acmg_classifier.pvs1.apc import evaluate_apc_pvs1
            apc = evaluate_apc_pvs1(pc)
            if apc is not None:
                strength, evidence = apc
                if strength == CriterionStrength.NOT_MET:
                    return CriteriaResult.not_met(ACMGCriterion.PVS1, evidence)
                return CriteriaResult.met(ACMGCriterion.PVS1, strength, evidence)

        # VCEP "G to non-G at the last exon nucleotide" rules (ATM / MLH1 /
        # CDH1): exonic changes VEP calls missense / synonymous, so like APC this
        # runs before the generic LoF-consequence gate.
        from acmg_classifier.pvs1.last_nt import evaluate_last_nt
        last_nt = evaluate_last_nt(pc, variant.chrom, variant.pos, variant.ref,
                                   self._fetch_seq)
        if last_nt is not None:
            strength, evidence = last_nt
            return CriteriaResult.met(ACMGCriterion.PVS1, strength, evidence)

        # Other VCEPs publish gene-specific PVS1 trees whose critical-region /
        # codon-range gates, initiation-codon strengths and canonical-splice /
        # whole-gene-deletion calls deviate from the generic tree (and would
        # otherwise be under-called — e.g. single-exon genes like RAG1/GP9 where
        # NMD is never predicted, or last-exon truncations the generic tree
        # withholds). Like APC, the handler returns a final strength and so runs
        # BEFORE the generic LoF gate and the ClinVar-count strength caps. It
        # returns None for consequences the gene's VCEP does not special-case.
        from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1
        vcep = evaluate_vcep_pvs1(pc, self._splice_overrides)
        if vcep is not None:
            strength, evidence = vcep
            if strength == CriterionStrength.NOT_MET:
                return CriteriaResult.not_met(ACMGCriterion.PVS1, evidence)
            return CriteriaResult.met(ACMGCriterion.PVS1, strength, evidence)

        if pc.consequence not in _LOF_CONSEQUENCES:
            return CriteriaResult.not_met(ACMGCriterion.PVS1, "Not a LoF consequence")

        # VCEP gate: a VCEP that declined PVS1 for the gene (loss-of-function is
        # not the disease mechanism — MYOC, the RASopathy / cardiomyopathy panels,
        # the activating PIK3 genes, RYR1, VWF, …) withholds it even for a
        # bona-fide null variant.
        if self._spec.is_not_applicable(pc.gene_symbol):
            return CriteriaResult.not_met(
                ACMGCriterion.PVS1,
                f"{pc.gene_symbol}: VCEP designates PVS1 not applicable "
                "(LoF not the disease mechanism)",
            )

        # The ClinGen SVI 2019 decision tree is large enough (NMD prediction,
        # last-exon position, biologically-relevant transcript check, etc.)
        # that it lives in its own sub-package. Import is local to keep the
        # criteria layer free of the heavy pvs1 dependencies until needed.
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        # A VCEP that explicitly applies PVS1 has established LoF as the disease
        # mechanism — pass that through so the decision tree skips the ClinVar/
        # LOEUF heuristic (which can miss under-represented genes).
        lof_established = self._spec.is_applicable(pc.gene_symbol) or None
        strength, evidence = evaluate_pvs1(
            variant, annotation, self._cfg, lof_established=lof_established,
        )
        if strength == CriterionStrength.NOT_MET:
            return CriteriaResult.not_met(ACMGCriterion.PVS1, evidence)
        strength, evidence = self._gene_validity_downgrade(pc.gene_symbol, strength, evidence)
        if strength == CriterionStrength.NOT_MET:
            return CriteriaResult.not_met(ACMGCriterion.PVS1, evidence)
        return CriteriaResult.met(ACMGCriterion.PVS1, strength, evidence)

    def _gene_validity_downgrade(
        self, gene: str | None, strength: CriterionStrength, evidence: str,
    ) -> tuple[CriterionStrength, str]:
        """ClinGen SVI PVS1 recommendation (Abou Tayoun et al. 2018, Table 1): the
        decision tree applies at full strength only to genes with a Definitive or
        Strong gene-disease relationship. For other genes without a VCEP
        specification (ClinGen gene-disease validity Moderate or weaker, or not
        curated), the final strength is decreased by two levels (Very Strong ->
        Moderate, Strong -> Supporting; Moderate / Supporting -> not applied).
        The one-level option of Table 1 needs a null mouse model, which HUVar
        does not evaluate. Genes with a VCEP specification follow the VCEP."""
        if not getattr(self._cfg, "pvs1_gdv_downgrade", True) is True or not gene:
            return strength, evidence
        from acmg_classifier.classification.gene_validity import _vcep_genes
        if gene in _vcep_genes(getattr(self._cfg, "disease_prevalence_tsv", None)):
            return strength, evidence
        path = getattr(self._cfg, "clingen_gene_validity_csv", None)
        if not isinstance(path, Path) or not path.exists():
            return strength, evidence
        from acmg_classifier.local_db.clingen_gene_db import load_gene_validity
        classes = load_gene_validity(path).get(gene, [])
        if any(c in ("definitive", "strong") for c in classes):
            return strength, evidence
        gdv = ("ClinGen gene-disease validity " + "/".join(sorted(set(classes)))
               if classes else "no ClinGen gene-disease validity curation")
        new = _LEVEL_DOWN2.get(strength, CriterionStrength.NOT_MET)
        note = (f"{evidence}; decreased two levels from {strength.value} ({gdv}; "
                "full strength requires a Definitive/Strong gene-disease relationship, "
                "ClinGen SVI PVS1 Table 1)")
        return new, note


_LEVEL_DOWN2 = {
    CriterionStrength.VERY_STRONG: CriterionStrength.MODERATE,
    CriterionStrength.STRONG: CriterionStrength.SUPPORTING,
}
