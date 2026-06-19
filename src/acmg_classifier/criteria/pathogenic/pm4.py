"""PM4 -- protein length change due to in-frame indel or stop-loss."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


def _load_pm4_columns(tsv_path) -> tuple[frozenset[str], dict[str, int]]:
    """Load PM4 per-gene config from ``disease_prevalence.tsv``:

    * the set of genes whose VCEP declined PM4 (``pm4`` == ``not_applicable``);
    * ``pm4_supporting_max_aa`` — the in-frame-indel size (amino acids) at or
      below which PM4 downgrades from the default Moderate to Supporting.

    A missing file/column degrades to "no gene declined" / "no size downgrade"."""
    import csv
    declined: set[str] = set()
    max_aa: dict[str, int] = {}
    try:
        if not tsv_path.exists():
            return frozenset(), {}
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                if (row.get("pm4") or "").strip().lower() == "not_applicable":
                    declined.add(gene)
                raw = (row.get("pm4_supporting_max_aa") or "").strip()
                if raw:
                    try:
                        max_aa[gene] = int(raw)
                    except ValueError:
                        pass
    except OSError:
        return frozenset(), {}
    return frozenset(declined), max_aa


class PM4Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._not_applicable, self._supporting_max_aa = _load_pm4_columns(
            cfg.disease_prevalence_tsv
        )

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.PM4, "No consequence")

        # VCEP gate: cancer panels (BRCA1/2, MMR, TP53, APC, PALB2) and the
        # PI3K-pathway specs decline PM4 — withhold it for those genes.
        if pc.gene_symbol in self._not_applicable:
            return CriteriaResult.not_met(
                ACMGCriterion.PM4, f"{pc.gene_symbol}: VCEP designates PM4 not applicable"
            )

        # PM4 fires for changes that alter protein length without truncation
        # (frameshift/stop-gain are PVS1). For in-frame indels inside a
        # repetitive region, the length change is expected/tolerated and the
        # variant is better captured by BP3 — explicitly redirect rather than
        # silently failing so the audit trail reflects the reasoning.
        if pc.consequence in (
            ConsequenceType.INFRAME_INSERTION,
            ConsequenceType.INFRAME_DELETION,
            ConsequenceType.STOP_LOST,
        ):
            if annotation.repeat and annotation.repeat.in_repeat:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM4,
                    f"In-frame indel in repeat ({annotation.repeat.repeat_class}); use BP3",
                )
            # Size-based downgrade: a VCEP may apply PM4 at Supporting for a small
            # in-frame indel (e.g. a single amino acid). Stop-loss is not size-
            # scoped, so it keeps the default Moderate. The amino-acid change is
            # the net codon-length difference of an in-frame indel.
            max_aa = self._supporting_max_aa.get(pc.gene_symbol)
            if (
                max_aa is not None
                and pc.consequence in (
                    ConsequenceType.INFRAME_INSERTION,
                    ConsequenceType.INFRAME_DELETION,
                )
                and _indel_aa_change(variant) is not None
                and _indel_aa_change(variant) <= max_aa
            ):
                return CriteriaResult.met(
                    ACMGCriterion.PM4,
                    strength=CriterionStrength.SUPPORTING,
                    evidence=(
                        f"{pc.consequence.value} of <={max_aa} aa "
                        f"({pc.gene_symbol}: VCEP PM4_Supporting)"
                    ),
                )
            return CriteriaResult.met(
                ACMGCriterion.PM4,
                evidence=f"{pc.consequence.value} outside repeat region",
            )
        return CriteriaResult.not_met(ACMGCriterion.PM4, "Not an in-frame indel or stop-loss")


def _indel_aa_change(variant: VariantRecord) -> int | None:
    """Net amino-acid length change of an in-frame indel = |len(ref)-len(alt)|/3.
    Returns ``None`` when the length difference is not a whole number of codons
    (defensive — an in-frame indel should always be a multiple of 3)."""
    diff = abs(len(variant.ref) - len(variant.alt))
    if diff == 0 or diff % 3 != 0:
        return None
    return diff // 3
