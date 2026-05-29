from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, computed_field

from acmg_classifier.models.enums import ACMGCriterion, Pathogenicity
from acmg_classifier.models.criteria import CriteriaResult


class ClassificationResult(BaseModel):
    """Final classification for a single variant."""

    model_config = {"arbitrary_types_allowed": True}

    variant_id: str
    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    filter: Optional[str] = None  # VCF FILTER column (e.g. PASS, ., LowQual)
    transcript_id: Optional[str] = None
    gene_symbol: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None

    annotation: Optional[object] = None  # AnnotationData (avoid circular import at runtime)

    criteria_results: list[CriteriaResult]

    # Rule-based (ACMG 2015)
    classification_2015: Pathogenicity
    classification_2015_rules: str  # e.g. "PVS1 + PM2"

    # Bayesian (Tavtigian 2020 + Bergquist 2024)
    bayesian_score: int
    classification_bayesian: Pathogenicity

    # Warnings (e.g. SQUIRLS thresholds not Walker-calibrated)
    warnings: list[str] = []

    @computed_field  # type: ignore[prop-decorator]
    @property
    def triggered_criteria(self) -> list[ACMGCriterion]:
        # The classifier emits one CriteriaResult per criterion (including
        # not_met) so the full audit trail is preserved. Reports usually only
        # need the "active" subset — triggered AND not suppressed.
        return [r.criterion for r in self.criteria_results if r.triggered and not r.suppressed]

    @computed_field  # type: ignore[prop-decorator]
    @property
    def suppressed_criteria(self) -> list[ACMGCriterion]:
        # Exposed separately so reviewers can see *why* a criterion was
        # excluded (typically to avoid double-counting evidence already
        # captured by a higher-strength criterion, e.g. PP3 under PVS1).
        return [r.criterion for r in self.criteria_results if r.suppressed]
