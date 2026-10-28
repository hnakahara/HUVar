"""Models for somatic oncogenicity classification (Horak et al., Genet Med 2022;
24:986-998, ClinGen/CGC/VICC SOP; erratum Genet Med 2022;24:1991)."""
from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, computed_field


class OncCriterion(str, Enum):
    OVS1 = "OVS1"
    OS1 = "OS1"
    OS2 = "OS2"
    OS3 = "OS3"
    OM1 = "OM1"
    OM2 = "OM2"
    OM3 = "OM3"
    OM4 = "OM4"
    OP1 = "OP1"
    OP2 = "OP2"
    OP3 = "OP3"
    OP4 = "OP4"
    SBVS1 = "SBVS1"
    SBS1 = "SBS1"
    SBS2 = "SBS2"
    SBP1 = "SBP1"
    SBP2 = "SBP2"


# SOP Table 1 point values.
POINTS: dict[OncCriterion, int] = {
    OncCriterion.OVS1: 8,
    OncCriterion.OS1: 4, OncCriterion.OS2: 4, OncCriterion.OS3: 4,
    OncCriterion.OM1: 2, OncCriterion.OM2: 2, OncCriterion.OM3: 2, OncCriterion.OM4: 2,
    OncCriterion.OP1: 1, OncCriterion.OP2: 1, OncCriterion.OP3: 1, OncCriterion.OP4: 1,
    OncCriterion.SBVS1: -8,
    OncCriterion.SBS1: -4, OncCriterion.SBS2: -4,
    OncCriterion.SBP1: -1, OncCriterion.SBP2: -1,
}

# Criteria that need information HUVar does not evaluate automatically; they are
# applied only from manual evidence (--evidence / --supplement).
MANUAL_ONLY = {OncCriterion.OP2}


class Oncogenicity(str, Enum):
    ONCOGENIC = "Oncogenic"
    LIKELY_ONCOGENIC = "LikelyOncogenic"
    VUS = "VUS"
    LIKELY_BENIGN = "LikelyBenign"
    BENIGN = "Benign"


def classify_points(score: int) -> Oncogenicity:
    """SOP Table 1: Oncogenic >= 10, Likely oncogenic 6-9, VUS 0-5,
    Likely benign -1 to -6, Benign <= -7."""
    if score >= 10:
        return Oncogenicity.ONCOGENIC
    if score >= 6:
        return Oncogenicity.LIKELY_ONCOGENIC
    if score >= 0:
        return Oncogenicity.VUS
    if score >= -6:
        return Oncogenicity.LIKELY_BENIGN
    return Oncogenicity.BENIGN


class OncCriterionResult(BaseModel):
    criterion: OncCriterion
    met: bool = False
    evidence: str = ""
    suppressed: bool = False
    manual: bool = False

    @property
    def points(self) -> int:
        return POINTS[self.criterion] if self.met and not self.suppressed else 0

    @property
    def active(self) -> bool:
        return self.met and not self.suppressed


class ManualOncEvidence(BaseModel):
    """Curator evidence for one somatic criterion (supplement row / --evidence)."""

    variant_id: str
    criterion: OncCriterion
    met: bool = True
    evidence: str = ""


class OncogenicityResult(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    variant_id: str
    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    filter: Optional[str] = None
    transcript_id: Optional[str] = None
    gene_symbol: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    consequence: Optional[str] = None
    gene_role: str = "not_cancer_gene"
    gene_role_source: str = ""
    gene_role_note: str = ""
    hotspot: str = ""
    criteria: list[OncCriterionResult]
    score: int
    classification: Oncogenicity
    annotation: Optional[object] = None
    warnings: list[str] = []

    @computed_field  # type: ignore[prop-decorator]
    @property
    def met_criteria(self) -> list[str]:
        return [r.criterion.value for r in self.criteria if r.active]

    def summary(self) -> str:
        """'OS3 (+4), OM1 (+2), OP4 (+1)'."""
        return ", ".join(f"{r.criterion.value} ({r.points:+d})" for r in self.criteria if r.active)
