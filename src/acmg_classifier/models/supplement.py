from __future__ import annotations
from typing import Optional
from pydantic import BaseModel

from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength


class SupplementEntry(BaseModel):
    """One row from the manual evidence supplement TSV."""

    variant_id: str        # CHROM:POS:REF:ALT
    criterion: ACMGCriterion
    strength: CriterionStrength
    evidence: str          # free-text PMID or description
    source: str = "supplement"

    # For criteria that accept variable strength (e.g. PS1, PVS1_Strong)
    override_strength: Optional[CriterionStrength] = None
