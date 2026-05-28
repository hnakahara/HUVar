"""SpliceAI predictor (optional; requires Illumina commercial license)."""
from __future__ import annotations
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.local_db.splice.base import SplicePredictor
from acmg_classifier.models.annotation import SpliceScore
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.utils.tabix import open_tabix, fetch_region

log = structlog.get_logger()


class SpliceAIPredictor(SplicePredictor):
    """
    Looks up precomputed SpliceAI scores from tabix-indexed VCF.

    Requires Illumina commercial license for the score file.
    Walker 2023 thresholds: >=0.2 -> PP3 Moderate, <=0.1 -> BP4 Supporting.
    Accepts separate SNV and indel VCF files; selects the appropriate file per variant.
    """

    def __init__(self, snv_vcf: Optional[Path], indel_vcf: Optional[Path] = None) -> None:
        self._snv_vcf = snv_vcf
        self._indel_vcf = indel_vcf

    def is_available(self) -> bool:
        snv_ok = self._snv_vcf is not None and self._snv_vcf.exists()
        indel_ok = self._indel_vcf is not None and self._indel_vcf.exists()
        return snv_ok or indel_ok

    def _vcf_for(self, variant: VariantRecord) -> Optional[Path]:
        if len(variant.ref) == 1 and len(variant.alt) == 1:
            return self._snv_vcf
        return self._indel_vcf

    def predict(self, variant: VariantRecord) -> SpliceScore:
        vcf_path = self._vcf_for(variant)
        if vcf_path is None or not vcf_path.exists():
            return SpliceScore(tool="spliceai", is_available=False)

        try:
            with open_tabix(vcf_path) as tf:
                for line in fetch_region(tf, variant.chrom, variant.pos, variant.pos):
                    fields = line.split("\t")
                    if len(fields) < 8:
                        continue
                    if fields[1] != str(variant.pos) or fields[3] != variant.ref or fields[4] != variant.alt:
                        continue
                    info = fields[7]
                    max_delta = _extract_spliceai_max_delta(info)
                    if max_delta is not None:
                        return SpliceScore(
                            tool="spliceai",
                            is_available=True,
                            max_delta=max_delta,
                            raw_score=max_delta,
                        )
        except Exception as exc:
            log.error("spliceai_error", error=str(exc))
        return SpliceScore(tool="spliceai", is_available=True, raw_score=None)


def _extract_spliceai_max_delta(info: str) -> Optional[float]:
    for token in info.split(";"):
        if token.startswith("SpliceAI="):
            payload = token[len("SpliceAI="):]
            parts = payload.split("|")
            try:
                scores = [float(parts[i]) for i in [2, 3, 4, 5] if i < len(parts)]
                return max(scores) if scores else None
            except ValueError:
                pass
    return None
