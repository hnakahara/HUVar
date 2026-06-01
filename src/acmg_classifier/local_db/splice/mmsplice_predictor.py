"""MMSplice runtime splice predictor — CURRENTLY DISABLED.

This module is intentionally left in the tree but is NOT wired into the
pipeline: the orchestrator's MMSplice branch, the SpliceTool.MMSPLICE enum
value, the criteria branches, the setup-data step, and the pyproject extra are
all commented out. Reason: mmsplice's dependency chain (numpy<2, cyvcf2<=0.30.x,
pyranges 0.0.x) conflicts with this project's cyvcf2/numpy. Re-enable by
uncommenting those references once the dependency conflict is resolved.

----------------------------------------------------------------------------
MMSplice runtime splice predictor (open-source default).

MMSplice (gagneurlab/MMSplice_MTSplice) computes splice-effect scores at runtime
from a gene-annotation GTF, a reference FASTA, and the batch VCF. Unlike the
tabix-backed SpliceAI/SQUIRLS predictors there is no precomputed score file, so
the whole batch is scored in a single model pass via precompute() and cached;
predict() then reads the cache.

The `mmsplice` package (TensorFlow/Keras) is an OPTIONAL dependency — when it is
not installed the predictor reports is_available()=False and the pipeline
silently skips splice evidence rather than crashing.

Score convention: delta_logit_psi is on a logit scale (MMSplice, Genome Biology
2019). |delta_logit_psi| > 2 is a strong splice effect; the sign indicates the
direction (negative = increased exon exclusion). We store it in
SpliceScore.raw_score with tool="mmsplice" and let the ACMG criteria interpret
the magnitude.
"""
from __future__ import annotations
import tempfile
from pathlib import Path

import structlog

from acmg_classifier.local_db.splice.base import SplicePredictor
from acmg_classifier.models.annotation import SpliceScore
from acmg_classifier.models.variant import VariantRecord

log = structlog.get_logger()


class MMSplicePredictor(SplicePredictor):
    """Runtime MMSplice scoring with batch precompute + per-variant lookup."""

    def __init__(self, gtf_path: Path, fasta_path: Path) -> None:
        self._gtf = gtf_path
        self._fasta = fasta_path
        # key (CHROM:POS:REF:ALT) -> delta_logit_psi. Populated by precompute().
        self._cache: dict[str, float] = {}

    def is_available(self) -> bool:
        """True only when the GTF/FASTA exist AND mmsplice is importable.

        The import check is what makes mmsplice an optional dependency: a user
        running AlphaMissense/SpliceAI never needs TensorFlow installed."""
        if not self._gtf.exists() or not self._fasta.exists():
            return False
        try:
            import mmsplice  # noqa: F401
        except Exception:
            return False
        return True

    def precompute(self, variants: list[VariantRecord]) -> None:
        """Score the entire batch in one MMSplice pass and cache the results.

        Called by the orchestrator on the MAIN thread (before the per-variant
        thread pool starts) because TensorFlow/Keras model inference is not
        thread-safe. Failures are logged and swallowed so a splice-tool problem
        cannot abort the whole annotation run — predict() then returns empty
        scores and the splice-dependent criteria treat them as "no data"."""
        if not self.is_available():
            log.warning(
                "mmsplice_unavailable",
                gtf=str(self._gtf),
                fasta=str(self._fasta),
                hint="pip install -e .[mmsplice] and run setup to fetch the GTF",
            )
            return
        annotatable = [v for v in variants if v.alt and v.alt != "."]
        if not annotatable:
            return
        try:
            self._run(annotatable)
        except Exception as exc:
            log.error("mmsplice_precompute_failed", error=str(exc))

    def _run(self, variants: list[VariantRecord]) -> None:
        import pysam
        from mmsplice.vcf_dataloader import SplicingVCFDataloader
        from mmsplice import MMSplice, predict_all_table
        from mmsplice.utils import max_varEff

        with tempfile.TemporaryDirectory() as tmp:
            vcf_path = Path(tmp) / "input.vcf"
            self._write_vcf(variants, vcf_path)
            # MMSplice's cyvcf2-based dataloader expects a bgzipped, tabix-indexed
            # VCF (the upstream examples use *.vcf.gz). Compress + index here.
            gz_path = Path(tmp) / "input.vcf.gz"
            pysam.tabix_compress(str(vcf_path), str(gz_path), force=True)
            pysam.tabix_index(str(gz_path), preset="vcf", force=True)

            dl = SplicingVCFDataloader(str(self._gtf), str(self._fasta), str(gz_path))
            model = MMSplice()
            df = predict_all_table(model, dl)
            # One variant can map to several exons; collapse to the maximum
            # absolute effect per variant (the upstream-recommended summary).
            df_max = max_varEff(df)
            self._ingest(df_max)
        log.info("mmsplice_precomputed", scored=len(self._cache), input=len(variants))

    def _ingest(self, df_max) -> None:
        """Pull (ID -> delta_logit_psi) out of the max_varEff dataframe.

        `ID` is the VCF ID column we set to variant.key; depending on the
        mmsplice version it is either a column or the dataframe index."""
        if "ID" in getattr(df_max, "columns", []):
            ids = df_max["ID"].tolist()
        else:
            ids = df_max.index.tolist()
        scores = df_max["delta_logit_psi"].tolist()
        for key, score in zip(ids, scores):
            try:
                self._cache[str(key)] = float(score)
            except (TypeError, ValueError):
                continue

    def _write_vcf(self, variants: list[VariantRecord], path: Path) -> None:
        """Write a minimal, coordinate-sorted VCF for MMSplice.

        CHROM is written WITHOUT the 'chr' prefix to match the Ensembl GTF/FASTA
        naming (same convention as vep_runner._write_input_vcf). The ID column
        carries variant.key so the prediction rows can be joined back."""
        rows = sorted(variants, key=lambda v: (v.chrom, v.pos))
        with path.open("w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for v in rows:
                chrom = v.chrom[3:] if v.chrom.startswith("chr") else v.chrom
                fh.write(f"{chrom}\t{v.pos}\t{v.key}\t{v.ref}\t{v.alt}\t.\t.\t.\n")

    def predict(self, variant: VariantRecord) -> SpliceScore:
        """Per-variant lookup into the precompute() cache (thread-safe read).

        Returns raw_score=None when the variant was not in an exon-adjacent
        region MMSplice scores, or when precompute() was skipped/failed."""
        if not self.is_available():
            return SpliceScore(tool="mmsplice", is_available=False)
        return SpliceScore(
            tool="mmsplice",
            is_available=True,
            raw_score=self._cache.get(variant.key),
        )
