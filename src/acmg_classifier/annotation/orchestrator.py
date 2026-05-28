"""
Annotation orchestrator: runs all local DB lookups in parallel and assembles AnnotationData.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed

import structlog

from acmg_classifier.config import Config
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.enums import SpliceTool
from acmg_classifier.models.variant import VariantRecord

log = structlog.get_logger()


class AnnotationOrchestrator:
    """Coordinates all local database queries for a batch of variants."""

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._vep = self._init_vep()
        self._gnomad = self._init_gnomad()
        self._clinvar_vcf_path = cfg.clinvar_vcf
        self._clinvar_sqlite_path = cfg.clinvar_sqlite
        self._am_path = cfg.alphamissense_tsv
        self._repeat_path = cfg.repeatmasker_bed
        self._splice = self._init_splice()

    def _init_vep(self):
        from acmg_classifier.local_db.vep_runner import LocalVEPRunner
        from acmg_classifier.setup.vep_installer import find_vep_cmd
        vep_cmd = find_vep_cmd()
        return LocalVEPRunner(
            vep_cmd=vep_cmd,
            cache_dir=self._cfg.vep_cache_dir,
            fasta=self._cfg.genome_fasta,
            assembly=self._cfg.assembly.value,
            workers=self._cfg.workers,
        )

    def _init_gnomad(self):
        from acmg_classifier.local_db.gnomad_db import GnomADDB
        return GnomADDB(self._cfg.gnomad_duckdb, self._cfg.gnomad_constraint_tsv)

    def _init_splice(self):
        if self._cfg.splice_tool == SpliceTool.SPLICEAI:
            from acmg_classifier.local_db.splice.spliceai_predictor import SpliceAIPredictor
            return SpliceAIPredictor(self._cfg.spliceai_vcf, self._cfg.spliceai_indel_vcf)
        from acmg_classifier.local_db.splice.squirls_predictor import SquirlsPredictor
        return SquirlsPredictor(self._cfg.squirls_db_dir)

    def annotate_batch(
        self,
        variants: list[VariantRecord],
    ) -> dict[str, AnnotationData]:
        # VEP runs in a single batch subprocess
        vep_results = self._vep.annotate_batch(variants, batch_size=self._cfg.vep_batch_size)

        results: dict[str, AnnotationData] = {}

        def _annotate_single(v: VariantRecord) -> tuple[str, AnnotationData]:
            return v.key, self._annotate_one(v, vep_results.get(v.key, []))

        with ThreadPoolExecutor(max_workers=self._cfg.workers) as pool:
            futures = {pool.submit(_annotate_single, v): v for v in variants}
            for future in as_completed(futures):
                try:
                    key, ann = future.result()
                    results[key] = ann
                except Exception as exc:
                    v = futures[future]
                    log.error("annotation_failed", variant=v.key, error=str(exc))

        return results

    def _annotate_one(self, variant: VariantRecord, consequences) -> AnnotationData:
        from acmg_classifier.local_db.clinvar_vcf import query_clinvar_vcf
        from acmg_classifier.local_db.alphamissense_db import query_alphamissense
        from acmg_classifier.local_db.repeatmasker_db import query_repeat

        gnomad = self._gnomad.query(variant.chrom, variant.pos, variant.ref, variant.alt)

        primary = None
        if consequences:
            primary = consequences[0]

        if gnomad and primary:
            gnomad = self._gnomad.enrich_with_constraint(gnomad, primary.gene_symbol)

        clinvar_vcf_recs = query_clinvar_vcf(
            self._clinvar_vcf_path,
            variant.chrom, variant.pos, variant.ref, variant.alt,
        )

        alphamissense = query_alphamissense(
            self._am_path,
            variant.chrom, variant.pos, variant.ref, variant.alt,
        )

        splice = self._splice.predict(variant)

        repeat = query_repeat(self._repeat_path, variant.chrom, variant.pos)

        return AnnotationData(
            consequences=consequences,
            gnomad=gnomad,
            alphamissense=alphamissense,
            splice=splice if splice.is_available else None,
            clinvar_vcf=clinvar_vcf_recs,
            repeat=repeat,
        )
