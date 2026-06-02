"""
Annotation orchestrator: runs all local DB lookups in parallel and assembles AnnotationData.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed

import structlog

from acmg_classifier.config import Config
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.enums import ConsequenceType, SpliceTool
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.utils.progress import progress_bar

log = structlog.get_logger()


def _resolve_uniprot_id(primary, consequences):
    """Resolve the UniProt accession used as the ESM1b lookup key.

    VEP's --merged cache only attaches swissprot/trembl xrefs to Ensembl
    (ENST) transcripts, not RefSeq (NM_). Because `primary_consequence`
    prefers the RefSeq MANE transcript, `primary.uniprot_id` is frequently
    None for perfectly scorable missense variants (e.g. TP53 R248W). Fall
    back to a sibling transcript describing the *same* substitution — same
    protein_position and amino_acids — so the variant can still be scored.

    Returns the primary's own accession when present, otherwise the first
    matching sibling's, or None when nothing carries a UniProt xref.
    """
    if primary is None:
        return None
    if primary.uniprot_id:
        return primary.uniprot_id
    for c in consequences or ():
        if c is primary or not c.uniprot_id:
            continue
        if c.consequence != ConsequenceType.MISSENSE:
            continue
        # Same residue position and same WT/alt pair guarantees the sibling
        # describes the identical protein substitution, so its UniProt
        # accession (and shared residue numbering) is safe to reuse.
        if c.protein_position != primary.protein_position:
            continue
        if c.amino_acids != primary.amino_acids:
            continue
        return c.uniprot_id
    return None


class AnnotationOrchestrator:
    """Coordinates all local database queries for a batch of variants."""

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._vep = self._init_vep()
        self._gnomad = self._init_gnomad()
        self._clinvar_vcf_path = cfg.clinvar_vcf
        self._clinvar_sqlite_path = cfg.clinvar_sqlite
        self._am_path = cfg.alphamissense_tsv
        self._esm1b_path = cfg.esm1b_sqlite
        self._repeat_path = cfg.repeatmasker_bed
        self._splice = self._init_splice()
        # SQUIRLS predictor for secondary reporting — always initialized so
        # squirls_score is available in the TSV regardless of which splice tool
        # is configured as primary. Reuses the same instance when SQUIRLS IS
        # the primary tool to avoid double DB initialization.
        if cfg.splice_tool == SpliceTool.SQUIRLS:
            self._squirls = self._splice
        else:
            from acmg_classifier.local_db.splice.squirls_predictor import SquirlsPredictor
            self._squirls = SquirlsPredictor(cfg.squirls_db_dir)

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
        if self._cfg.splice_tool == SpliceTool.SQUIRLS:
            # Retained for when the SQUIRLS DB is downloadable again. Not the
            # default and not currently selectable from the CLI.
            from acmg_classifier.local_db.splice.squirls_predictor import SquirlsPredictor
            return SquirlsPredictor(self._cfg.squirls_db_dir)
        # MMSplice DISABLED (dependency conflict). Code retained, commented out:
        # if self._cfg.splice_tool == SpliceTool.MMSPLICE:
        #     from acmg_classifier.local_db.splice.mmsplice_predictor import MMSplicePredictor
        #     return MMSplicePredictor(self._cfg.mmsplice_gtf, self._cfg.genome_fasta)
        # Default (SpliceTool.NONE): splice evaluation is disabled — no splice
        # predictor, so no splice-based evidence is contributed.
        from acmg_classifier.local_db.splice.base import NullSplicePredictor
        return NullSplicePredictor()

    def annotate_batch(
        self,
        variants: list[VariantRecord],
    ) -> dict[str, AnnotationData]:
        """Annotate a batch: one VEP subprocess + parallel per-variant DB lookups.

        VEP is run once for the whole batch (subprocess startup cost is the
        bottleneck — see vep_runner.annotate_batch). The remaining lookups
        (gnomAD/ClinVar/AlphaMissense/etc.) are I/O-bound and embarrassingly
        parallel per variant, so a thread pool is the right model: GIL is
        released during SQLite/tabix calls, and we want to overlap network/
        disk waits rather than CPU work."""
        vep_results = self._vep.annotate_batch(variants, batch_size=self._cfg.vep_batch_size)

        # MMSplice batch precompute DISABLED (MMSplice integration is off).
        # The tabix-backed predictors (SpliceAI/SQUIRLS) look up precomputed
        # scores per variant and need no batch step. Re-enable with MMSplice:
        # self._splice.precompute(variants)

        results: dict[str, AnnotationData] = {}

        def _annotate_single(v: VariantRecord) -> tuple[str, AnnotationData]:
            # vep_results.get(..., []) gives an empty list when VEP dropped
            # the variant (rare but possible — see vep_runner under-count
            # log); the rest of the annotation still proceeds.
            return v.key, self._annotate_one(v, vep_results.get(v.key, []))

        with ThreadPoolExecutor(max_workers=self._cfg.workers) as pool:
            futures = {pool.submit(_annotate_single, v): v for v in variants}
            # Wrap as_completed with progress_bar so the user sees one tick
            # per completed variant rather than a long silent wait. The
            # context manager is a no-op when stderr is not a tty.
            with progress_bar("Annotating variants", total=len(variants)) as advance:
                for future in as_completed(futures):
                    try:
                        key, ann = future.result()
                        results[key] = ann
                    except Exception as exc:
                        # Per-variant failures are isolated: log and continue so
                        # a single bad record cannot abort the whole batch.
                        v = futures[future]
                        log.error("annotation_failed", variant=v.key, error=str(exc))
                    finally:
                        advance()

        return results

    def _annotate_one(self, variant: VariantRecord, consequences) -> AnnotationData:
        """Build the AnnotationData for a single variant by querying every DB.

        Lookups are sequential here (one variant at a time) because the
        thread pool in annotate_batch already provides the parallelism
        across variants. Doing nested parallelism would oversubscribe the
        DB clients without throughput gain."""
        from acmg_classifier.local_db.clinvar_vcf import query_clinvar_vcf
        from acmg_classifier.local_db.alphamissense_db import query_alphamissense
        from acmg_classifier.local_db.repeatmasker_db import query_repeat

        gnomad = self._gnomad.query(variant.chrom, variant.pos, variant.ref, variant.alt)

        # `consequences` is pre-sorted by vep_runner._parse_vep_record so the
        # first entry is already the clinically-preferred (MANE > canonical)
        # transcript — no need to re-pick primary here.
        primary = None
        if consequences:
            primary = consequences[0]

        # pLI/LOEUF/missense-Z are gene-level (not variant-level) and only
        # meaningful once we know which gene the variant belongs to. Enrich
        # *after* the primary consequence is resolved.
        if gnomad and primary:
            gnomad = self._gnomad.enrich_with_constraint(gnomad, primary.gene_symbol)

        clinvar_vcf_recs = query_clinvar_vcf(
            self._clinvar_vcf_path,
            variant.chrom, variant.pos, variant.ref, variant.alt,
        )

        # Both missense predictors are always queried so the TSV always carries
        # both scores for manual review. Only the tool selected by insilico_tool
        # config is used inside the ACMG criteria (PP3/BP4) to avoid inflating
        # evidence by combining tools that share training data.
        alphamissense = query_alphamissense(
            self._am_path,
            variant.chrom, variant.pos, variant.ref, variant.alt,
        )
        esm1b = self._lookup_esm1b(primary, consequences)

        splice = self._splice.predict(variant)
        # When SQUIRLS is the primary splice tool, reuse the same result.
        # When SpliceAI is primary, run SQUIRLS separately for TSV reporting.
        squirls = splice if self._squirls is self._splice else self._squirls.predict(variant)

        repeat = query_repeat(self._repeat_path, variant.chrom, variant.pos)

        return AnnotationData(
            consequences=consequences,
            gnomad=gnomad,
            alphamissense=alphamissense,
            esm1b=esm1b,
            # Dropping the splice record entirely when the predictor was
            # unavailable lets downstream criteria treat "no data" as
            # "splice unknown" rather than "splice = 0".
            splice=splice if splice.is_available else None,
            squirls=squirls if squirls.is_available else None,
            clinvar_vcf=clinvar_vcf_recs,
            repeat=repeat,
        )

    def _lookup_esm1b(self, primary, consequences):
        """Lookup ESM1b LLR using UniProt accession + protein position + alt AA.

        Brandes 2023 archives are keyed by UniProt (SwissProt/TrEMBL), so we
        depend on VEP --uniprot. The UniProt accession is resolved by
        `_resolve_uniprot_id`, which falls back to a sibling Ensembl
        transcript when the primary (RefSeq) one carries no UniProt xref.
        """
        # Cascade of gates — each `return None` represents a precondition the
        # ESM1b archive needs. We bail early rather than fall through so the
        # actual DB lookup is reached only when every required field is
        # present and well-formed:
        #   - missense only (Brandes archive is missense-substitution scored)
        #   - protein_position + alt amino acid identify the substitution
        #   - alt_aa must be a single 1-letter code (Brandes uses 1-letter)
        #   - UniProt ID is the file-name key in the archive
        if primary is None:
            return None
        if primary.consequence != ConsequenceType.MISSENSE:
            return None
        if primary.protein_position is None or not primary.amino_acids:
            return None
        parts = primary.amino_acids.split("/")
        if len(parts) != 2:
            return None
        alt_aa = parts[1].strip()
        if not alt_aa or len(alt_aa) != 1:
            return None
        uniprot_id = _resolve_uniprot_id(primary, consequences)
        if not uniprot_id:
            return None

        from acmg_classifier.local_db.esm1b_db import query_esm1b
        return query_esm1b(
            self._esm1b_path,
            uniprot_id,
            primary.protein_position,
            alt_aa,
        )
