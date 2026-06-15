"""Main pipeline: reads VCF → annotates → classifies → writes output."""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import structlog

from acmg_classifier.config import Config
from acmg_classifier.io.vcf_reader import read_vcf, detect_assembly_from_header
from acmg_classifier.io.supplement_reader import read_supplement
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.classification import ClassificationResult
from acmg_classifier.utils.progress import track

log = structlog.get_logger()


def run_pipeline(
    vcf_path: Path,
    cfg: Config,
    output_path: Optional[Path] = None,
    supplement_path: Optional[Path] = None,
    limit: Optional[int] = None,
    profile_path: Optional[Path] = None,
) -> list[ClassificationResult]:
    """Classify all variants in *vcf_path* and write TSV to *output_path*.

    `limit` truncates the variant list before annotation — used together
    with `profile_path` for fast performance profiling on a representative
    subset rather than waiting hours for a full run.

    `profile_path` enables cProfile around the annotation + classification
    phases and dumps the binary stats to that path; top functions are also
    written to the log so the user can spot hot paths without re-loading.

    Top-level orchestration. The deferred imports keep startup fast for
    CLI subcommands that don't classify (validate/status) by avoiding
    eager imports of heavy dependencies (pysam, duckdb, VEP cache loaders).
    """
    from acmg_classifier.annotation.orchestrator import AnnotationOrchestrator
    from acmg_classifier.criteria.registry import CriteriaRegistry
    from acmg_classifier.classification.classifier_2015 import Classifier2015
    from acmg_classifier.classification.classifier_bayesian import ClassifierBayesian
    from acmg_classifier.io.tsv_writer import write_tsv, write_skipped

    # Auto-detect assembly from ##contig lines if the user didn't pass it.
    # Most clinical VCFs encode GRCh37/GRCh38 in the header, so we don't
    # force the user to repeat it on the CLI.
    if cfg.assembly is None:
        detected = detect_assembly_from_header(vcf_path)
        if detected:
            cfg = cfg.model_copy(update={"assembly": detected})

    supplement = {}
    if supplement_path:
        supplement = read_supplement(supplement_path)

    orchestrator = AnnotationOrchestrator(cfg)
    registry = CriteriaRegistry(cfg)
    clf_2015 = Classifier2015()
    clf_bay = ClassifierBayesian()

    from acmg_classifier.models.annotation import AnnotationData

    results: list[ClassificationResult] = []
    all_records = list(read_vcf(vcf_path, cfg.assembly))
    # ALT='.' records carry no alternate allele — they are not variants and
    # cannot be classified. Exclude them from the main output and write them
    # to a `.not_annotated` sidecar so users can see *which* rows were
    # skipped (silent dropping would surprise reviewers).
    variants = [v for v in all_records if v.alt and v.alt != "."]
    skipped = [v for v in all_records if not (v.alt and v.alt != ".")]
    log.info(
        "variants_loaded",
        total=len(all_records), classifiable=len(variants), skipped_no_alt=len(skipped),
    )

    # --limit truncates AFTER skip filtering so the user gets the requested
    # number of *classifiable* variants in the profile, not N rows that
    # might be mostly ALT='.' placeholders.
    if limit is not None and limit < len(variants):
        log.info("variants_limited", from_count=len(variants), to_count=limit)
        variants = variants[:limit]

    # Optional cProfile wrap. Both annotation (subprocess + I/O heavy) and
    # the per-variant classification loop run inside the profile so the
    # caller can see whether the bottleneck is annotation or evaluation.
    profiler = None
    if profile_path is not None:
        import cProfile
        profiler = cProfile.Profile()
        profiler.enable()

    annotations = orchestrator.annotate_batch(variants)
    log.info("annotations_returned", count=len(annotations))

    missing_keys: list[str] = []
    annotated_count = 0

    # Classification is CPU-cheap relative to annotation but still benefits
    # from a progress bar when there are many variants — without it the user
    # only sees the annotation bar finish and then a silent gap until the
    # final summary line.
    for variant in track(variants, "Classifying", total=len(variants)):
        ann = annotations.get(variant.key)
        # Annotation failures are non-fatal — we emit a row with an empty
        # AnnotationData and surface the warning in the output so the user
        # sees there *was* a variant rather than a silent gap. This matches
        # clinical-lab expectations of "explain every input row".
        if ann is None:
            missing_keys.append(variant.key)
            ann = AnnotationData()
            extra_warnings = ["VEP returned no annotation for this variant"]
        else:
            annotated_count += 1
            extra_warnings = []

        sup_entries = supplement.get(variant.key, [])
        criteria_results = registry.evaluate_all(variant, ann, sup_entries)
        classification_2015, rules = clf_2015.classify(criteria_results)
        score, classification_bay = clf_bay.classify(criteria_results)

        pc = ann.primary_consequence
        result = ClassificationResult(
            variant_id=variant.key,
            chrom=variant.chrom,
            pos=variant.pos,
            ref=variant.ref,
            alt=variant.alt,
            filter=variant.filter,
            transcript_id=pc.transcript_id if pc else None,
            gene_symbol=pc.gene_symbol if pc else None,
            hgvs_c=pc.hgvs_c if pc else None,
            hgvs_p=pc.hgvs_p if pc else None,
            annotation=ann,
            criteria_results=criteria_results,
            classification_2015=classification_2015,
            classification_2015_rules=rules,
            bayesian_score=score,
            classification_bayesian=classification_bay,
            warnings=extra_warnings,
        )
        results.append(result)

    # Stop the profiler BEFORE the TSV write so output I/O doesn't pollute
    # the top-N. Dump the raw stats for offline `snakeviz` / `pstats` use
    # and emit the top 30 cumulative-time entries to stderr immediately so
    # users can see the bottleneck without leaving the terminal.
    if profiler is not None and profile_path is not None:
        import io as _io
        import pstats
        profiler.disable()
        profile_path.parent.mkdir(parents=True, exist_ok=True)
        profiler.dump_stats(str(profile_path))
        buf = _io.StringIO()
        pstats.Stats(profiler, stream=buf).sort_stats("cumulative").print_stats(30)
        log.info("profile_dumped", path=str(profile_path))
        # Use sys.stderr directly: structlog would re-serialise the multi-
        # line stats output as one giant field, which is unreadable.
        import sys as _sys
        _sys.stderr.write("\n===== cProfile top 30 (cumulative) =====\n")
        _sys.stderr.write(buf.getvalue())
        _sys.stderr.write("========================================\n")

    log.info(
        "pipeline_summary",
        vcf_variants=len(variants),
        annotated=annotated_count,
        missing_annotation=len(missing_keys),
        output_rows=len(results),
    )
    if missing_keys:
        log.warning(
            "annotation_missing_samples",
            count=len(missing_keys),
            first_10=missing_keys[:10],
        )

    write_tsv(results, output_path)

    if skipped:
        if output_path is not None:
            skipped_path = output_path.with_name(
                output_path.stem + ".not_annotated" + output_path.suffix
            )
            write_skipped(skipped, skipped_path)
            log.info("skipped_file_written", path=str(skipped_path), count=len(skipped))
        else:
            log.warning(
                "skipped_records_not_written",
                count=len(skipped),
                reason="main output is stdout; pass -o to also write the not-annotated file",
            )

    return results


def run_single(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    cfg: Config,
    supplement: Optional[list[SupplementEntry]] = None,
) -> ClassificationResult:
    """Classify a single variant and print a rich report to the terminal.

    ``supplement`` is an optional list of manual curator evidence entries for
    this variant (from ``explain --evidence`` / ``--supplement``); they are
    combined with the automated calls per ``cfg.supplement_mode``."""
    from acmg_classifier.annotation.orchestrator import AnnotationOrchestrator
    from acmg_classifier.criteria.registry import CriteriaRegistry
    from acmg_classifier.classification.classifier_2015 import Classifier2015
    from acmg_classifier.classification.classifier_bayesian import ClassifierBayesian
    from acmg_classifier.io.report_writer import print_report

    variant = VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=cfg.assembly)
    orchestrator = AnnotationOrchestrator(cfg)
    registry = CriteriaRegistry(cfg)

    ann = orchestrator.annotate_batch([variant])[variant.key]
    criteria_results = registry.evaluate_all(variant, ann, supplement)

    clf_2015 = Classifier2015()
    clf_bay = ClassifierBayesian()
    classification_2015, rules = clf_2015.classify(criteria_results)
    score, classification_bay = clf_bay.classify(criteria_results)

    result = ClassificationResult(
        variant_id=variant.key,
        criteria_results=criteria_results,
        classification_2015=classification_2015,
        classification_2015_rules=rules,
        bayesian_score=score,
        classification_bayesian=classification_bay,
    )
    print_report(result, variant, ann)
    return result
