"""Main pipeline: reads VCF → annotates → classifies → writes output."""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import structlog

from acmg_classifier.config import Config
from acmg_classifier.io.vcf_reader import read_vcf, detect_assembly_from_header
from acmg_classifier.io.supplement_reader import read_supplement
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.classification import ClassificationResult

log = structlog.get_logger()


def run_pipeline(
    vcf_path: Path,
    cfg: Config,
    output_path: Optional[Path] = None,
    supplement_path: Optional[Path] = None,
) -> list[ClassificationResult]:
    """Classify all variants in *vcf_path* and write TSV to *output_path*."""
    from acmg_classifier.annotation.orchestrator import AnnotationOrchestrator
    from acmg_classifier.criteria.registry import CriteriaRegistry
    from acmg_classifier.classification.classifier_2015 import Classifier2015
    from acmg_classifier.classification.classifier_bayesian import ClassifierBayesian
    from acmg_classifier.io.tsv_writer import write_tsv, write_skipped

    # Detect assembly if not explicit
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
    # ALT='.' records carry no alternate allele -> they are not variants and cannot
    # be classified. Exclude them from the main output and report them separately.
    variants = [v for v in all_records if v.alt and v.alt != "."]
    skipped = [v for v in all_records if not (v.alt and v.alt != ".")]
    log.info(
        "variants_loaded",
        total=len(all_records), classifiable=len(variants), skipped_no_alt=len(skipped),
    )

    annotations = orchestrator.annotate_batch(variants)
    log.info("annotations_returned", count=len(annotations))

    missing_keys: list[str] = []
    annotated_count = 0

    for variant in variants:
        ann = annotations.get(variant.key)
        if ann is None:
            missing_keys.append(variant.key)
            ann = AnnotationData()  # empty annotation -> still output the row
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
) -> ClassificationResult:
    """Classify a single variant and print a rich report to the terminal."""
    from acmg_classifier.annotation.orchestrator import AnnotationOrchestrator
    from acmg_classifier.criteria.registry import CriteriaRegistry
    from acmg_classifier.classification.classifier_2015 import Classifier2015
    from acmg_classifier.classification.classifier_bayesian import ClassifierBayesian
    from acmg_classifier.io.report_writer import print_report

    variant = VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=cfg.assembly)
    orchestrator = AnnotationOrchestrator(cfg)
    registry = CriteriaRegistry(cfg)

    ann = orchestrator.annotate_batch([variant])[variant.key]
    criteria_results = registry.evaluate_all(variant, ann)

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
