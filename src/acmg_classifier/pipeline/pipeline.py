"""Main pipeline: reads VCF → annotates → classifies → writes output."""
from __future__ import annotations
from pathlib import Path
from typing import TYPE_CHECKING, Optional
import structlog

from acmg_classifier.config import Config
from acmg_classifier.io.vcf_reader import read_vcf, detect_assembly_from_header
from acmg_classifier.io.supplement_reader import read_supplement
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.classification import ClassificationResult
from acmg_classifier.utils.progress import track

if TYPE_CHECKING:
    from acmg_classifier.models.annotation import AnnotationData

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
        # Curated evidence is in use: also import the criteria a ClinGen expert
        # panel applied to the variant in ClinVar (>=3 stars).
        cfg = cfg.model_copy(update={"curated_evidence": True})

    orchestrator = AnnotationOrchestrator(cfg)
    registry = CriteriaRegistry(cfg)
    clf_2015 = Classifier2015()
    clf_bay = ClassifierBayesian()

    # Alleles are normalised (trimmed + left-aligned against the genome FASTA) so
    # that indels written with padding bases still match ClinVar / gnomAD / SpliceAI.
    all_records = list(read_vcf(vcf_path, cfg.assembly, fasta=cfg.genome_fasta))
    # ALT='.' records carry no alternate allele — they are not variants and
    # cannot be classified. Exclude them from the main output and write them
    # to a `.not_annotated` sidecar so users can see *which* rows were
    # skipped (silent dropping would surprise reviewers).
    # Symbolic (<DEL>, <DUP>, <CNV>...), breakend and '*' alleles are structural
    # variants or placeholders that HUVar cannot classify; they go to the same
    # sidecar with the reason (Reviewer 1, Genome Medicine revision).
    from acmg_classifier.models.variant import unsupported_allele_reason
    skip_reason = {v.key: unsupported_allele_reason(v.ref, v.alt) for v in all_records}
    variants = [v for v in all_records if skip_reason[v.key] is None]
    skipped = [v for v in all_records if skip_reason[v.key] is not None]
    log.info(
        "variants_loaded",
        total=len(all_records), classifiable=len(variants), skipped=len(skipped),
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

    missing_keys = [v.key for v in variants if v.key not in annotations]
    annotated_count = len(variants) - len(missing_keys)

    if profiler is not None:
        # Per-evaluator wall time next to the cProfile output; the profile run
        # classifies in this process so cProfile sees the evaluators.
        registry.enable_timing()
    ctx = _ClassifyContext(cfg, registry, clf_2015, clf_bay, annotations, supplement)
    results = _classify_all(variants, ctx, serial=profiler is not None)

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
        _sys.stderr.write(_format_timings(registry.timings(), len(variants)))

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

    # Data provenance (HUVar version + resource versions) next to the output,
    # and a warning when the local ClinVar release is stale.
    try:
        from acmg_classifier.provenance import collect, write_sidecar
        if output_path is not None:
            side = write_sidecar(cfg, output_path)
            log.info("provenance_written", path=str(side))
        for w in collect(cfg)["warnings"]:
            log.warning("data_version_warning", message=w)
    except Exception as exc:  # noqa: BLE001 — provenance must never break a run
        log.warning("provenance_failed", error=str(exc))

    if skipped:
        if output_path is not None:
            skipped_path = output_path.with_name(
                output_path.stem + ".not_annotated" + output_path.suffix
            )
            write_skipped(skipped, skipped_path, reasons=skip_reason)
            log.info("skipped_file_written", path=str(skipped_path), count=len(skipped))
        else:
            log.warning(
                "skipped_records_not_written",
                count=len(skipped),
                reason="main output is stdout; pass -o to also write the not-annotated file",
            )

    return results


class _ClassifyContext:
    """Everything the per-variant classification needs (shared with forked
    worker processes copy-on-write)."""

    def __init__(self, cfg, registry, clf_2015, clf_bay, annotations, supplement):
        self.cfg = cfg
        self.registry = registry
        self.clf_2015 = clf_2015
        self.clf_bay = clf_bay
        self.annotations = annotations
        self.supplement = supplement


def _classify_record(variant: VariantRecord, ctx: _ClassifyContext) -> ClassificationResult:
    """Evaluate the criteria for one annotated VCF record and build its result."""
    from acmg_classifier.models.annotation import AnnotationData
    from acmg_classifier.classification.gene_validity import apply_gene_validity
    from acmg_classifier.models.variant import large_indel_warning
    from acmg_classifier.local_db.erepo_db import vcep_record, vcep_warning

    cfg, registry, clf_2015, clf_bay = ctx.cfg, ctx.registry, ctx.clf_2015, ctx.clf_bay
    ann = ctx.annotations.get(variant.key)
    # Annotation failures are non-fatal — we emit a row with an empty
    # AnnotationData and surface the warning in the output so the user
    # sees there *was* a variant rather than a silent gap. This matches
    # clinical-lab expectations of "explain every input row".
    if ann is None:
        ann = AnnotationData()
        extra_warnings = ["VEP returned no annotation for this variant"]
    else:
        extra_warnings = []

    # Curator evidence may be keyed by the normalised or by the VCF's own representation
    sup_entries = (ctx.supplement.get(variant.key)
                   or (ctx.supplement.get(variant.input_variant) if variant.input_variant else None)
                   or [])
    criteria_results = registry.evaluate_all(variant, ann, sup_entries)
    classification_2015, rules = clf_2015.classify(criteria_results)
    score, classification_bay = clf_bay.classify(criteria_results)

    def _evaluate_on(ann_tx, _variant=variant, _sup=sup_entries):
        crit = registry.evaluate_all(_variant, ann_tx, _sup)
        c2015, r2015 = clf_2015.classify(crit)
        sc, cbay = clf_bay.classify(crit)
        tpc = ann_tx.primary_consequence
        c2015, cbay, _w = apply_gene_validity(
            tpc.gene_symbol if tpc else None, c2015, cbay, cfg)
        return crit, c2015, r2015, sc, cbay

    pc = ann.primary_consequence
    classification_2015, classification_bay, gdv_warn = apply_gene_validity(
        pc.gene_symbol if pc else None, classification_2015, classification_bay, cfg,
    )
    extra_warnings = list(extra_warnings) + gdv_warn
    big = large_indel_warning(variant.ref, variant.alt)
    if big:
        extra_warnings.append(big)
    extra_warnings += vcep_warning(cfg.assembly.value, variant.chrom, variant.pos,
                                   variant.ref, variant.alt)
    return ClassificationResult(
        variant_id=variant.key,
        input_variant=variant.input_variant,
        vcep=vcep_record(cfg.assembly.value, variant.chrom, variant.pos,
                         variant.ref, variant.alt),
        chrom=variant.chrom,
        pos=variant.pos,
        ref=variant.ref,
        alt=variant.alt,
        filter=variant.filter,
        transcript_id=pc.transcript_id if pc else None,
        gene_symbol=pc.gene_symbol if pc else None,
        hgvs_c=pc.hgvs_c if pc else None,
        hgvs_p=pc.hgvs_p if pc else None,
        mane_plus_clinical=plus_clinical_evaluations(ann, cfg, criteria_results, _evaluate_on),
        annotation=ann,
        criteria_results=criteria_results,
        classification_2015=classification_2015,
        classification_2015_rules=rules,
        bayesian_score=score,
        classification_bayesian=classification_bay,
        warnings=extra_warnings,
    )


# Below this many variants the fork / pickling overhead outweighs the gain.
_MIN_PARALLEL_VARIANTS = 200
# Variants classified in the parent before forking: the first record of each
# gene, so lazily loaded per-gene tables (and the large shared ones) are in
# memory once and inherited by every worker instead of loaded per worker.
_MAX_WARMUP = 256

# Set in the parent just before the pool forks; workers read it (fork start
# method only, so it is inherited, never pickled).
_FORK_CTX: dict = {}


def classify_process_count(cfg: Config, n_variants: int) -> int:
    """Worker processes for classifying ``n_variants`` (1 = in-process).

    Parallel classification needs the ``fork`` start method (Linux), a
    non-daemonic process and a batch large enough to pay for the pool. Inside a
    Celery prefork worker (a daemonic billiard process, e.g. the HUVar_app batch
    jobs) the automatic setting stays in-process so a web job does not
    multiply its CPU use unasked; an explicit ``classify_workers`` still applies
    there."""
    import multiprocessing as mp
    import sys

    explicit = int(getattr(cfg, "classify_workers", 0) or 0)
    want = explicit or int(cfg.workers or 1)
    if want <= 1 or n_variants < _MIN_PARALLEL_VARIANTS:
        return 1
    if not sys.platform.startswith("linux") or "fork" not in mp.get_all_start_methods():
        return 1
    if mp.current_process().daemon:  # stdlib refuses children of daemons
        return 1
    if not explicit and "billiard" in sys.modules:
        try:
            from billiard.process import current_process as _bp  # type: ignore
            if _bp().daemon:
                return 1
        except Exception:  # noqa: BLE001
            pass
    return min(want, max(1, n_variants // 50))


def _worker_init() -> None:
    from acmg_classifier.utils import progress
    # Only the parent draws the "Classifying" bar.
    progress.set_enabled(False)


def _classify_range(bounds: tuple[int, int]) -> list[tuple[int, ClassificationResult]]:
    ctx: _ClassifyContext = _FORK_CTX["ctx"]
    variants = _FORK_CTX["variants"]
    done = _FORK_CTX["done"]
    out = []
    for i in range(*bounds):
        if i in done:
            continue
        res = _classify_record(variants[i], ctx)
        # The parent holds the identical AnnotationData (evaluators never
        # mutate it) and re-attaches it, so it is not pickled back.
        res.annotation = None
        out.append((i, res))
    return out


def _warmup_indices(variants: list[VariantRecord], annotations: dict) -> list[int]:
    seen: set = set()
    picked: list[int] = []
    for i, v in enumerate(variants):
        ann = annotations.get(v.key)
        pc = ann.primary_consequence if ann is not None else None
        gene = pc.gene_symbol if pc is not None else None
        if gene in seen:
            continue
        seen.add(gene)
        picked.append(i)
        if len(picked) >= _MAX_WARMUP:
            break
    return picked


def _classify_all(variants: list[VariantRecord], ctx: _ClassifyContext,
                  serial: bool = False) -> list[ClassificationResult]:
    """Classify ``variants`` in order, in worker processes when worthwhile.

    The criteria evaluation is pure Python (GIL-bound), so threads do not help;
    forked processes share the loaded reference tables copy-on-write. Results
    are identical to the in-process loop and returned in input order."""
    from acmg_classifier.models.annotation import AnnotationData
    from acmg_classifier.utils.progress import progress_bar

    n = len(variants)
    n_proc = 1 if serial else classify_process_count(ctx.cfg, n)
    if n_proc <= 1:
        # Classification is CPU-cheap relative to annotation but still benefits
        # from a progress bar when there are many variants — without it the user
        # only sees the annotation bar finish and then a silent gap until the
        # final summary line.
        return [_classify_record(v, ctx) for v in track(variants, "Classifying", total=n)]

    import gc
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from concurrent.futures.process import BrokenProcessPool

    slots: list[Optional[ClassificationResult]] = [None] * n
    warm = _warmup_indices(variants, ctx.annotations)
    for i in track(warm, "Classifying (one variant per gene, pre-loading tables)",
                   total=len(warm)):
        slots[i] = _classify_record(variants[i], ctx)
    # Small chunks keep the bar moving and balance genes of unequal cost.
    chunk = max(16, min(512, n // (n_proc * 16) or 16))
    ranges = [(s, min(s + chunk, n)) for s in range(0, n, chunk)]
    log.info("classify_parallel", processes=n_proc, chunks=len(ranges),
             chunk_size=chunk, warmup=len(warm))
    _FORK_CTX.update(ctx=ctx, variants=variants, done=frozenset(warm))
    # Freeze the objects loaded so far so the workers' garbage collector does
    # not touch (and thereby copy) the shared pages.
    gc.collect()
    gc.freeze()
    try:
        with ProcessPoolExecutor(max_workers=n_proc, mp_context=mp.get_context("fork"),
                                 initializer=_worker_init) as pool:
            # Every task is submitted (and so every worker forked) BEFORE the
            # progress bar starts: the bar's refresh thread must not exist at
            # fork time, or a lock it holds would stay locked in the children.
            futures = [pool.submit(_classify_range, r) for r in ranges]
            with progress_bar("Classifying", total=n) as advance:
                advance(len(warm))
                try:
                    for fut in as_completed(futures):
                        part = fut.result()
                        for i, res in part:
                            ann = ctx.annotations.get(variants[i].key)
                            res.annotation = ann if ann is not None else AnnotationData()
                            slots[i] = res
                        advance(len(part))
                except BrokenProcessPool as exc:
                    # A worker died (e.g. killed for memory). Finish in-process
                    # below rather than losing the run.
                    log.warning("classify_parallel_failed", error=str(exc),
                                hint="set ACMG_CLASSIFY_WORKERS lower (or 1) if workers "
                                     "run out of memory")
    finally:
        _FORK_CTX.clear()
        gc.unfreeze()
    missing = [i for i, r in enumerate(slots) if r is None]
    if missing:
        for i in track(missing, "Classifying (in-process)", total=len(missing)):
            slots[i] = _classify_record(variants[i], ctx)
    return slots  # type: ignore[return-value]


def _format_timings(timings: dict[str, float], n_variants: int) -> str:
    """Per-evaluator wall-time table for ``classify --profile``."""
    if not timings:
        return ""
    total = sum(timings.values()) or 1.0
    lines = ["\n===== criteria evaluation time per evaluator =====",
             f"{'evaluator':<34}{'total s':>10}{'ms/variant':>12}{'share':>8}"]
    for name, sec in timings.items():
        per = sec * 1000.0 / max(1, n_variants)
        lines.append(f"{name:<34}{sec:>10.2f}{per:>12.3f}{sec / total:>8.1%}")
    lines.append("=" * 64 + "\n")
    return "\n".join(lines)


def annotate_one(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    cfg: Config,
) -> tuple["VariantRecord", "AnnotationData"]:
    """Annotate a single variant — the heavy, CSpec-independent step.

    The returned annotation is reusable across several CSpec evaluations via
    :func:`classify_annotated`, so the expensive VEP / splice / in-silico work
    runs only once even when a variant is classified under multiple CSpecs."""
    from acmg_classifier.annotation.orchestrator import AnnotationOrchestrator
    from acmg_classifier.models.variant import unsupported_allele_reason

    why = unsupported_allele_reason(ref, alt)
    if why is not None:
        raise ValueError(f"{chrom}:{pos}:{ref}:{alt}: {why}")
    variant = VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=cfg.assembly)
    # Same normalisation as the VCF path, so a padded indel typed on the command line
    # (e.g. GAAA GAA) is looked up as GA>G.
    from acmg_classifier.annotation.normalizer import normalize_record, open_fasta
    fasta = open_fasta(cfg.genome_fasta)
    try:
        variant = normalize_record(variant, fasta)
    finally:
        if fasta is not None:
            fasta.close()
    orchestrator = AnnotationOrchestrator(cfg)
    ann = orchestrator.annotate_batch([variant])[variant.key]
    return variant, ann


def _large_indel(variant) -> list[str]:
    from acmg_classifier.models.variant import large_indel_warning
    w = large_indel_warning(variant.ref, variant.alt)
    return [w] if w else []


def classify_annotated(
    variant: "VariantRecord",
    ann: "AnnotationData",
    cfg: Config,
    supplement: Optional[list[SupplementEntry]] = None,
) -> ClassificationResult:
    """Evaluate criteria and classify a pre-annotated variant under ``cfg``.

    The CSpec-dependent step: thresholds come from ``cfg.disease_prevalence_tsv``,
    so the *same* annotation can be classified under different CSpecs by passing
    a ``cfg`` whose ``disease_prevalence_tsv_override`` points at the desired
    per-CSpec overlay table."""
    from acmg_classifier.criteria.registry import CriteriaRegistry
    from acmg_classifier.classification.classifier_2015 import Classifier2015
    from acmg_classifier.classification.classifier_bayesian import ClassifierBayesian

    from acmg_classifier.classification.gene_validity import apply_gene_validity

    registry = CriteriaRegistry(cfg)
    clf_2015, clf_bay = Classifier2015(), ClassifierBayesian()

    def _evaluate_on(ann_tx):
        crit = registry.evaluate_all(variant, ann_tx, supplement)
        c2015, r2015 = clf_2015.classify(crit)
        sc, cbay = clf_bay.classify(crit)
        tpc = ann_tx.primary_consequence if ann_tx is not None else None
        c2015, cbay, warn = apply_gene_validity(
            tpc.gene_symbol if tpc else None, c2015, cbay, cfg)
        return crit, c2015, r2015, sc, cbay, warn

    criteria_results, classification_2015, rules, score, classification_bay, gdv_warn = (
        _evaluate_on(ann))
    pc = ann.primary_consequence if ann is not None else None

    # Plus Clinical evaluations are reported next to the default (MANE Select)
    # evaluation only; a result that is itself evaluated on a Plus Clinical
    # transcript just lists the Plus Clinical notation.
    if ann is not None and not ann.evaluation_transcript:
        pcs = plus_clinical_evaluations(
            ann, cfg, criteria_results, lambda a: _evaluate_on(a)[:5])
    else:
        pcs = ann.mane_plus_clinical_summary() if ann is not None else []

    from acmg_classifier.local_db.erepo_db import vcep_record
    return ClassificationResult(
        variant_id=variant.key,
        vcep=vcep_record(cfg.assembly.value, variant.chrom, variant.pos,
                         variant.ref, variant.alt),
        mane_plus_clinical=pcs,
        criteria_results=criteria_results,
        classification_2015=classification_2015,
        classification_2015_rules=rules,
        bayesian_score=score,
        classification_bayesian=classification_bay,
        warnings=gdv_warn + _large_indel(variant) + _vcep_warning(cfg, variant)
        + _cspec_note(cfg, pc) + transcript_note(ann),
    )


def _criterion_state(r) -> str:
    return r.strength.value if (r.triggered and not r.suppressed) else "NotMet"


def criteria_changes(base, other) -> str:
    """Criteria whose outcome differs between two evaluations, e.g.
    ``"PP3:NotMet->Moderate; PM1:Moderate->NotMet"`` (order of ``base``)."""
    b = {r.criterion: _criterion_state(r) for r in base}
    o = {r.criterion: _criterion_state(r) for r in other}
    order = [r.criterion for r in base] + [c for c in o if c not in b]
    seen, out = set(), []
    for c in order:
        if c in seen:
            continue
        seen.add(c)
        x, y = b.get(c, "NotMet"), o.get(c, "NotMet")
        if x != y:
            out.append(f"{c.value}:{x}->{y}")
    return "; ".join(out)


def plus_clinical_evaluations(ann, cfg: Config, base_criteria, evaluate) -> list:
    """MANE Plus Clinical summary of ``ann``, each entry classified with the
    criteria evaluated on that transcript (``evaluate(ann_tx)`` returns
    ``(criteria, classification_2015, rules, score, classification_bayesian)``).

    Skipped (notation only) when ``cfg.evaluate_plus_clinical`` is False."""
    summary = ann.mane_plus_clinical_summary() if ann is not None else []
    if not summary or getattr(cfg, "evaluate_plus_clinical", True) is not True:
        return summary
    for info in summary:
        try:
            crit, c2015, rules, score, cbay = evaluate(ann.for_transcript(info.transcript_id))
        except Exception as exc:  # never fail the MANE Select result
            log.warning("plus_clinical_eval_failed", transcript=info.transcript_id,
                        error=str(exc))
            continue
        info.classification_2015 = c2015.value
        info.rules_2015 = rules
        info.bayesian_score = score
        info.classification_bayesian = cbay.value
        info.criteria_changes = criteria_changes(base_criteria, crit)
    return summary


def transcript_note(ann) -> list[str]:
    """Warning for a result evaluated on a non-default transcript: which
    numbering the VCEP position rules / ClinVar codon comparisons used."""
    if ann is None or not getattr(ann, "evaluation_transcript", None):
        return []
    pc, rc = ann.primary_consequence, ann.rule_consequence
    if pc is None:
        return []
    if pc.transcript_id.split(".")[0] != ann.evaluation_transcript.split(".")[0]:
        return [f"TRANSCRIPT: {ann.evaluation_transcript} has no VEP consequence for this "
                f"variant; evaluated on {pc.transcript_id}"]
    if rc is None or rc is pc:
        return [f"TRANSCRIPT: evaluated on {pc.transcript_id}"]
    note = (f"TRANSCRIPT: evaluated on {pc.transcript_id} (MANE Plus Clinical); VCEP "
            "position rules and ClinVar codon comparisons (PM1, PS1, PM5, BA1/BS1/BP4 "
            f"variant lists, VCEP PVS1 bands, PP3/BP4 domains) use MANE Select "
            f"{rc.transcript_id} ({(rc.hgvs_c or '').split(':')[-1] or rc.consequence.value})")
    if pc.protein_position is not None and rc.protein_position is None:
        note += ("; this residue has no MANE Select position, so those rules cannot be "
                 "met — review against the Plus Clinical isoform manually")
    return [note]


def _cspec_note(cfg: Config, pc) -> list[str]:
    """Name the CSpec chosen with --cspec for this variant's gene."""
    sel = getattr(cfg, "cspec_selection", None)
    gene = getattr(pc, "gene_symbol", None) if pc is not None else None
    if isinstance(sel, dict) and gene in sel:
        return [f"CSPEC: {gene} evaluated under the '{sel[gene]}' specification (--cspec)"]
    return []


def _vcep_warning(cfg: Config, variant) -> list[str]:
    from acmg_classifier.local_db.erepo_db import vcep_warning
    try:
        return vcep_warning(cfg.assembly.value, variant.chrom, variant.pos,
                            variant.ref, variant.alt)
    except Exception:  # never fail a classification on the lookup
        return []


def run_single(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    cfg: Config,
    supplement: Optional[list[SupplementEntry]] = None,
    transcript: Optional[str] = None,
) -> ClassificationResult:
    """Classify a single variant and print a rich report to the terminal.

    ``transcript`` evaluates the criteria on that transcript (e.g. the gene's
    MANE Plus Clinical) instead of the default MANE Select.

    ``supplement`` is an optional list of manual curator evidence entries for
    this variant (from ``explain --evidence`` / ``--supplement``); they are
    combined with the automated calls per ``cfg.supplement_mode``."""
    from acmg_classifier.io.report_writer import print_report

    if supplement:
        # Manual evidence supplied: expert-panel criteria in ClinVar are used too.
        cfg = cfg.model_copy(update={"curated_evidence": True})
    variant, ann = annotate_one(chrom, pos, ref, alt, cfg)
    if transcript:
        ann = ann.for_transcript(transcript)
    result = classify_annotated(variant, ann, cfg, supplement)
    print_report(result, variant, ann)
    return result
