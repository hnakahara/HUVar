"""Somatic oncogenicity pipeline: VCF -> annotation -> Horak 2022 criteria -> TSV."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.config import Config
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.somatic.models import (
    ManualOncEvidence, OncCriterion, OncogenicityResult, classify_points,
)

log = structlog.get_logger()

_NOT_MET = {"not_met", "notmet", "no", "false", "0"}


def parse_inline_evidence(items: list[str], variant_id: str) -> list[ManualOncEvidence]:
    """``CRITERION[:not_met][:NOTE]`` items (e.g. ``OP2:retinoblastoma``,
    ``OS2:PMID 12345``, ``OM1:not_met:outside the kinase domain``)."""
    out = []
    for item in items:
        code, _, rest = item.partition(":")
        try:
            crit = OncCriterion(code.strip().upper())
        except ValueError as exc:
            raise ValueError(f"unknown somatic criterion '{code}' "
                             f"(choices: {', '.join(c.value for c in OncCriterion)})") from exc
        flag, _, note = rest.partition(":")
        met = True
        if flag.strip().lower() in _NOT_MET:
            met = False
        elif flag:
            note = rest
        out.append(ManualOncEvidence(variant_id=variant_id, criterion=crit, met=met,
                                     evidence=note.strip()))
    return out


def read_somatic_supplement(path: Path) -> dict[str, list[ManualOncEvidence]]:
    """TSV with columns variant_id, criterion, evidence and optional met (yes/no)."""
    out: dict[str, list[ManualOncEvidence]] = {}
    with path.open(encoding="utf-8", newline="") as fh:
        rd = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        need = {"variant_id", "criterion"}
        if rd.fieldnames is None or not need.issubset(rd.fieldnames):
            raise ValueError(f"somatic supplement {path}: needs columns {sorted(need)}")
        for i, r in enumerate(rd, start=2):
            try:
                crit = OncCriterion(r["criterion"].strip().upper())
            except ValueError as exc:
                raise ValueError(f"{path} line {i}: unknown criterion {r['criterion']!r}") from exc
            met = (r.get("met") or "yes").strip().lower() not in _NOT_MET
            e = ManualOncEvidence(variant_id=r["variant_id"].strip(), criterion=crit, met=met,
                                  evidence=(r.get("evidence") or "").strip())
            out.setdefault(e.variant_id, []).append(e)
    return out


class SomaticClassifier:
    """Germline HUVar criteria (reused for in silico / PM1 / functional data) plus
    the somatic oncogenicity criteria for one variant at a time."""

    def __init__(self, cfg: Config) -> None:
        from acmg_classifier.criteria.registry import CriteriaRegistry
        from acmg_classifier.somatic.evaluator import SomaticEvaluator
        self._cfg = cfg
        self._registry = CriteriaRegistry(cfg)
        self._eval = SomaticEvaluator(cfg)

    def classify(self, variant: VariantRecord, ann, manual=None,
                 warnings: Optional[list[str]] = None) -> OncogenicityResult:
        germline = self._registry.evaluate_all(variant, ann, None)
        crit, ctx = self._eval.evaluate(variant, ann, germline, manual)
        score = sum(r.points for r in crit)
        pc = ann.primary_consequence
        role = ctx["role"]
        warn = list(warnings or [])
        if pc is not None and role is None:
            warn.append(f"GENE_ROLE: {pc.gene_symbol} is not a canonical cancer gene in NCG; "
                        "OVS1/OM2/OS2 are not applied")
        elif role is not None and role.role in ("dual", "unresolved"):
            warn.append(f"GENE_ROLE: {role.describe()} - OVS1/OS2 not applied; review the role in "
                        "this tumour type")
        if ann.gnomad is not None and (ann.gnomad.ac or 0) > 0:
            warn.append("POPULATION: present in gnomAD; for haematological genes consider clonal "
                        "haematopoiesis in the controls before relying on frequency evidence")
        return OncogenicityResult(
            variant_id=variant.key, chrom=variant.chrom, pos=variant.pos, ref=variant.ref,
            alt=variant.alt, filter=variant.filter,
            transcript_id=pc.transcript_id if pc else None,
            gene_symbol=pc.gene_symbol if pc else None,
            hgvs_c=pc.hgvs_c if pc else None, hgvs_p=pc.hgvs_p if pc else None,
            consequence=pc.consequence.value if pc else None,
            gene_role=role.role if role else "not_cancer_gene",
            gene_role_source=role.source if role else "",
            gene_role_note=role.note if role else "",
            hotspot=ctx["hotspot"],
            criteria=crit, score=score, classification=classify_points(score),
            annotation=ann, warnings=warn,
        )


def run_somatic_pipeline(vcf_path: Path, cfg: Config, output_path: Optional[Path] = None,
                         supplement_path: Optional[Path] = None) -> list[OncogenicityResult]:
    from acmg_classifier.annotation.orchestrator import AnnotationOrchestrator
    from acmg_classifier.io.vcf_reader import detect_assembly_from_header, read_vcf
    from acmg_classifier.models.annotation import AnnotationData
    from acmg_classifier.models.variant import large_indel_warning, unsupported_allele_reason
    from acmg_classifier.somatic.report import write_somatic_tsv
    from acmg_classifier.utils.progress import track

    if cfg.assembly is None:
        detected = detect_assembly_from_header(vcf_path)
        if detected:
            cfg = cfg.model_copy(update={"assembly": detected})
    manual = read_somatic_supplement(supplement_path) if supplement_path else {}
    records = list(read_vcf(vcf_path, cfg.assembly))
    variants = [v for v in records if unsupported_allele_reason(v.ref, v.alt) is None]
    skipped = len(records) - len(variants)
    if skipped:
        log.warning("somatic_skipped_records", count=skipped,
                    reason="no ALT / symbolic / breakend alleles are not classified")
    annotations = AnnotationOrchestrator(cfg).annotate_batch(variants)
    clf = SomaticClassifier(cfg)
    results = []
    for v in track(variants, "Classifying (somatic)", total=len(variants)):
        ann = annotations.get(v.key)
        warn = []
        if ann is None:
            ann, warn = AnnotationData(), ["VEP returned no annotation for this variant"]
        big = large_indel_warning(v.ref, v.alt)
        if big:
            warn.append(big)
        results.append(clf.classify(v, ann, manual.get(v.key), warn))
    write_somatic_tsv(results, output_path)
    log.info("somatic_summary", variants=len(results))
    return results


def run_somatic_single(chrom: str, pos: int, ref: str, alt: str, cfg: Config,
                       manual=None, transcript: Optional[str] = None,
                       print_result: bool = True) -> OncogenicityResult:
    from acmg_classifier.pipeline.pipeline import annotate_one
    from acmg_classifier.somatic.report import print_somatic_report

    variant, ann = annotate_one(chrom, pos, ref, alt, cfg)
    if transcript:
        ann = ann.for_transcript(transcript)
    result = SomaticClassifier(cfg).classify(variant, ann, manual)
    if print_result:
        print_somatic_report(result)
    return result
