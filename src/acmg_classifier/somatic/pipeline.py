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

    def __init__(self, cfg: Config, tumor_type: Optional[str] = None) -> None:
        from acmg_classifier.criteria.registry import CriteriaRegistry
        from acmg_classifier.somatic.civic_clinical import load_clinical
        from acmg_classifier.somatic.evaluator import SomaticEvaluator
        from acmg_classifier.somatic.tumor import disease_ontology
        from acmg_classifier.somatic.tumor import tumor_type as _tumor
        self._cfg = cfg
        self._registry = CriteriaRegistry(cfg)
        self._eval = SomaticEvaluator(cfg)
        self._tumor = _tumor(tumor_type, cfg)          # ValueError for an unknown code
        self._do = disease_ontology(cfg)
        path = getattr(cfg, "civic_clinical_tsv", None)
        self._civic_path = path
        self._civic_clinical = load_clinical(path) if path is not None else {}

    def _amp(self, variant: VariantRecord, ann, crit, classification):
        from acmg_classifier.local_db.clinvar_sqlite import query_somatic_impact
        from acmg_classifier.somatic.amp import assign, civic_calls, clinvar_calls
        from acmg_classifier.somatic.civic_clinical import VariantFacts, match_records
        from acmg_classifier.somatic.evaluator import _NULL
        from acmg_classifier.somatic.models import OncCriterion, Oncogenicity
        from acmg_classifier.somatic.resources import short_protein_change

        pc = ann.primary_consequence
        notes = []
        calls = []
        if pc is not None:
            ref, _, alt = (pc.amino_acids or "/").partition("/")
            if pc.consequence.value == "missense_variant" and ref and alt and pc.protein_position:
                change = f"{ref}{pc.protein_position}{alt}"
            else:
                change = short_protein_change(pc.hgvs_p)

            def _num(x):
                try:
                    return int(str(x).split("/")[0])
                except (TypeError, ValueError):
                    return None
            facts = VariantFacts(
                gene=pc.gene_symbol, protein_change=change, ref_aa=ref or None,
                position=pc.protein_position, consequence=pc.consequence.value,
                exon=_num(pc.exon), intron=_num(pc.intron), is_null=pc.consequence in _NULL,
                oncogenic=classification in (Oncogenicity.ONCOGENIC, Oncogenicity.LIKELY_ONCOGENIC))
            if self._civic_clinical:
                calls += civic_calls(match_records(self._civic_clinical, facts), self._tumor, self._do)
            else:
                notes.append("CIViC clinical table not installed (setup_data.py --only civic)")
        sci = query_somatic_impact(self._cfg.clinvar_sqlite, variant.chrom, variant.pos,
                                   variant.ref, variant.alt)
        if sci is None:
            notes.append("ClinVar DB without somatic clinical impact (rebuild with "
                         "--only clinvar-sqlite --force-clinvar)")
        else:
            calls += clinvar_calls(sci, self._tumor, self._do)
        by = {r.criterion: r for r in crit}
        common = any(by[c].active for c in (OncCriterion.SBVS1, OncCriterion.SBS1))
        benign = classification in (Oncogenicity.BENIGN, Oncogenicity.LIKELY_BENIGN)
        res = assign(calls, self._tumor, self._do, benign, common)
        if notes:
            res.note = "; ".join([res.note] + notes if res.note else notes)
        return res

    def classify(self, variant: VariantRecord, ann, manual=None,
                 warnings: Optional[list[str]] = None,
                 germline: Optional[list] = None) -> OncogenicityResult:
        """``germline``: the HUVar criteria results of this variant and annotation
        when already evaluated (without curator evidence) -- reused instead of
        evaluating them again."""
        if germline is None:
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
        classification = classify_points(score)
        amp = self._amp(variant, ann, crit, classification)
        return OncogenicityResult(
            tumor_type=amp.tumor, amp_tier=amp.tier, amp_level=amp.level,
            amp_evidence=amp.evidence_text(), amp_note=amp.note, amp_candidates=amp.candidates,
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
            criteria=crit, score=score, classification=classification,
            annotation=ann, warnings=warn,
        )


def _germline_config(cfg: Config, curated: bool) -> Config:
    # Curator evidence in use: also import the criteria a ClinGen expert panel
    # applied in ClinVar (>=3 stars), as `classify --supplement` does.
    return cfg.model_copy(update={"curated_evidence": True}) if curated else cfg


def run_somatic_pipeline(vcf_path: Path, cfg: Config, output_path: Optional[Path] = None,
                         supplement_path: Optional[Path] = None,
                         tumor_type: Optional[str] = None,
                         germline: bool = True,
                         somatic_supplement_path: Optional[Path] = None,
                         limit: Optional[int] = None) -> list[OncogenicityResult]:
    """Oncogenicity (and, with ``germline``, the germline ACMG/AMP classification of
    `classify`) of every classifiable VCF record, from one annotation pass.

    ``supplement_path`` is the ACMG curator TSV of `classify --supplement`, which
    applies to the germline classification only; ``somatic_supplement_path`` the
    somatic manual evidence (OVS1..SBP2). ``limit`` keeps the first N classifiable
    records."""
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
    manual = read_somatic_supplement(somatic_supplement_path) if somatic_supplement_path else {}
    records = list(read_vcf(vcf_path, cfg.assembly, fasta=cfg.genome_fasta))
    variants = [v for v in records if unsupported_allele_reason(v.ref, v.alt) is None]
    skipped = len(records) - len(variants)
    if skipped:
        log.warning("somatic_skipped_records", count=skipped,
                    reason="no ALT / symbolic / breakend alleles are not classified")
    if limit is not None and limit < len(variants):
        log.info("variants_limited", from_count=len(variants), to_count=limit)
        variants = variants[:limit]
    clf = SomaticClassifier(cfg, tumor_type)
    annotations = AnnotationOrchestrator(cfg).annotate_batch(variants)
    germ_results = None
    if germline:
        from acmg_classifier.io.supplement_reader import read_supplement
        from acmg_classifier.pipeline.pipeline import classify_variants
        germ_sup = read_supplement(supplement_path) if supplement_path else {}
        germ_results = classify_variants(
            variants, annotations, _germline_config(cfg, bool(supplement_path)), germ_sup)
    # Without curator evidence the germline criteria are exactly what the somatic
    # evaluation reads (in silico, PM1, functional data), so they are reused.
    reuse = germ_results is not None and not supplement_path
    results = []
    for i, v in enumerate(track(variants, "Classifying (somatic)", total=len(variants))):
        ann = annotations.get(v.key)
        warn = []
        if ann is None:
            ann, warn = AnnotationData(), ["VEP returned no annotation for this variant"]
        big = large_indel_warning(v.ref, v.alt)
        if big:
            warn.append(big)
        entries = manual.get(v.key) or (manual.get(v.input_variant) if v.input_variant else None)
        if v.input_variant:
            warn.append(f"NORMALIZED: VCF allele {v.input_variant} evaluated as {v.key}")
        g = germ_results[i] if germ_results is not None else None
        res = clf.classify(v, ann, entries, warn,
                           germline=g.criteria_results if reuse else None)
        res.input_variant = v.input_variant
        res.sample_calls = list(v.sample_calls)
        res.germline = g
        results.append(res)
    write_somatic_tsv(results, output_path, germline=germline)
    log.info("somatic_summary", variants=len(results))
    return results


def run_somatic_single(chrom: str, pos: int, ref: str, alt: str, cfg: Config,
                       manual=None, transcript: Optional[str] = None,
                       print_result: bool = True,
                       tumor_type: Optional[str] = None,
                       germline: bool = False,
                       supplement=None) -> OncogenicityResult:
    """``manual``: somatic evidence (ManualOncEvidence). ``germline`` also classifies
    the variant under ACMG/AMP (as `explain`; with the curator entries
    ``supplement``) and prints that report after the oncogenicity one."""
    from acmg_classifier.pipeline.pipeline import annotate_one, classify_annotated
    from acmg_classifier.somatic.report import print_somatic_report

    clf = SomaticClassifier(cfg, tumor_type)   # validates the tumour type before annotating
    variant, ann = annotate_one(chrom, pos, ref, alt, cfg)
    if transcript:
        ann = ann.for_transcript(transcript)
    g = None
    if germline:
        g = classify_annotated(variant, ann, _germline_config(cfg, bool(supplement)),
                               supplement or None)
    result = clf.classify(variant, ann, manual,
                          germline=g.criteria_results if g is not None and not supplement
                          else None)
    result.germline = g
    if print_result:
        print_somatic_report(result)
        if g is not None:
            from acmg_classifier.io.report_writer import print_report
            print_report(g, variant, ann)
    return result
