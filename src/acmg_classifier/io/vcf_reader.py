from __future__ import annotations
import re
from pathlib import Path
from typing import Iterator, Optional
from cyvcf2 import VCF  # type: ignore
import structlog

from acmg_classifier.models.enums import Assembly
from acmg_classifier.models.variant import VariantRecord

log = structlog.get_logger()


_ASSEMBLY_PATTERNS = {
    Assembly.GRCH38: re.compile(r"GRCh38|hg38|GRCh38\.p", re.IGNORECASE),
    Assembly.GRCH37: re.compile(r"GRCh37|hg19|GRCh37\.p|b37", re.IGNORECASE),
}


def detect_assembly_from_header(vcf_path: Path) -> Optional[Assembly]:
    """Scan VCF header for assembly information."""
    vcf = VCF(str(vcf_path))
    for line in vcf.raw_header.splitlines():
        if line.startswith("##reference") or line.startswith("##contig"):
            for assembly, pattern in _ASSEMBLY_PATTERNS.items():
                if pattern.search(line):
                    vcf.close()
                    return assembly
    vcf.close()
    return None


_MISSING_INT = -2147483648   # cyvcf2's missing-value marker for integer FORMAT fields


def _format(variant, name: str):
    try:
        return variant.format(name)
    except Exception:  # noqa: BLE001 -- field absent from this record
        return None


def _int_or_none(value) -> Optional[int]:
    try:
        v = int(value)
    except (TypeError, ValueError):
        return None
    return None if v < 0 else v


def sample_calls_for(variant, samples: list[str], alt_index: int) -> list:
    """SampleCall for every sample whose GT contains ALT allele ``alt_index``
    (1-based). depth = FORMAT/DP, else the sum of AD; allele_count = AD[alt], else
    FORMAT/AO (freebayes); VAF = AD[alt] / sum(AD) when AD is present, else
    FORMAT/AF, else allele_count / depth."""
    from acmg_classifier.models.variant import SampleCall

    if not samples:
        return []
    try:
        genotypes = variant.genotypes
    except Exception:  # noqa: BLE001 -- no GT field
        return []
    dp, ad, af = _format(variant, "DP"), _format(variant, "AD"), _format(variant, "AF")
    ao = _format(variant, "AO")
    calls = []
    for s, gt in enumerate(genotypes):
        if alt_index not in gt[:-1]:
            continue
        depth = _int_or_none(dp[s][0]) if dp is not None else None
        count = vaf = None
        if ad is not None:
            counts = [_int_or_none(x) for x in ad[s]]
            total = sum(c for c in counts if c is not None)
            if alt_index < len(counts) and counts[alt_index] is not None:
                count = counts[alt_index]
                if total > 0:
                    vaf = count / total
            if depth is None and total > 0:
                depth = total
        if count is None and ao is not None and alt_index - 1 < len(ao[s]):
            count = _int_or_none(ao[s][alt_index - 1])
        if vaf is None and af is not None and alt_index - 1 < len(af[s]):
            value = float(af[s][alt_index - 1])
            vaf = None if value != value else value          # NaN = missing
        if vaf is None and count is not None and depth:
            vaf = count / depth
        calls.append(SampleCall(sample=samples[s], depth=depth, allele_count=count,
                                vaf=None if vaf is None else round(vaf, 4)))
    return calls


def read_vcf(
    vcf_path: Path,
    assembly: Optional[Assembly] = None,
    sample_id: Optional[str] = None,
    fasta: Optional[Path] = None,
) -> Iterator[VariantRecord]:
    """Yield VariantRecord for every ALT allele in the VCF (multi-allelic split).

    Each allele is normalised (minimal representation, left-aligned against the
    reference ``fasta`` when given) so that it matches ClinVar / gnomAD / SpliceAI,
    which are looked up by exact CHROM:POS:REF:ALT. Splitting a multi-allelic
    record keeps the shared padding (``GAAA>GAA,G`` -> ``GAAA>GAA``), and such an
    allele would otherwise miss every database. The original representation is kept
    in ``VariantRecord.input_variant``. Unsupported alleles (ALT='.', symbolic,
    breakend, '*') are passed through untouched."""
    from acmg_classifier.annotation.normalizer import normalize_record, open_fasta
    from acmg_classifier.models.variant import unsupported_allele_reason
    resolved_assembly = assembly or detect_assembly_from_header(vcf_path)
    if resolved_assembly is None:
        raise ValueError(
            f"Cannot determine assembly for {vcf_path}. "
            "Pass --assembly explicitly."
        )

    ref_fasta = open_fasta(fasta)
    if fasta is not None and ref_fasta is None:
        log.warning("normalisation_fasta_unavailable", path=str(fasta),
                    hint="indels are trimmed but not left-aligned")
    vcf = VCF(str(vcf_path))
    samples = list(vcf.samples or [])
    n_normalised = 0       # alleles whose representation changed
    n_records = 0          # VCF data lines seen
    n_no_alt = 0           # records with empty ALT (reference / '.' sites)
    n_failed = 0           # ALT alleles that failed VariantRecord construction
    n_yielded = 0          # VariantRecords actually emitted
    try:
        for variant in vcf:
            n_records += 1
            # cyvcf2 returns the FILTER string (or None for PASS); don't char-join it.
            filter_val = variant.FILTER if variant.FILTER else "PASS"
            # ALT='.' (no alternate allele) -> cyvcf2 yields an empty ALT list.
            # Emit a placeholder row so the variant still appears in the output.
            alts = variant.ALT if variant.ALT else ["."]
            if variant.ALT == [] or variant.ALT is None:
                n_no_alt += 1
                # Expected for no-variant sites; the pipeline reports these in a
                # dedicated "not_annotated" file, so keep this at debug level.
                log.debug(
                    "vcf_record_no_alt",
                    chrom=variant.CHROM, pos=variant.POS, ref=variant.REF,
                )
            for alt_index, alt in enumerate(alts, start=1):
                try:
                    record = VariantRecord(
                        chrom=variant.CHROM,
                        pos=variant.POS,
                        ref=variant.REF,
                        alt=alt,
                        assembly=resolved_assembly,
                        vcf_id=variant.ID,
                        qual=variant.QUAL,
                        filter=filter_val,
                        sample_id=sample_id,
                    )
                except Exception as exc:
                    n_failed += 1
                    log.error(
                        "vcf_record_failed",
                        chrom=variant.CHROM, pos=variant.POS,
                        ref=variant.REF, alt=alt, error=str(exc),
                    )
                    continue
                if samples and variant.ALT:
                    try:
                        record.sample_calls = sample_calls_for(variant, samples, alt_index)
                    except Exception as exc:  # noqa: BLE001 -- genotypes are informative only
                        log.warning("vcf_sample_fields_failed", chrom=variant.CHROM,
                                    pos=variant.POS, error=str(exc))
                if unsupported_allele_reason(record.ref, record.alt) is None:
                    normalised = normalize_record(record, ref_fasta)
                    if normalised is not record:
                        n_normalised += 1
                        record = normalised
                n_yielded += 1
                yield record
    finally:
        vcf.close()
        if ref_fasta is not None:
            ref_fasta.close()
        log.info(
            "vcf_read_summary",
            records=n_records, no_alt=n_no_alt, failed=n_failed, yielded=n_yielded,
            normalised=n_normalised,
        )
