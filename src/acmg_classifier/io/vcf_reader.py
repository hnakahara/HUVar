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


def read_vcf(
    vcf_path: Path,
    assembly: Optional[Assembly] = None,
    sample_id: Optional[str] = None,
) -> Iterator[VariantRecord]:
    """Yield VariantRecord for every ALT allele in the VCF (multi-allelic split)."""
    resolved_assembly = assembly or detect_assembly_from_header(vcf_path)
    if resolved_assembly is None:
        raise ValueError(
            f"Cannot determine assembly for {vcf_path}. "
            "Pass --assembly explicitly."
        )

    vcf = VCF(str(vcf_path))
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
            for alt in alts:
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
                n_yielded += 1
                yield record
    finally:
        vcf.close()
        log.info(
            "vcf_read_summary",
            records=n_records, no_alt=n_no_alt, failed=n_failed, yielded=n_yielded,
        )
