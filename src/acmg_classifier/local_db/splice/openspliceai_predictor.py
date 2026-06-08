"""OpenSpliceAI predictor (runtime splice scoring; GPL-3.0, commercial use OK).

Runs `openspliceai variant` as a subprocess on a temporary VCF.
Install: pip install openspliceai
Models:  download OSAI_MANE from https://github.com/Kuanhao-Chao/openspliceai
         and place under data/<assembly>/openspliceai/<flanking_size>nt/
"""
from __future__ import annotations
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.local_db.splice.base import SplicePredictor
from acmg_classifier.models.annotation import SpliceScore
from acmg_classifier.models.enums import Assembly
from acmg_classifier.models.variant import VariantRecord

log = structlog.get_logger()

# Standard chromosome lengths required for valid ##contig VCF headers.
# pysam (used internally by openspliceai) raises a BCF write error when
# contig definitions are missing from the input VCF header.
_CONTIG_LENGTHS: dict[Assembly, dict[str, int]] = {
    Assembly.GRCH38: {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
        "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
        "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
        "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
        "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
        "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
        "chrM": 16569,
    },
    Assembly.GRCH37: {
        "chr1": 249250621, "chr2": 243199373, "chr3": 198022430,
        "chr4": 191154276, "chr5": 180915260, "chr6": 171115067,
        "chr7": 159138663, "chr8": 146364022, "chr9": 141213431,
        "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
        "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
        "chr16": 90354753, "chr17": 81195210, "chr18": 78077248,
        "chr19": 59128983, "chr20": 63025520, "chr21": 48129895,
        "chr22": 51304566, "chrX": 155270560, "chrY": 59373566,
        "chrM": 16571,
    },
}

_ASSEMBLY_ANNOTATION = {
    Assembly.GRCH38: "grch38",
    Assembly.GRCH37: "grch37",
}


class OpenSpliceAIPredictor(SplicePredictor):
    """Runtime splice predictor using OpenSpliceAI.

    Runs `openspliceai variant` as a subprocess in batch mode for efficiency.
    Requires: `pip install openspliceai` and OSAI_MANE model files.

    Walker 2023 thresholds (same convention as SpliceAI):
      max_delta >= 0.2 → PP3 Moderate, max_delta <= 0.1 → BP4 Supporting.
    """

    def __init__(
        self,
        model_dir: Optional[Path],
        assembly: Assembly,
        flanking_size: int = 2000,
        ref_genome: Optional[Path] = None,
        annotation_file: Optional[Path] = None,
    ) -> None:
        self._model_dir = model_dir
        self._ref_genome = ref_genome
        # The `-A grch38` / `grch37` KEYWORDS are unusable: openspliceai's
        # Annotator maps them to the relative path "./data/vcf/<asm>.txt", which
        # only resolves if the CWD happens to be the openspliceai source tree.
        # We therefore pass an explicit annotation-table path (downloaded by
        # setup_data.step_openspliceai). The keyword is kept only for the
        # warning hint below.
        self._annotation_keyword = _ASSEMBLY_ANNOTATION[assembly]
        self._annotation_file = annotation_file
        self._flanking_size = flanking_size
        self._assembly = assembly
        self._cache: dict[str, SpliceScore] = {}
        self._available: bool = self._check_available()

        if not self._available:
            log.warning(
                "openspliceai_unavailable",
                model_dir=str(model_dir) if model_dir else None,
                ref_genome=str(ref_genome) if ref_genome else None,
                annotation_file=str(annotation_file) if annotation_file else None,
                hint=(
                    "OpenSpliceAI scoring will be skipped. Requires: `pip install "
                    "openspliceai`; OSAI_MANE model files under "
                    "data/<assembly>/openspliceai/<flanking_size>nt/; the reference "
                    "genome FASTA (`-R`); and the gene annotation table "
                    f"{self._annotation_keyword}.txt (`-A`) — both staged by "
                    "`setup_data.py`."
                ),
            )

    def _check_available(self) -> bool:
        if shutil.which("openspliceai") is None:
            return False
        if self._model_dir is None or not self._model_dir.exists():
            return False
        # `openspliceai variant` REQUIRES a reference genome FASTA (-R). Without
        # it the subprocess exits non-zero and every variant silently scores as
        # unavailable, so treat a missing reference as "tool unavailable".
        if self._ref_genome is None or not self._ref_genome.exists():
            return False
        # The annotation table (-A) is likewise required: the built-in keyword
        # resolves to a non-existent relative path, so without an explicit file
        # openspliceai crashes on every batch.
        if self._annotation_file is None or not self._annotation_file.exists():
            return False
        return True

    def is_available(self) -> bool:
        return self._available

    def precompute(self, variants: list[VariantRecord]) -> None:
        """Batch-predict splice scores for all variants via openspliceai variant."""
        if not self._available or not variants:
            return

        with tempfile.TemporaryDirectory() as tmpdir:
            in_vcf = Path(tmpdir) / "input.vcf"
            out_vcf = Path(tmpdir) / "output.vcf"

            _write_temp_vcf(in_vcf, variants, self._assembly)

            cmd = [
                "openspliceai", "variant",
                "-R", str(self._ref_genome),
                "-A", str(self._annotation_file),
                "-I", str(in_vcf),
                "-O", str(out_vcf),
                "-m", str(self._model_dir),
                "--model-type", "pytorch",
                "-f", str(self._flanking_size),
            ]
            try:
                subprocess.run(cmd, capture_output=True, text=True, check=True)
            except subprocess.CalledProcessError as exc:
                log.error(
                    "openspliceai_failed",
                    returncode=exc.returncode,
                    stderr=(exc.stderr or "")[-2000:],
                )
                return
            except FileNotFoundError:
                log.error("openspliceai_not_found")
                self._available = False
                return

            self._cache = _parse_output_vcf(out_vcf)

    def predict(self, variant: VariantRecord) -> SpliceScore:
        return self._cache.get(
            variant.key,
            SpliceScore(tool="openspliceai", is_available=False),
        )


def _write_temp_vcf(path: Path, variants: list[VariantRecord], assembly: Assembly) -> None:
    """Write a pysam-valid VCF with ##contig headers for all standard chromosomes.

    openspliceai uses pysam internally; without ##contig definitions pysam
    raises 'Invalid BCF, CONTIG id=N not present in the header' and the
    output write fails entirely.
    """
    contigs = _CONTIG_LENGTHS[assembly]
    with path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for chrom, length in contigs.items():
            f.write(f"##contig=<ID={chrom},length={length}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in variants:
            f.write(f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\t.\t.\n")


def _parse_output_vcf(path: Path) -> dict[str, SpliceScore]:
    """Parse openspliceai variant output VCF and extract max delta scores."""
    scores: dict[str, SpliceScore] = {}
    if not path.exists():
        return scores

    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            key = f"{chrom}:{pos}:{ref}:{alt}"
            max_delta = _extract_max_delta(fields[7])
            scores[key] = SpliceScore(
                tool="openspliceai",
                is_available=True,
                max_delta=max_delta,
                raw_score=max_delta,
            )
    return scores


def _extract_max_delta(info: str) -> Optional[float]:
    """Parse max delta score from an OpenSpliceAI INFO field.

    Format: OpenSpliceAI=ALLELE|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|...
    Multiple transcripts are comma-separated; we return the overall maximum
    of DS_AG/DS_AL/DS_DG/DS_DL — same convention as SpliceAI (Walker 2023).
    """
    for token in info.split(";"):
        if not token.startswith("OpenSpliceAI="):
            continue
        payload = token[len("OpenSpliceAI="):]
        best: Optional[float] = None
        for entry in payload.split(","):
            parts = entry.split("|")
            try:
                deltas = [float(parts[i]) for i in (2, 3, 4, 5) if i < len(parts)]
            except ValueError:
                continue
            if deltas:
                val = max(deltas)
                if best is None or val > best:
                    best = val
        return best
    return None
