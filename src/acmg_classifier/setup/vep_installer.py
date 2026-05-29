"""Detect or help install Ensembl VEP (Docker / conda / native)."""
from __future__ import annotations
import shutil
import subprocess

import structlog

log = structlog.get_logger()

_VEP_MIN_VERSION = 111


def find_vep_cmd() -> str:
    """Return the VEP command string (vep path or docker run invocation).

    Native install is preferred because docker startup adds ~1-2s per VEP
    invocation, which dominates batch annotation throughput. Docker is only
    used as a fallback so users who can't install VEP system-wide are still
    able to run the pipeline."""
    native = shutil.which("vep")
    if native:
        log.info("vep_found", method="native", path=native)
        return native

    # Docker fallback: verify the specific release image exists locally
    # before returning "docker" so a missing image surfaces a clear pull
    # instruction rather than a cryptic docker-run failure later.
    if shutil.which("docker"):
        try:
            result = subprocess.run(
                ["docker", "image", "inspect", "ensemblorg/ensembl-vep:release_111"],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                log.info("vep_found", method="docker")
                return "docker"
        except Exception:
            pass

    raise RuntimeError(
        "VEP not found. Install via one of:\n"
        "  docker pull ensemblorg/ensembl-vep:release_111\n"
        "  conda install -c bioconda ensembl-vep=111\n"
        "  See: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html"
    )


def check_vep_version(vep_cmd: str) -> bool:
    """Verify VEP version is >= 111."""
    try:
        result = subprocess.run([vep_cmd, "--help"], capture_output=True, text=True)
        for line in result.stdout.splitlines():
            if "version" in line.lower() and str(_VEP_MIN_VERSION) in line:
                return True
    except Exception:
        pass
    return False
