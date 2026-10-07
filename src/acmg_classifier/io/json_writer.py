"""Write classification results to structured JSON."""
from __future__ import annotations
import json
import sys
from pathlib import Path
from typing import Optional

from acmg_classifier.models.classification import ClassificationResult


def write_json(results: list[ClassificationResult], output_path: Optional[Path],
               provenance: Optional[dict] = None) -> None:
    """Write results as JSON. With ``provenance`` (see acmg_classifier.provenance),
    the output is ``{"metadata": provenance, "results": [...]}``."""
    data = [r.model_dump() for r in results]
    if provenance is not None:
        data = {"metadata": provenance, "results": data}
    text = json.dumps(data, indent=2, default=str)
    if output_path:
        output_path.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text + "\n")
