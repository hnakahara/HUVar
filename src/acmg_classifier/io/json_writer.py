"""Write classification results to structured JSON."""
from __future__ import annotations
import json
import sys
from pathlib import Path
from typing import Optional

from acmg_classifier.models.classification import ClassificationResult


def write_json(results: list[ClassificationResult], output_path: Optional[Path]) -> None:
    data = [r.model_dump() for r in results]
    text = json.dumps(data, indent=2, default=str)
    if output_path:
        output_path.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text + "\n")
