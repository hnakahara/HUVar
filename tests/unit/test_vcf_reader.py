"""Unit tests for VCF reader."""
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

FIXTURE_VCF = Path(__file__).parent.parent / "fixtures" / "sample.vcf"


def test_fixture_vcf_exists():
    assert FIXTURE_VCF.exists()


def test_detect_assembly_grch38():
    from acmg_classifier.io.vcf_reader import detect_assembly_from_header
    from acmg_classifier.models.enums import Assembly
    result = detect_assembly_from_header(FIXTURE_VCF)
    assert result == Assembly.GRCH38
