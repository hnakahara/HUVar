"""Data provenance manifest / sidecar (Genome Medicine revision, Reviewer 3 comment 1)."""
import datetime as dt
import gzip
import json

from acmg_classifier.config import Config
from acmg_classifier.models.enums import Assembly
from acmg_classifier import provenance


def _mk(tmp_path, date="20260530"):
    d = tmp_path / "data"
    cv = d / "GRCh38" / "clinvar"
    cv.mkdir(parents=True)
    with gzip.open(cv / "clinvar_GRCh38.vcf.gz", "wt") as fh:
        fh.write(f"##fileformat=VCFv4.1\n##fileDate={date}\n##source=ClinVar\n#CHROM\tPOS\n")
    with gzip.open(cv / "ClinVarRCVRelease.xml.gz", "wt") as fh:
        fh.write('<?xml version="1.0"?>\n<ReleaseSet Type="full" Dated="2026-05-30">\n')
    return d


def test_collect_detects_clinvar_and_manifest(tmp_path):
    d = _mk(tmp_path)
    provenance.record_resource(d, "clingen_dosage", source="https://x", downloaded="2026-06-01")
    provenance.record_resource(d, "clinvar_vcf_GRCh38", source="https://y", downloaded="2026-06-01")
    prov = provenance.collect(Config(data_dir=d, assembly=Assembly.GRCH38))
    r = prov["resources"]
    assert r["clinvar_vcf"]["release"] == "2026-05-30"
    assert r["clinvar_vcf"]["downloaded"] == "2026-06-01"
    assert r["clinvar_rcv_xml"]["release"] == "2026-05-30"
    assert r["clingen_dosage"]["source"] == "https://x"
    assert "clinvar_vcf_GRCh38" not in r
    assert prov["assembly"] == "GRCh38" and prov["huvar_version"]


def test_stale_clinvar_warning(tmp_path):
    old = (dt.date.today() - dt.timedelta(days=200)).strftime("%Y%m%d")
    d = _mk(tmp_path, old)
    prov = provenance.collect(Config(data_dir=d, assembly=Assembly.GRCH38))
    assert prov["warnings"] and "days old" in prov["warnings"][0]


def test_sidecar(tmp_path):
    d = _mk(tmp_path)
    out = tmp_path / "res.tsv"
    out.write_text("x")
    side = provenance.write_sidecar(Config(data_dir=d, assembly=Assembly.GRCH38), out)
    assert side.name == "res.tsv.provenance.json"
    assert json.loads(side.read_text())["resources"]["clinvar_vcf"]["release"] == "2026-05-30"
