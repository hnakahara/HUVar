"""Data provenance: which versions of the annotation resources a run used.

Two sources of information are combined:

* ``<data_dir>/data_manifest.json`` — written by ``scripts/setup_data.py`` as
  each resource is downloaded (source URL, download date, optional version);
* auto-detection from the files themselves (ClinVar VCF ``##fileDate``, ClinVar
  RCV XML ``Dated``, VEP cache ``info.txt``, the gnomAD build names used by
  :class:`~acmg_classifier.config.Config`, CSpec snapshot date, MANE release,
  eRepo supplement), so that existing installations without a manifest still
  report their versions.

The result is shown by ``acmg-classify status``, written next to every
classification output as ``<output>.provenance.json`` and embedded in JSON
output, so each classification can be traced to the data it was based on.
"""
from __future__ import annotations

import datetime as _dt
import gzip
import json
import re
from pathlib import Path
from typing import Any, Optional

MANIFEST_NAME = "data_manifest.json"
# ClinVar older than this (days, relative to today) triggers a warning.
CLINVAR_MAX_AGE_DAYS = 90

# parents: [0]=acmg_classifier [1]=src [2]=repo
_REPO = Path(__file__).resolve().parents[2]


def software_version() -> str:
    try:
        from importlib.metadata import version
        return version("acmg-classifier")
    except Exception:  # noqa: BLE001
        return "unknown"


def _today() -> str:
    return _dt.date.today().isoformat()


def load_manifest(data_dir: Path) -> dict[str, Any]:
    p = Path(data_dir) / MANIFEST_NAME
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {}


def record_resource(data_dir: Path, resource: str, path: Optional[Path] = None, **fields) -> None:
    """Add/replace one resource entry in the data manifest."""
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    man = load_manifest(data_dir)
    entry = {k: (str(v) if isinstance(v, Path) else v) for k, v in fields.items()}
    entry.setdefault("downloaded", _today())
    if path is not None:
        entry["path"] = str(path)
    man[resource] = entry
    (data_dir / MANIFEST_NAME).write_text(json.dumps(man, indent=2, sort_keys=True), encoding="utf-8")


# ---------------------------------------------------------------------------
# Auto-detection helpers (cheap: read only file headers)
# ---------------------------------------------------------------------------

def _clinvar_vcf_date(path: Path) -> Optional[str]:
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for _ in range(50):
                line = fh.readline()
                if not line.startswith("##"):
                    break
                if line.startswith("##fileDate="):
                    d = line.split("=", 1)[1].strip()
                    if re.fullmatch(r"\d{8}", d):
                        d = f"{d[:4]}-{d[4:6]}-{d[6:]}"
                    return d
    except Exception:  # noqa: BLE001
        return None
    return None


def _clinvar_xml_date(path: Path) -> Optional[str]:
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            head = fh.read(4096)
        m = re.search(r'Dated="([0-9-]+)"', head)
        return m.group(1) if m else None
    except Exception:  # noqa: BLE001
        return None


def _vep_info(cache_dir: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    infos = sorted(cache_dir.glob("*/*/info.txt")) if cache_dir.exists() else []
    for info in infos:
        try:
            for line in info.read_text(encoding="utf-8", errors="replace").splitlines():
                parts = line.split("\t", 1)
                if len(parts) == 2 and parts[0] in (
                    "source_gencode", "source_refseq", "source_assembly", "assembly",
                ):
                    out[parts[0]] = parts[1].strip()
            out["cache"] = f"{info.parent.parent.name}/{info.parent.name}"
            break
        except Exception:  # noqa: BLE001
            continue
    return out


def _cspec_snapshot() -> Optional[str]:
    p = _REPO / "resources" / "clingen" / "cspec_json" / "cspec_summary.json"
    try:
        meta = json.loads(p.read_text(encoding="utf-8")).get("metadata", {})
        return (meta.get("rendered") or {}).get("when")
    except Exception:  # noqa: BLE001
        return None


def _bundled_versions() -> dict[str, dict[str, str]]:
    """Versions of curated resources shipped in resources/ (tracked in git)."""
    p = _REPO / "resources" / "shared" / "bundled_data_versions.tsv"
    out: dict[str, dict[str, str]] = {}
    try:
        lines = [ln for ln in p.read_text(encoding="utf-8").splitlines()
                 if ln and not ln.startswith("#")]
        header = lines[0].split("\t")
        for ln in lines[1:]:
            row = dict(zip(header, ln.split("\t")))
            out[row.pop("resource")] = row
    except Exception:  # noqa: BLE001
        pass
    return out


def _first_comment(path: Path) -> Optional[str]:
    try:
        with path.open(encoding="utf-8") as fh:
            line = fh.readline().strip()
        return line.lstrip("#").strip() if line.startswith("#") else None
    except Exception:  # noqa: BLE001
        return None


def collect(cfg) -> dict[str, Any]:
    """Provenance record for a run configured by ``cfg``."""
    data_dir = Path(cfg.data_dir)
    man = load_manifest(data_dir)
    res: dict[str, Any] = {}
    asm = cfg.assembly.value
    used: set[str] = set()

    def put(name: str, **kv) -> None:
        key = f"{name}_{asm}" if f"{name}_{asm}" in man else name
        used.add(key)
        entry = dict(man.get(key, {}))
        entry.update({k: v for k, v in kv.items() if v is not None})
        if entry:
            res[name] = entry

    vcf_date = _clinvar_vcf_date(cfg.clinvar_vcf) if Path(cfg.clinvar_vcf).exists() else None
    xml = cfg.assembly_dir / "clinvar" / "ClinVarRCVRelease.xml.gz"
    put("clinvar_vcf", release=vcf_date, file=Path(cfg.clinvar_vcf).name)
    put("clinvar_rcv_xml", release=_clinvar_xml_date(xml) if xml.exists() else None)
    put("vep", release=str(getattr(cfg, "ensembl_release", "")) or None, **_vep_info(cfg.vep_cache_dir))
    put("gnomad", build=Path(cfg.gnomad_duckdb).name)
    nc = getattr(cfg, "gnomad_noncancer_duckdb", None)
    if nc is not None and Path(nc).exists():
        put("gnomad_non_cancer", build=Path(nc).name)
    bundled = _bundled_versions()
    for name, row in bundled.items():
        put(name, **row)
    put("clingen_cspec", rendered=_cspec_snapshot())
    sup = data_dir / "shared" / f"erepo_manual_criteria_{'hg38' if cfg.assembly.value == 'GRCh38' else 'hg19'}.tsv"
    put("clingen_erepo_supplement", file=sup.name if sup.exists() else None)
    for key, attr in (("clingen_dosage", "clingen_dosage_tsv"),
                      ("clingen_gene_validity", "clingen_gene_validity_csv")):
        p = getattr(cfg, attr, None)
        put(key, file=Path(p).name if p is not None and Path(p).exists() else None)
    pp2 = getattr(cfg, "pp2_gene_stats_tsv", None)
    if isinstance(pp2, Path) and pp2.exists():
        with pp2.open(encoding="utf-8") as fh:
            head = fh.readline().lstrip("# ").strip()
        put("pp2_gene_stats", file=pp2.name, note=head or None)
    put("run_options",
        exclude_self_expert_panel=bool(getattr(cfg, "exclude_self_expert_panel", False) is True),
        pm1_heuristic=bool(getattr(cfg, "pm1_heuristic", False) is True),
        tp53_functional=bool(getattr(cfg, "use_tp53_functional", False) is True),
        pm1_window=getattr(cfg, "pm1_window", None),
        pm1_min_path_variants=getattr(cfg, "pm1_min_path_variants", None),
        pp2_min_path=getattr(cfg, "pp2_min_path", None),
        pp2_max_benign_frac=getattr(cfg, "pp2_max_benign_frac", None))
    put("insilico", missense=getattr(cfg.insilico_tool, "value", str(cfg.insilico_tool)),
        splice=getattr(cfg.splice_tool, "value", str(cfg.splice_tool)))
    other_asm = "GRCh37" if asm == "GRCh38" else "GRCh38"
    for k, v in man.items():  # anything else recorded at setup (this assembly)
        if k not in used and not k.endswith(f"_{other_asm}"):
            res.setdefault(k, v)

    warnings = []
    rel = res.get("clinvar_vcf", {}).get("release")
    if rel:
        try:
            age = (_dt.date.today() - _dt.date.fromisoformat(rel[:10])).days
            res["clinvar_vcf"]["age_days"] = age
            if age > CLINVAR_MAX_AGE_DAYS:
                warnings.append(
                    f"ClinVar release {rel} is {age} days old (> {CLINVAR_MAX_AGE_DAYS}); "
                    "consider `setup_data.py --force-clinvar --only clinvar-vcf clinvar-sqlite`"
                )
        except ValueError:
            pass
    return {
        "huvar_version": software_version(),
        "run_date": _today(),
        "assembly": cfg.assembly.value,
        "resources": res,
        "warnings": warnings,
    }


def write_sidecar(cfg, output_path: Path) -> Path:
    """Write ``<output>.provenance.json`` next to a classification output."""
    prov = collect(cfg)
    side = output_path.with_name(output_path.name + ".provenance.json")
    side.write_text(json.dumps(prov, indent=2, sort_keys=True), encoding="utf-8")
    return side
