#!/usr/bin/env python3
"""Download every ClinGen CSpec criteria specification (GNxxx) from the CSpec API
and keep it as a reviewable snapshot.

    https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/GNxxx

Each specification is saved byte-for-byte as returned by the API
(``<out>/GNxxx.json``). ``manifest.tsv`` lists every saved specification (title,
VCEP, status, version, last update, genes, SHA-256) and how it compares with the
June 2026 snapshot the engine tables were built from
(``resources/clingen/cspec_json``): ``new`` / ``version`` (version changed) /
``content`` (same version, different criteria) / ``same``.

The API has no listing endpoint, so ids are tried from GN001 upwards until
``--stop-after`` consecutive ids fail (non-existent ids return HTTP 500/404).

Usage:
    python scripts/fetch_cspec_api.py                       # refresh resources/clingen/cspec_api
    python scripts/fetch_cspec_api.py --out DIR --workers 4
    python scripts/fetch_cspec_api.py --manifest-only       # rebuild manifest.tsv from DIR
    python scripts/fetch_cspec_api.py --pdf                 # + PDFs of the Released specs
    python scripts/fetch_cspec_api.py --pdf-only [--pdf-ids GN043,GN129] [--pdf-all]

Note: the API JSON carries the per-strength criteria text only. The general
"VCEP Specifications" paragraph and the comments shown on the CSpec web page
(e.g. the RASopathy PM1/PM5 rule) are not part of it. They are in the PDF the
web page offers ("PDF" button), which ``--pdf`` downloads from the same service
the page uses (``/cspec/ui/svi/pdf?url=<page>?layout=print``) into
``resources/clingen/cspec_pdf/GNxxx.pdf`` with its own ``manifest.tsv``. A PDF
is rendered on the server (about 10 s each), so only Released specifications
are fetched by default, and a PDF whose specification version is unchanged is
not downloaded again.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
API = "https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/{}"
DEFAULT_OUT = ROOT / "resources" / "clingen" / "cspec_api"
DEFAULT_PDF_OUT = ROOT / "resources" / "clingen" / "cspec_pdf"
UI_PAGE = "https://cspec.genome.network/cspec/ui/svi/doc/{}?layout=print"
PDF_SERVICE = "https://cspec.genome.network/cspec/ui/svi/pdf"
SNAPSHOT = ROOT / "resources" / "clingen" / "cspec_json"


def _gn(n: int) -> str:
    return f"GN{n:03d}"


def _fetch(gn: str, retries: int = 3, timeout: int = 60) -> tuple[str, int, bytes | None]:
    """(gn, http status, body or None). A body is returned only for a JSON
    specification whose @id names this GN id."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(API.format(gn), headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read()
            doc = json.loads(body)
            if str(doc.get("@id", "")).endswith("/" + gn):
                return gn, 200, body
            return gn, 200, None
        except urllib.error.HTTPError as exc:
            if exc.code in (404, 500):          # the API answers 500 for unknown ids
                return gn, exc.code, None
            status = exc.code
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            status = -1
            print(f"  [retry] {gn}: {exc}", file=sys.stderr)
        time.sleep(2 * (attempt + 1))
    return gn, status, None


def fetch_all(out: Path, workers: int, stop_after: int, first: int) -> list[str]:
    out.mkdir(parents=True, exist_ok=True)
    saved: list[str] = []
    n, misses = first, 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        while misses < stop_after:
            batch = [_gn(i) for i in range(n, n + workers * 2)]
            n += len(batch)
            for gn, status, body in pool.map(_fetch, batch):
                if body is None:
                    misses += 1
                    if status not in (404, 500):
                        print(f"  [warn] {gn}: HTTP {status}", file=sys.stderr)
                    continue
                misses = 0
                (out / f"{gn}.json").write_bytes(body)
                saved.append(gn)
    return saved


def _genes(doc: dict) -> list[str]:
    out = []
    for rs in doc.get("ruleSets") or []:
        for g in rs.get("genes") or []:
            label = g.get("label")
            if label and label not in out:
                out.append(label)
    return out


def _criteria(doc: dict) -> list:
    """Order-independent view of the criteria texts, for the snapshot comparison."""
    rows = []
    for rs in doc.get("ruleSets") or []:
        genes = ",".join(sorted(g.get("label", "") for g in rs.get("genes") or []))
        for cc in rs.get("criteriaCodes") or []:
            for es in cc.get("evidenceStrengths") or []:
                rows.append((genes, cc.get("label"), es.get("label"),
                             (es.get("applicability") or "").lower(),
                             (es.get("description") or "").strip()))
    return sorted(rows)


def _version(v) -> tuple:
    """Version as an int tuple without trailing zeros: the API now writes "2.0"
    where the June snapshot had "2.0.0"."""
    parts = []
    for x in str(v or "").split("."):
        try:
            parts.append(int(x))
        except ValueError:
            parts.append(x)
    while parts and parts[-1] == 0:
        parts.pop()
    return tuple(parts)


def build_manifest(out: Path) -> Path:
    rows = []
    for path in sorted(out.glob("GN*.json")):
        raw = path.read_bytes()
        doc = json.loads(raw)
        gn = path.stem
        snap = SNAPSHOT / path.name
        if not snap.exists():
            diff = "new"
        else:
            old = json.loads(snap.read_bytes())
            if _version(old.get("version")) != _version(doc.get("version")):
                diff = f"version {old.get('version')} -> {doc.get('version')}"
            elif _criteria(old) != _criteria(doc):
                diff = "content"
            else:
                diff = "same"
        n_codes = sum(len(rs.get("criteriaCodes") or []) for rs in doc.get("ruleSets") or [])
        rows.append({
            "gn": gn,
            "genes": ",".join(_genes(doc)),
            "status": doc.get("cspecStatus") or doc.get("currentStatus") or "",
            "version": doc.get("version") or "",
            "last_updated": (doc.get("lastUpdated") or "")[:10],
            "vcep": (doc.get("affiliation") or {}).get("label", ""),
            "n_criteria": n_codes,
            "vs_2026-06-05": diff,
            "title": doc.get("label") or "",
            "sha256": hashlib.sha256(raw).hexdigest(),
        })
    man = out / "manifest.tsv"
    with man.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["gn"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return man


def _fetch_pdf(gn: str, retries: int = 3, timeout: int = 150) -> tuple[str, bytes | None, str]:
    from urllib.parse import urlencode
    url = PDF_SERVICE + "?" + urlencode({"url": UI_PAGE.format(gn), "filename": f"ClinGen_{gn}"})
    err = ""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                body = resp.read()
            if body[:5] == b"%PDF-":
                return gn, body, ""
            err = "not a PDF"
        except (urllib.error.URLError, TimeoutError) as exc:
            err = str(exc)
        time.sleep(5 * (attempt + 1))
    return gn, None, err


def _pdf_texts(out: Path) -> None:
    """Plain-text copies (``text/GNxxx.txt``, poppler ``pdftotext``) for grep and
    ``git diff``; refreshed when the PDF is newer. Skipped without pdftotext."""
    import shutil
    import subprocess
    exe = shutil.which("pdftotext")
    if exe is None:
        print("  [info] pdftotext not found; text copies skipped")
        return
    tdir = out / "text"
    tdir.mkdir(exist_ok=True)
    n = 0
    for pdf in sorted(out.glob("GN*.pdf")):
        txt = tdir / f"{pdf.stem}.txt"
        if txt.exists() and txt.stat().st_mtime >= pdf.stat().st_mtime:
            continue
        subprocess.run([exe, "-enc", "UTF-8", str(pdf), str(txt)], check=False)
        n += 1
    print(f"  text copies written: {n}")


def fetch_pdfs(api_dir: Path, out: Path, ids: list[str] | None, all_status: bool,
               workers: int, budget: float | None) -> None:
    """Download the web-page PDFs into *out*; resumable (``manifest.tsv`` keeps
    the version each PDF was made from) and stoppable after *budget* seconds."""
    man = {r["gn"]: r for r in csv.DictReader(
        (api_dir / "manifest.tsv").open(encoding="utf-8"), delimiter="\t")}
    if ids is None:
        ids = [gn for gn, r in man.items() if all_status or r["status"] == "Released"]
    out.mkdir(parents=True, exist_ok=True)
    pman_path = out / "manifest.tsv"
    done: dict[str, dict] = {}
    if pman_path.exists():
        done = {r["gn"]: r for r in csv.DictReader(pman_path.open(encoding="utf-8"),
                                                    delimiter="\t")}
    todo = [gn for gn in ids if gn in man and not (
        (out / f"{gn}.pdf").exists() and gn in done
        and _version(done[gn]["version"]) == _version(man[gn]["version"]))]
    print(f"PDFs: {len(ids)} selected, {len(todo)} to download")
    t0 = time.time()
    failed = []

    def _write_manifest() -> None:
        cols = ["gn", "genes", "status", "version", "fetched_at", "bytes", "sha256", "title"]
        with pman_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            for gn in sorted(done):
                w.writerow({c: done[gn].get(c, "") for c in cols})

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i in range(0, len(todo), workers):
            if budget is not None and time.time() - t0 > budget:
                print(f"  time budget reached; {len(todo) - i} PDFs left (run again to resume)")
                break
            for gn, body, err in pool.map(_fetch_pdf, todo[i:i + workers]):
                if body is None:
                    failed.append(gn)
                    print(f"  [warn] {gn}: {err}", file=sys.stderr)
                    continue
                (out / f"{gn}.pdf").write_bytes(body)
                r = man[gn]
                done[gn] = {"gn": gn, "genes": r["genes"], "status": r["status"],
                            "version": r["version"],
                            "fetched_at": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
                            "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(),
                            "title": r["title"]}
            _write_manifest()
    _write_manifest()
    _pdf_texts(out)
    if failed:
        print(f"  failed: {', '.join(failed)}")
    print(f"PDFs in {out}: {len(done)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--first", type=int, default=1, help="first GN number to try")
    ap.add_argument("--stop-after", type=int, default=20,
                    help="stop after this many consecutive missing ids")
    ap.add_argument("--manifest-only", action="store_true")
    ap.add_argument("--pdf", action="store_true", help="also download the web-page PDFs")
    ap.add_argument("--pdf-only", action="store_true", help="only the PDFs (uses the manifest)")
    ap.add_argument("--pdf-out", type=Path, default=DEFAULT_PDF_OUT)
    ap.add_argument("--pdf-ids", help="comma-separated GN ids instead of the Released ones")
    ap.add_argument("--pdf-all", action="store_true", help="PDFs of every status")
    ap.add_argument("--pdf-workers", type=int, default=2)
    ap.add_argument("--pdf-budget", type=float, default=None,
                    help="stop starting new PDFs after this many seconds")
    args = ap.parse_args()
    if args.pdf_only:
        ids = [x.strip() for x in args.pdf_ids.split(",")] if args.pdf_ids else None
        fetch_pdfs(args.out, args.pdf_out, ids, args.pdf_all, args.pdf_workers, args.pdf_budget)
        return
    if not args.manifest_only:
        started = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        saved = fetch_all(args.out, args.workers, args.stop_after, args.first)
        (args.out / "FETCHED_AT").write_text(
            f"{started}\t{len(saved)} specifications\t{API.format('GNxxx')}\n", encoding="utf-8")
        print(f"saved {len(saved)} specifications ({saved[0]}..{saved[-1]}) -> {args.out}")
    man = build_manifest(args.out)
    print(f"manifest -> {man}")
    if args.pdf:
        ids = [x.strip() for x in args.pdf_ids.split(",")] if args.pdf_ids else None
        fetch_pdfs(args.out, args.pdf_out, ids, args.pdf_all, args.pdf_workers, args.pdf_budget)


if __name__ == "__main__":
    main()
