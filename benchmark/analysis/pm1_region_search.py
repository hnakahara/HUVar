#!/usr/bin/env python3
"""Reviewer 3, comment 7: calibrate the ClinVar PM1 hotspot heuristic against the
VCEP-defined PM1 regions / residues.

For every gene in pm1_hotspots.tsv, each residue of the MANE Select protein is
labelled positive when it lies in a VCEP PM1 region or residue (77 genes) and
negative otherwise; genes whose VCEP declares PM1 not applicable (48 genes)
contribute only negatives. The heuristic marks residue p as a hotspot when the
ClinVar P/LP missense changes (distinct amino-acid changes, >= `stars`) within
+/- W residues reach N and the benign rule is satisfied:

  W        in {0, 1, 3, 5, 10, 15, 25}
  N        in {1, 2, 3, 4, 5, 6}
  benign   block (no B/LB missense in the window, current rule) |
           frac<=0.1 | frac<=0.2 of B/LB among missense in the window | ignore
  stars    1 | 2
  enrich   off | window P/LP density >= 3x or 5x the gene's mean density
  self     include | exclude the changes at residue p itself
           ("exclude" mimics a variant at a residue with no other report)

Residue-level precision / recall / F1 are pooled over all genes (the N/A genes
add false positives). Thresholds are selected on a random half of the genes
(stratified: defined / N/A) and evaluated on the other half, 200 times; the
configuration selected most often is then applied, variant by variant (query
variant excluded, as in HUVar), to the eRepo missense variants in the genes
where HUVar actually uses the heuristic (no VCEP PM1 table; independent of the
calibration genes) and to all eRepo PM1-eligible genes.

Run on the analysis server:
  python pm1_region_search.py --data-dir $DATA \
     --mane-exons <HUVar>/resources/shared/mane_exons.tsv \
     --hotspots <HUVar>/resources/shared/pm1_hotspots.tsv \
     --results erepo_benchmark/inputs/erepo_hg38_huhvar_default_results.tsv \
     --erepo <HUVar>/resources/GRCh38/erepo_classifications_hg38.tsv \
     --out results/pm1region
"""
from __future__ import annotations

import argparse
import csv
import itertools
import random
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from acmg_classifier.config import Config
from acmg_classifier.criteria.pm1_hotspots import PM1Hotspots, _parse_ranges, _parse_residues
from acmg_classifier.local_db.clinvar_sqlite import _PP2_BENIGN, _PP2_PATH, _is_missense_p

WINDOWS = (0, 1, 3, 5, 10, 15, 25)
COUNTS = (1, 2, 3, 4, 5, 6)
BENIGN = ("block", "frac0.1", "frac0.2", "ignore")
STARS = (1, 2)
ENRICH = (0, 3, 5)
SELF = ("include", "exclude")
DEFAULT = (25, 3, "block", 1, 0, "include")
_MIS = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?$")


def strip_chr(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def cfg_name(c):
    w, n, b, s, e, sf = c
    return f"W={w};N={n};benign={b};stars>={s};enrich={'off' if not e else str(e) + 'x'};self={sf}"


def prf(tp, fp, fn):
    p = tp / (tp + fp) if tp + fp else 0.0
    r = tp / (tp + fn) if tp + fn else 0.0
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)


def load_lengths(path: Path) -> dict[str, int]:
    out = {}
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            if r.get("mane_tag") != "MANE Select":
                continue
            out[r["gene"]] = (int(r["cds_end"]) - int(r["cds_start"]) + 1) // 3 - 1
    return out


def load_labels(path: Path, lengths):
    """{gene: bool array (index = residue, 1-based; [0] unused)}; N/A -> all False."""
    raw = defaultdict(list)
    na = set()
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["gene_symbol"].strip()
            st = (r.get("strength") or "").strip().lower()
            if st == "not_applicable":
                na.add(g)
            elif st != "cys_creating":
                raw[g].append((_parse_ranges(r.get("regions") or ""), _parse_residues(r.get("residues") or "")))
    labels, kind = {}, {}
    for g in set(raw) | na:
        if g not in lengths:
            continue
        L = lengths[g]
        lab = np.zeros(L + 1, dtype=bool)
        for ranges, residues in raw.get(g, ()):
            for a, b in ranges:
                lab[max(1, a):min(L, b) + 1] = True
            for p in residues:
                if 1 <= p <= L:
                    lab[p] = True
        if g in raw and not lab.any():
            continue
        labels[g] = lab
        kind[g] = "defined" if g in raw else "not_applicable"
    return labels, kind


def gene_rows(con, gene):
    q = ("SELECT chrom, pos, ref, alt, hgvs_p, amino_acid_change, clinical_significance, star_rating, "
         "codon_position FROM variants WHERE gene_symbol = ? AND star_rating >= 1 "
         "AND clinical_significance IN (?,?,?,?,?,?)")
    out = []
    for r in con.execute(q, (gene, *_PP2_PATH, *_PP2_BENIGN)).fetchall():
        if _is_missense_p(r[4]) and r[8] is not None:
            out.append((strip_chr(str(r[0])), r[1], r[2], r[3], r[5] or r[4], r[6] in _PP2_PATH, r[7], int(r[8])))
    return out


def residue_counts(rows, L, stars):
    P = np.zeros(L + 2)
    B = np.zeros(L + 2)
    seen = set()
    for chrom, pos, ref, alt, key, path, st, cod in rows:
        if st < stars or not 1 <= cod <= L or (key, path) in seen:
            continue
        seen.add((key, path))
        (P if path else B)[cod] += 1
    return P, B


def window_sum(x, w):
    c = np.concatenate([[0.0], np.cumsum(x)])
    n = len(x)
    idx = np.arange(n)
    lo = np.clip(idx - w, 0, n)
    hi = np.clip(idx + w + 1, 0, n)
    return c[hi] - c[lo]


def predict_residues(P, B, L, conf):
    w, n, benign, _, enrich, self_ = conf
    Pw, Bw = window_sum(P, w), window_sum(B, w)
    if self_ == "exclude":
        Pw, Bw = Pw - P, Bw - B
    ok = Pw >= n
    if benign == "block":
        ok &= Bw == 0
    elif benign.startswith("frac"):
        f = float(benign[4:])
        tot = Pw + Bw
        ok &= np.where(tot > 0, Bw / np.maximum(tot, 1), 0) <= f
    if enrich:
        mean = P[1:L + 1].sum() / L if L else 0
        ok &= Pw / (2 * w + 1) >= enrich * mean if mean > 0 else False
    pred = ok[: L + 1].copy()
    pred[0] = False
    return pred


def variant_eval(con, results, erepo, hot, conf, cache, lengths):
    """Variant-level evaluation of one configuration on eRepo missense variants."""
    w, n, benign, stars, enrich, _ = conf
    truth = {}
    with erepo.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            truth[r["variant_id"]] = "PM1" in {c.strip().split("_")[0] for c in (r.get("met_codes") or "").split(",")}
    cnt = {"heuristic_genes": Counter(), "all_pm1_genes": Counter()}
    with results.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = _MIS.match(r.get("hgvs_p") or "")
            vid, gene = r["variant_id"], r.get("gene") or ""
            if not m or m.group(3) == "Ter" or m.group(1) == m.group(3) or vid not in truth or not gene:
                continue
            if hot.is_not_applicable(gene):
                continue
            heur = not (hot.has_gene(gene) or hot.has_cys_creating(gene))
            if gene not in cache:
                cache[gene] = gene_rows(con, gene)
            rows = cache[gene]
            aa = int(m.group(2))
            chrom, pos, ref, alt = vid.split(":")
            chrom, pos = strip_chr(chrom), int(pos)
            pk, bk, allp = set(), set(), set()
            for rc, rp, rr, ra, key, path, st, cod in rows:
                if st < stars:
                    continue
                if path:
                    allp.add(key)
                if abs(cod - aa) > w:
                    continue
                if rc == chrom and rp == pos and rr == ref and ra == alt:
                    continue
                (pk if path else bk).add(key)
            npw, nbw = len(pk), len(bk)
            ok = npw >= n
            if benign == "block":
                ok = ok and nbw == 0
            elif benign.startswith("frac"):
                ok = ok and (nbw / (npw + nbw) if npw + nbw else 0) <= float(benign[4:])
            if enrich and ok:
                L = lengths.get(gene) or max([c for *_, c in rows] + [aa])
                ok = npw / (2 * w + 1) >= enrich * (len(allp) / L)
            t = truth[vid]
            k = "tp" if ok and t else "fp" if ok else "fn" if t else "tn"
            cnt["all_pm1_genes"][k] += 1
            if heur:
                cnt["heuristic_genes"][k] += 1
    out = {}
    for grp, c in cnt.items():
        p, r_, f = prf(c["tp"], c["fp"], c["fn"])
        out[grp] = dict(n=sum(c.values()), erepo_pm1=c["tp"] + c["fn"], predicted=c["tp"] + c["fp"],
                        precision=round(p, 3), sensitivity=round(r_, 3), f1=round(f, 3))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--mane-exons", type=Path, required=True)
    ap.add_argument("--hotspots", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--erepo", type=Path, required=True)
    ap.add_argument("--out", default="pm1region")
    ap.add_argument("--repeats", type=int, default=200)
    a = ap.parse_args()

    cfg = Config(data_dir=a.data_dir, assembly="GRCh38")
    con = sqlite3.connect("file:" + cfg.clinvar_sqlite.as_posix() + "?mode=ro&immutable=1", uri=True)
    lengths = load_lengths(a.mane_exons)
    labels, kind = load_labels(a.hotspots, lengths)
    genes = sorted(labels)
    cache = {g: gene_rows(con, g) for g in genes}
    counts = {(g, s): residue_counts(cache[g], len(labels[g]) - 1, s) for g in genes for s in STARS}

    confs = list(itertools.product(WINDOWS, COUNTS, BENIGN, STARS, ENRICH, SELF))
    gi = {g: i for i, g in enumerate(genes)}
    # S[c, g] = (tp, fp, fn): per-gene residue counts for each configuration.
    S = np.zeros((len(confs), len(genes), 3), dtype=np.int64)
    for ci, c in enumerate(confs):
        for g in genes:
            pr, lab = predict_residues(*counts[(g, c[3])], len(labels[g]) - 1, c), labels[g]
            S[ci, gi[g]] = ((pr & lab).sum(), (pr & ~lab).sum(), (~pr & lab).sum())
    cidx = {c: i for i, c in enumerate(confs)}

    def evaluate(gs, c):
        idx = [gi[g] for g in gs]
        tp, fp, fn = (int(v) for v in S[cidx[c], idx].sum(axis=0))
        p, r, f = prf(tp, fp, fn)
        jac = [S[cidx[c], gi[g]] for g in gs if labels[g].any()]
        jac = [a / (a + b + d) if a + b + d else 0.0 for a, b, d in jac]
        return dict(precision=round(p, 3), recall=round(r, 3), f1=round(f, 3),
                    macro_jaccard_defined=round(float(np.mean(jac)), 3) if jac else 0.0, tp=tp, fp=fp, fn=fn)

    def f1_all(gs):
        idx = [gi[g] for g in gs]
        t = S[:, idx].sum(axis=1).astype(float)
        p = np.divide(t[:, 0], t[:, 0] + t[:, 1], out=np.zeros(len(confs)), where=(t[:, 0] + t[:, 1]) > 0)
        r = np.divide(t[:, 0], t[:, 0] + t[:, 2], out=np.zeros(len(confs)), where=(t[:, 0] + t[:, 2]) > 0)
        return np.divide(2 * p * r, p + r, out=np.zeros(len(confs)), where=(p + r) > 0)

    grid = []
    for c in confs:
        grid.append(dict(config=cfg_name(c), **evaluate(genes, c)))
    grid.sort(key=lambda r: -r["f1"])
    with open(f"{a.out}_grid.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(grid[0]), delimiter="\t")
        w.writeheader()
        w.writerows(grid)

    defined = [g for g in genes if kind[g] == "defined"]
    na = [g for g in genes if kind[g] == "not_applicable"]
    lines = [f"genes: {len(defined)} with VCEP PM1 regions/residues, {len(na)} PM1 not applicable; "
             f"residues labelled positive {sum(int(labels[g].sum()) for g in genes):,} of "
             f"{sum(len(labels[g]) - 1 for g in genes):,}",
             "", "## Residue-level overlap, all genes (top 10 configurations)"]
    lines += [f"  {r['config']}: P {r['precision']} R {r['recall']} F1 {r['f1']} Jaccard {r['macro_jaccard_defined']}"
              for r in grid[:10]]
    d = next(r for r in grid if r["config"] == cfg_name(DEFAULT))
    lines.append(f"  current default {d['config']}: P {d['precision']} R {d['recall']} F1 {d['f1']} "
                 f"Jaccard {d['macro_jaccard_defined']} (rank {grid.index(d) + 1} of {len(grid)})")

    sel, rows = Counter(), []
    for i in range(a.repeats):
        rng = random.Random(3000 + i)
        dd, nn = defined[:], na[:]
        rng.shuffle(dd)
        rng.shuffle(nn)
        tr = dd[: len(dd) // 2] + nn[: len(nn) // 2]
        te = dd[len(dd) // 2:] + nn[len(nn) // 2:]
        f = f1_all(tr)
        top = np.flatnonzero(np.isclose(f, f.max()))
        best = DEFAULT if cidx[DEFAULT] in top else confs[int(top[0])]
        sel[best] += 1
        m, dm = evaluate(te, best), evaluate(te, DEFAULT)
        rows.append(dict(rep=i, selected=cfg_name(best), test_f1=m["f1"], test_precision=m["precision"],
                         test_recall=m["recall"], default_f1=dm["f1"]))
    with open(f"{a.out}_splits.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    mean = lambda k: round(float(np.mean([r[k] for r in rows])), 3)  # noqa: E731
    lines += ["", f"## Gene-level split ({a.repeats} repeats): residue-level F1 on the test half",
              f"  selected configuration: mean F1 {mean('test_f1')} (precision {mean('test_precision')}, "
              f"recall {mean('test_recall')}); current default: mean F1 {mean('default_f1')}",
              "  most often selected:"]
    lines += [f"    {k:>4}  {cfg_name(c)}" for c, k in sel.most_common(5)]

    best = sel.most_common(1)[0][0]
    hot = PM1Hotspots(a.hotspots)
    lines += ["", "## Variant-level check on eRepo missense variants (reference: VCEP PM1), query variant excluded"]
    for label, c in (("current default", DEFAULT), ("selected", best)):
        v = variant_eval(con, a.results, a.erepo, hot, c, cache, lengths)
        for grp, m in v.items():
            lines.append(f"  {label:<15} {grp:<16} {m}")
    lines.append("  (heuristic_genes = genes without a VCEP PM1 table, not used for calibration; "
                 "all_pm1_genes overlaps the calibration genes)")

    txt = "\n".join(lines) + "\n"
    Path(f"{a.out}_summary.txt").write_text(txt, encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
