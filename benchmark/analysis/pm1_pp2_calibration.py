#!/usr/bin/env python3
"""Reviewer 3, comments 7 (PM1) and 9 (PP2): calibration of HUVar's ClinVar heuristics.

Run on the analysis server (needs the ClinVar SQLite and the gnomAD constraint
table under $DATA, and the HUVar package):

  python pm1_pp2_calibration.py --data-dir $DATA \
      --results erepo_benchmark/inputs/erepo_hg38_huhvar_default_results.tsv \
      --erepo <HUVar>/resources/GRCh38/erepo_classifications_hg38.tsv \
      --out results/calib

PM1 (comment 7)
  The fallback heuristic (>= N other P/LP missense changes within +/- W amino
  acids, no B/LB missense in the window, >= 1 star, query variant excluded) is
  evaluated on eRepo missense variants for a grid W in {10, 25, 50} x
  N in {2, 3, 4, 5}; reference = PM1 met by the VCEP (any strength).
    group "heuristic_genes": genes where HUVar actually uses the heuristic
                             (no VCEP hotspot table, PM1 not declared N/A);
    group "all_pm1_genes":   all eRepo genes whose VCEP allows PM1 (heuristic
                             applied as if no VCEP table existed).

PP2 (comment 9)
  Gene-level reference = VCEP PP2 decision (applicable / not applicable) in
  disease_prevalence.tsv. Per gene: P/LP missense, ClinVar B/LB missense,
  gnomAD common missense (pp2_gene_stats.tsv), missense Z. Grid:
  min P/LP in {5, 10, 15, 20} x max benign fraction in {0.05, 0.10, 0.15,
  0.20, 0.30} x Z branch {off, Z>=3.09 with benign <= 0.15, Z>=3.09 with
  benign <= 0.30}, with and without the gnomAD common-missense correction.
  Thresholds are chosen on a stratified random half of the genes (balanced
  accuracy; ties -> the current default, then the simpler rule) and evaluated
  on the other half (seed 2026); 500 repeated splits summarise stability.

Outputs <out>_pm1_grid.tsv, <out>_pp2_genes.tsv, <out>_pp2_grid.tsv,
<out>_pp2_splits.tsv and <out>_summary.txt.
"""
from __future__ import annotations

import argparse
import csv
import random
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from acmg_classifier.config import Config
from acmg_classifier.criteria.pm1_hotspots import PM1Hotspots
from acmg_classifier.criteria.pp2_genes import PP2Applicability
from acmg_classifier.local_db.clinvar_sqlite import _PP2_BENIGN, _PP2_PATH, _is_missense_p
from acmg_classifier.local_db.gnomad_db import _load_constraint

_MIS = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?$")
WINDOWS = (10, 25, 50)
COUNTS = (2, 3, 4, 5)
PP2_MIN_PATH = (5, 10, 15, 20)
PP2_MAX_FRAC = (0.05, 0.10, 0.15, 0.20, 0.30)
PP2_Z = (None, (3.09, 0.15), (3.09, 0.30))
DEFAULT_PP2 = (10, 0.05, (3.09, 0.15), True)


def strip_chr(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def connect(path: Path) -> sqlite3.Connection:
    return sqlite3.connect("file:" + path.as_posix() + "?mode=ro&immutable=1", uri=True)


def metrics(tp, fp, fn, tn=None):
    prec = tp / (tp + fp) if tp + fp else float("nan")
    sens = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * prec * sens / (prec + sens) if (prec + sens) and prec == prec and sens == sens else float("nan")
    out = dict(tp=tp, fp=fp, fn=fn, precision=round(prec, 3), sensitivity=round(sens, 3), f1=round(f1, 3))
    if tn is not None:
        spec = tn / (tn + fp) if tn + fp else float("nan")
        out.update(tn=tn, specificity=round(spec, 3), balanced_acc=round((sens + spec) / 2, 3))
    return out


# ----------------------------------------------------------------------------- PM1
def run_pm1(cfg: Config, results: Path, erepo: Path, con, out: str, lines: list[str]) -> None:
    truth = {}
    with erepo.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            codes = {c.strip().split("_")[0] for c in (r.get("met_codes") or "").split(",")}
            truth[r["variant_id"]] = "PM1" in codes
    hot = PM1Hotspots(cfg.pm1_hotspots_tsv)
    variants = []  # (gene, aa_pos, chrom, pos, ref, alt, truth, group_heuristic)
    with results.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = _MIS.match(r.get("hgvs_p") or "")
            if not m or m.group(3) in ("Ter",) or m.group(1) == m.group(3):
                continue
            vid, gene = r["variant_id"], r.get("gene") or ""
            if vid not in truth or not gene or hot.is_not_applicable(gene):
                continue
            heur = not (hot.has_gene(gene) or hot.has_cys_creating(gene))
            chrom, pos, ref, alt = vid.split(":")
            variants.append((gene, int(m.group(2)), strip_chr(chrom), int(pos), ref, alt, truth[vid], heur))

    by_gene = defaultdict(list)
    for v in variants:
        by_gene[v[0]].append(v)
    q = ("SELECT chrom, pos, ref, alt, hgvs_p, amino_acid_change, clinical_significance, codon_position "
         "FROM variants WHERE gene_symbol = ? AND star_rating >= 1 AND clinical_significance IN (?,?,?,?,?,?)")
    counts = {g: Counter() for g in ("heuristic_genes", "all_pm1_genes")}
    for gene, vs in by_gene.items():
        rows = [r for r in con.execute(q, (gene, *_PP2_PATH, *_PP2_BENIGN)).fetchall()
                if _is_missense_p(r[4]) and r[7] is not None]
        for (_, aa, chrom, pos, ref, alt, t, heur) in vs:
            for w in WINDOWS:
                path, benign = set(), False
                for rc, rp, rr, ra, hp, aac, sig, cod in rows:
                    if abs(int(cod) - aa) > w:
                        continue
                    if strip_chr(str(rc)) == chrom and rp == pos and rr == ref and ra == alt:
                        continue
                    if sig in _PP2_PATH:
                        path.add(aac or hp)
                    else:
                        benign = True
                for n in COUNTS:
                    pred = (not benign) and len(path) >= n
                    key = (w, n, "tp" if pred and t else "fp" if pred else "fn" if t else "tn")
                    counts["all_pm1_genes"][key] += 1
                    if heur:
                        counts["heuristic_genes"][key] += 1

    rows_out = []
    for grp, c in counts.items():
        nvar = sum(1 for v in variants if grp == "all_pm1_genes" or v[7])
        ngene = len({v[0] for v in variants if grp == "all_pm1_genes" or v[7]})
        npos = sum(1 for v in variants if (grp == "all_pm1_genes" or v[7]) and v[6])
        for w in WINDOWS:
            for n in COUNTS:
                m = metrics(c[(w, n, "tp")], c[(w, n, "fp")], c[(w, n, "fn")], c[(w, n, "tn")])
                rows_out.append(dict(group=grp, window=w, min_plp=n, variants=nvar, genes=ngene,
                                     erepo_pm1=npos, predicted=m["tp"] + m["fp"], **m))
    with open(f"{out}_pm1_grid.tsv", "w", newline="", encoding="utf-8") as fh:
        w_ = csv.DictWriter(fh, fieldnames=list(rows_out[0]), delimiter="\t")
        w_.writeheader()
        w_.writerows(rows_out)
    lines.append("## PM1 heuristic grid (reference: VCEP PM1 in the eRepo, any strength)")
    for grp in counts:
        sub = [r for r in rows_out if r["group"] == grp]
        lines.append(f"{grp}: {sub[0]['variants']} missense variants, {sub[0]['genes']} genes, "
                     f"eRepo PM1 {sub[0]['erepo_pm1']}")
        lines.append("  window min_plp predicted precision sensitivity F1")
        for r in sub:
            flag = "  <- default" if (r["window"], r["min_plp"]) == (25, 3) else ""
            lines.append(f"  {r['window']:>6} {r['min_plp']:>7} {r['predicted']:>9} {r['precision']:>9} "
                         f"{r['sensitivity']:>11} {r['f1']:>5}{flag}")
        best = max(sub, key=lambda r: (r["f1"] if r["f1"] == r["f1"] else -1))
        lines.append(f"  best F1: window {best['window']}, min {best['min_plp']} (F1 {best['f1']})")
    lines.append("")


# ----------------------------------------------------------------------------- PP2
def pp2_predict(g, min_path, max_frac, z, use_common):
    benign = g["benign_clinvar"] + (g["common"] if use_common else 0)
    total = g["path"] + benign
    if g["path"] < min_path:
        return False
    frac = benign / total if total else 0.0
    if frac <= max_frac:
        return True
    if z is not None and g["mis_z"] is not None:
        return g["mis_z"] >= z[0] and frac <= z[1]
    return False


def pp2_eval(genes, rule):
    tp = fp = fn = tn = 0
    for g in genes:
        p = pp2_predict(g, *rule)
        if p and g["truth"]:
            tp += 1
        elif p:
            fp += 1
        elif g["truth"]:
            fn += 1
        else:
            tn += 1
    return metrics(tp, fp, fn, tn)


def rule_name(rule):
    mp, mf, z, uc = rule
    zs = "off" if z is None else f"Z>={z[0]}&benign<={z[1]}"
    return f"minPLP={mp};maxBenign={mf};Zbranch={zs};gnomADcommon={'yes' if uc else 'no'}"


def choose(train, rules):
    def key(rule):
        m = pp2_eval(train, rule)
        ba = m["balanced_acc"] if m["balanced_acc"] == m["balanced_acc"] else -1
        return (ba, rule == DEFAULT_PP2, rule[2] is None)
    return max(rules, key=key)


def split(genes, rng):
    pos = [g for g in genes if g["truth"]]
    neg = [g for g in genes if not g["truth"]]
    rng.shuffle(pos)
    rng.shuffle(neg)
    tr = pos[: len(pos) // 2] + neg[: len(neg) // 2]
    te = pos[len(pos) // 2:] + neg[len(neg) // 2:]
    return tr, te


def run_pp2(cfg: Config, con, out: str, lines: list[str]) -> None:
    app = PP2Applicability(cfg.disease_prevalence_tsv)
    decided = {g: v for g, v in app._by_gene.items() if v in ("applicable", "not_applicable")}
    common = {}
    with Path(cfg.pp2_gene_stats_tsv).open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            common[r["gene"]] = int(r["common_missense"] or 0)
    constraint = _load_constraint(cfg.gnomad_constraint_tsv)
    q = ("SELECT hgvs_p, clinical_significance FROM variants WHERE gene_symbol = ? AND star_rating >= 1 "
         "AND clinical_significance IN (?,?,?,?,?,?)")
    genes = []
    for gene, dec in sorted(decided.items()):
        rows = con.execute(q, (gene, *_PP2_PATH, *_PP2_BENIGN)).fetchall()
        genes.append(dict(
            gene=gene, truth=dec == "applicable",
            path=sum(1 for hp, s in rows if s in _PP2_PATH and _is_missense_p(hp)),
            benign_clinvar=sum(1 for hp, s in rows if s in _PP2_BENIGN and _is_missense_p(hp)),
            common=common.get(gene, 0), mis_z=constraint.get(gene, (None, None, None))[2]))
    with open(f"{out}_pp2_genes.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(genes[0]) + ["default_rule"], delimiter="\t")
        w.writeheader()
        for g in genes:
            w.writerow({**g, "default_rule": pp2_predict(g, *DEFAULT_PP2)})

    rules = [(mp, mf, z, uc) for uc in (True, False) for mp in PP2_MIN_PATH for mf in PP2_MAX_FRAC for z in PP2_Z]
    grid = [dict(rule=rule_name(r), **pp2_eval(genes, r)) for r in rules]
    with open(f"{out}_pp2_grid.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(grid[0]), delimiter="\t")
        w.writeheader()
        w.writerows(grid)

    npos = sum(g["truth"] for g in genes)
    lines.append(f"## PP2 (reference: VCEP PP2 decision; {len(genes)} genes, applicable {npos}, "
                 f"not applicable {len(genes) - npos})")
    for label, r in (("current default (with gnomAD common)", DEFAULT_PP2),
                     ("current default without gnomAD common", DEFAULT_PP2[:3] + (False,))):
        lines.append(f"{label}, all genes: {pp2_eval(genes, r)}")

    rng = random.Random(2026)
    tr, te = split(genes, rng)
    best = choose(tr, rules)
    lines.append(f"primary split (seed 2026): train {len(tr)} genes ({sum(g['truth'] for g in tr)} applicable), "
                 f"test {len(te)} ({sum(g['truth'] for g in te)} applicable)")
    lines.append(f"  selected on train: {rule_name(best)}  train {pp2_eval(tr, best)}")
    lines.append(f"  test: {pp2_eval(te, best)}")
    lines.append(f"  current default on test: {pp2_eval(te, DEFAULT_PP2)}")

    sel, rows = Counter(), []
    for i in range(500):
        tr, te = split(genes, random.Random(10_000 + i))
        b = choose(tr, rules)
        sel[rule_name(b)] += 1
        m, d = pp2_eval(te, b), pp2_eval(te, DEFAULT_PP2)
        rows.append(dict(rep=i, selected=rule_name(b), test_sens=m["sensitivity"], test_spec=m["specificity"],
                         test_ba=m["balanced_acc"], default_sens=d["sensitivity"],
                         default_spec=d["specificity"], default_ba=d["balanced_acc"]))
    with open(f"{out}_pp2_splits.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    def mean(k):
        v = [r[k] for r in rows if r[k] == r[k]]
        return round(sum(v) / len(v), 3) if v else float("nan")
    lines.append("500 repeated stratified splits: mean test sensitivity / specificity / balanced accuracy")
    lines.append(f"  selected rule: {mean('test_sens')} / {mean('test_spec')} / {mean('test_ba')}")
    lines.append(f"  current default: {mean('default_sens')} / {mean('default_spec')} / {mean('default_ba')}")
    lines.append("  most often selected rules:")
    for r, c in sel.most_common(5):
        lines.append(f"    {c:>4}  {r}")
    lines.append("")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--erepo", type=Path, required=True)
    ap.add_argument("--out", default="calib")
    ap.add_argument("--skip-pm1", action="store_true")
    ap.add_argument("--skip-pp2", action="store_true")
    a = ap.parse_args()
    cfg = Config(data_dir=a.data_dir, assembly="GRCh38")
    con = connect(cfg.clinvar_sqlite)
    lines: list[str] = []
    if not a.skip_pm1:
        run_pm1(cfg, a.results, a.erepo, con, a.out, lines)
    if not a.skip_pp2:
        run_pp2(cfg, con, a.out, lines)
    txt = "\n".join(lines) + "\n"
    Path(f"{a.out}_summary.txt").write_text(txt, encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
