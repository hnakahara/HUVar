"""PS3 / BS3 from calibrated functional assays.

Two sources, chosen per gene by ``resources/shared/functional_evidence_rules.tsv``:

* ``enigma_table9`` (BRCA1/BRCA2): the ClinGen ENIGMA VCEP Specifications Table 9,
  which gives, per variant, the PS3 or BS3 code the VCEP assigns from calibrated
  assays (splicing already handled by the VCEP). A variant coded ``None`` or not in
  the table gets no functional code.
* ``oddspath`` (default for genes without a VCEP specification): MaveDB score sets
  with a *clinical* calibration (not research-use-only; ExCALIBR is excluded). The
  variant's score is placed in the calibration's functional ranges; the strength is
  re-derived from the range's OddsPath with the Brnich et al. 2019 thresholds
  (Genome Med 12:3, Table 3) and capped per gene. Calibration choice per score set:
  the MaveDB *primary* calibration (``Config.mavedb_calibration = "primary"``, the
  one the score-set team designates), or the most conservative clinical one.
  Several score sets for one variant are not stacked: the strongest call is used,
  and opposite calls (PS3 vs BS3) of equal strength cancel (Brnich 2019).

Genes whose VCEP restricts PS3/BS3 to an approved assay list or decision tree that
does not include multiplexed assays (``not_applicable``; by default every gene with
a VCEP specification that is not listed in the rules table) get nothing from MaveDB.

Research-use-only calibrations (opt-in, ``Config.mavedb_research_use`` /
``--with-mavedb-research``; e.g. ExCALIBR). MaveDB flags these as not reviewed for
clinical use, so they are never read by default. When enabled:

* a score set's clinical calibrations always take precedence; its RUO
  calibrations are used only when it has no clinical one, and then
  conservatively (all RUO calibrations must agree on PS3 vs BS3; the weakest
  strength is used; any disagreement gives no call);
* for ENIGMA-governed BRCA1/BRCA2, Table 9 stays authoritative — a variant listed
  there (including "None") keeps the VCEP code — and only variants absent from
  Table 9 are evaluated with MaveDB (clinical calibrations first, then RUO);
* ``not_applicable`` genes are unaffected (their VCEP excludes multiplexed assays);
* the evidence text marks every RUO-based call "research use only".

Safeguards: null consequences (PVS1 territory) are skipped; protein-level (cDNA)
assays do not measure splicing, so BS3 from them is withheld when the splice
predictor gives >= 0.2 (ENIGMA's SpliceAI cutoff) or no prediction.
"""
from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional

from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength

SRC_MAVEDB = "[MaveDB]"
SRC_ENIGMA_T9 = "[ENIGMA Table 9]"

S = CriterionStrength
_RANK = {S.SUPPORTING: 1, S.MODERATE: 2, S.THREE_POINT: 3, S.STRONG: 4, S.VERY_STRONG: 8}
_FROM_TEXT = {
    "supporting": S.SUPPORTING, "moderate": S.MODERATE, "moderate_plus": S.MODERATE,
    "strong": S.STRONG, "very_strong": S.VERY_STRONG, "verystrong": S.VERY_STRONG,
}
_NULL = {
    ConsequenceType.FRAMESHIFT, ConsequenceType.STOP_GAINED, ConsequenceType.SPLICE_ACCEPTOR,
    ConsequenceType.SPLICE_DONOR, ConsequenceType.START_LOST, ConsequenceType.TRANSCRIPT_ABLATION,
}
_PROTEIN_ELIGIBLE = {
    ConsequenceType.MISSENSE, ConsequenceType.SYNONYMOUS,
    ConsequenceType.INFRAME_INSERTION, ConsequenceType.INFRAME_DELETION,
}
SPLICE_CUTOFF = 0.2
_DEFAULT_MAX = (S.STRONG, S.STRONG)
_REPO = Path(__file__).resolve().parents[3] / "resources" / "shared"
_TABLE9_TX = {"BRCA1": "NM_007294", "BRCA2": "NM_000059"}


def oddspath_strength(criterion: str, odds: float) -> Optional[CriterionStrength]:
    """Brnich et al. 2019 Table 3: OddsPath -> evidence strength (None = indeterminate)."""
    if criterion == "PS3":
        for cut, st in ((350, S.VERY_STRONG), (18.7, S.STRONG), (4.3, S.MODERATE), (2.1, S.SUPPORTING)):
            if odds > cut:
                return st
        return None
    if criterion == "BS3":
        for cut, st in ((0.053, S.STRONG), (0.23, S.MODERATE), (0.48, S.SUPPORTING)):
            if odds < cut:
                return st
    return None


def _cap(st: Optional[CriterionStrength], cap: Optional[CriterionStrength]) -> Optional[CriterionStrength]:
    if st is None or cap is None:
        return None
    return st if _RANK[st] <= _RANK[cap] else cap


@dataclass
class FunctionalCall:
    """Outcome of the functional-assay lookup for one variant."""
    criterion: Optional[ACMGCriterion]       # PS3 / BS3, or None (assessed, no code)
    strength: Optional[CriterionStrength]
    evidence: str


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _rows(path: Path):
    with _open(path) as fh:
        yield from csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")


def _resolve(path, name: str) -> Path:
    """Configured path when it exists, else the packaged resources copy."""
    if isinstance(path, (str, Path)) and Path(path).exists():
        return Path(path)
    return _REPO / name


@lru_cache(maxsize=4)
def load_rules(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    return {r["gene"]: r for r in _rows(path)}


def _is_ruo(r: dict) -> bool:
    return (r.get("ruo") or "0") == "1"


@lru_cache(maxsize=4)
def load_mavedb(path: Path, research_use: bool = False) -> dict[str, list[dict]]:
    """key (chrN:pos:ref:alt or NP_x:p.X) -> rows (one per score set x calibration).

    Research-use-only rows (``ruo=1``) are kept only when ``research_use`` is set,
    so the default run neither reads nor holds them in memory. Tables built before
    the ``ruo`` column existed contain clinical rows only."""
    out: dict[str, list[dict]] = {}
    if not path.exists():
        return out
    from sys import intern
    for r in _rows(path):
        if _is_ruo(r) and not research_use:
            continue
        # Most fields repeat across hundreds of thousands of rows (gene, score-set and
        # calibration URNs, labels): interning them roughly halves the table's memory.
        for f in ("gene", "urn", "level", "cal_urn", "primary", "criterion",
                  "mavedb_strength", "label", "ruo", "cal_title"):
            v = r.get(f)
            if v is not None:
                r[f] = intern(v)
        out.setdefault(r["key"], []).append(r)
    return out


@lru_cache(maxsize=4)
def load_scoresets(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    return {r["urn"]: r for r in _rows(path)}


@lru_cache(maxsize=4)
def load_table9(path: Path) -> dict[tuple[str, str], dict]:
    if not path.exists():
        return {}
    return {(r["gene"], r["hgvs_c"]): r for r in _rows(path)}


@lru_cache(maxsize=8)
def _has_ruo_column(path: Path) -> bool:
    """Whether a MaveDB table carries the ``ruo`` column (built since RUO support)."""
    if not path.exists():
        return False
    with _open(path) as fh:
        for ln in fh:
            if not ln.startswith("#"):
                return "ruo" in ln.rstrip("\n").split("\t")
    return False


def _strength_of(name: str) -> Optional[CriterionStrength]:
    return _FROM_TEXT.get((name or "").strip().lower().replace(" ", "_"))


class FunctionalAssays:
    """Per-variant PS3/BS3 call from ENIGMA Table 9 / MaveDB (see module docstring)."""

    def __init__(self, cfg) -> None:
        self._cfg = cfg
        self._enabled = getattr(cfg, "functional_assays", True) is True
        self._mode = getattr(cfg, "mavedb_calibration", "primary")
        self._ruo = getattr(cfg, "mavedb_research_use", False) is True
        self._rules = load_rules(_resolve(getattr(cfg, "functional_rules_tsv", None),
                                          "functional_evidence_rules.tsv"))
        self._mavedb_path = _resolve(getattr(cfg, "mavedb_functional_tsv", None), "mavedb_functional.tsv.gz")
        self._sets_path = _resolve(getattr(cfg, "mavedb_scoresets_tsv", None), "mavedb_scoresets.tsv")
        if self._ruo and not _has_ruo_column(self._mavedb_path):
            # A data/shared copy built before RUO support holds clinical rows only,
            # which would make the option a silent no-op: use the packaged table
            # (same pinned MaveDB release) and say so.
            packaged = _REPO / "mavedb_functional.tsv.gz"
            if packaged != self._mavedb_path and _has_ruo_column(packaged):
                import structlog
                structlog.get_logger().warning(
                    "mavedb_table_without_ruo", path=str(self._mavedb_path),
                    using=str(packaged),
                    hint="rebuild with setup_data.py --only mavedb --force-mavedb")
                self._mavedb_path = packaged
                self._sets_path = _REPO / "mavedb_scoresets.tsv"
        self._t9_path = _resolve(getattr(cfg, "enigma_table9_tsv", None), "enigma_brca12_table9.tsv")
        from acmg_classifier.classification.gene_validity import _vcep_genes
        # VCEP genes decide the not_applicable default; fall back to the packaged
        # table so a missing data copy never opens MaveDB up to every VCEP gene.
        self._vcep = _vcep_genes(_resolve(getattr(cfg, "disease_prevalence_tsv", None),
                                          "disease_prevalence.tsv"))

    # -- policy -------------------------------------------------------------------
    def policy(self, gene: str) -> tuple[str, Optional[CriterionStrength], Optional[CriterionStrength], str]:
        r = self._rules.get(gene)
        if r is not None:
            return (r["policy"], _strength_of(r["ps3_max"]), _strength_of(r["bs3_max"]),
                    (r.get("spec") or "").strip())
        if gene in self._vcep:
            return "not_applicable", None, None, self._vcep[gene]
        return ("oddspath",) + _DEFAULT_MAX + ("no VCEP specification",)

    # -- entry point ----------------------------------------------------------------
    def call(self, variant, annotation) -> Optional[FunctionalCall]:
        """None = no functional-assay data for this variant (callers fall back)."""
        if not self._enabled:
            return None
        pc = annotation.primary_consequence
        if pc is None:
            return None
        policy, ps3_max, bs3_max, spec = self.policy(pc.gene_symbol)
        if policy == "enigma_table9":
            # VCEP-curated per variant, splicing included: applied as listed.
            t9 = self._table9(annotation, pc.gene_symbol)
            if t9 is not None or not self._ruo or pc.consequence in _NULL:
                return t9
            # Opt-in (research use): a variant absent from Table 9 is evaluated
            # with MaveDB, labelled as outside the ENIGMA specification.
            return self._mavedb(variant, annotation, pc, ps3_max, bs3_max,
                                f"{spec}; not in ENIGMA Table 9", beyond_vcep=True)
        if policy != "oddspath" or pc.consequence in _NULL:
            return None
        return self._mavedb(variant, annotation, pc, ps3_max, bs3_max, spec)

    # -- ENIGMA Table 9 ---------------------------------------------------------------
    def _table9(self, annotation, gene: str) -> Optional[FunctionalCall]:
        rc = annotation.rule_consequence or annotation.primary_consequence
        acc, _, c = (rc.hgvs_c or "").partition(":")
        if not c or acc.split(".")[0] != _TABLE9_TX.get(gene):
            return None
        r = load_table9(self._t9_path).get((gene, c))
        if r is None:
            return None
        why = f"{SRC_ENIGMA_T9} {gene} {c}: {r['text']}"
        if r["results"]:
            why += f" [{r['results']}]"
        code = r["code"].strip()
        if code in ("PS3", "BS3"):
            st = _strength_of(r["strength"]) or S.STRONG
            return FunctionalCall(ACMGCriterion(code), st, why)
        return FunctionalCall(None, None, why + " -> no PS3/BS3 (ENIGMA Table 9: None)")

    # -- MaveDB ---------------------------------------------------------------------
    def _keys(self, variant, annotation, pc) -> list[tuple[str, str]]:
        keys: list[tuple[str, str]] = []
        if getattr(variant.assembly, "value", variant.assembly) == "GRCh38":
            chrom = variant.chrom if variant.chrom.startswith("chr") else "chr" + variant.chrom
            keys.append(("genomic", f"{chrom}:{variant.pos}:{variant.ref}:{variant.alt}"))
        if pc.consequence in _PROTEIN_ELIGIBLE:
            seen = set()
            for c in annotation.consequences:
                if c.gene_symbol != pc.gene_symbol or not c.hgvs_p or ":" not in c.hgvs_p:
                    continue
                acc, _, p = c.hgvs_p.partition(":")
                k = f"{acc.split('.')[0]}:{p.replace('%3D', '=')}"
                if k not in seen:
                    seen.add(k)
                    keys.append(("protein", k))
        return keys

    def _per_set(self, rows: list[dict]) -> tuple[Optional[str], Optional[CriterionStrength], dict]:
        """Calibration choice for one score set -> (criterion, Brnich strength, row used).

        Clinical calibrations take precedence; research-use-only ones (present only
        when opted in) are used when the set has no clinical calibration, and then
        always conservatively (none of them is a MaveDB primary for clinical use)."""
        clinical = [r for r in rows if not _is_ruo(r)]
        if clinical:
            return self._choose(clinical)
        return self._conservative(rows)

    def _choose(self, rows: list[dict]) -> tuple[Optional[str], Optional[CriterionStrength], dict]:
        """Clinical calibrations: the MaveDB primary, or the most conservative."""
        primary = [r for r in rows if r["primary"] == "1"]
        if self._mode == "primary" and primary:
            c, st = self._derive(primary[0])
            return c, st, primary[0]
        return self._conservative(rows)

    @staticmethod
    def _derive(r: dict) -> tuple[Optional[str], Optional[CriterionStrength]]:
        crit = r["criterion"]
        if crit not in ("PS3", "BS3"):
            return None, None
        if r["oddspath"]:
            st = oddspath_strength(crit, float(r["oddspath"]))
        else:
            st = _strength_of(r["mavedb_strength"])
        return (crit, st) if st else (None, None)

    def _conservative(self, rows: list[dict]) -> tuple[Optional[str], Optional[CriterionStrength], dict]:
        """All calibrations must agree on the criterion; the weakest strength wins."""
        derive = self._derive
        derived = [(derive(r), r) for r in rows]
        crits = {c for (c, _), _ in derived}
        if len(crits) != 1 or None in crits:
            return None, None, rows[0]
        (c, _), r = min(derived, key=lambda x: _RANK[x[0][1]])
        return c, min((st for (_, st), _ in derived), key=lambda s: _RANK[s]), r

    def _mavedb(self, variant, annotation, pc, ps3_max, bs3_max, spec,
                beyond_vcep: bool = False) -> Optional[FunctionalCall]:
        table = load_mavedb(self._mavedb_path, self._ruo)
        hits: dict[str, list[dict]] = {}
        level: dict[str, str] = {}
        for lvl, k in self._keys(variant, annotation, pc):
            for r in table.get(k, ()):
                if r["gene"] != pc.gene_symbol or r["level"] != lvl:
                    continue
                hits.setdefault(r["urn"], []).append(r)
                level[r["urn"]] = lvl
        if not hits:
            return None
        sets = load_scoresets(self._sets_path)
        sp = annotation.splice
        splice = sp.max_delta if (sp is not None and sp.is_available) else None
        calls, notes = [], []
        used_ruo = used_label = False
        for urn, rows in sorted(hits.items()):
            crit, st, r = self._per_set(rows)
            meta = sets.get(urn, {})
            ref = (f"{meta.get('first_author', '')} {meta.get('year', '')}".strip()
                   or meta.get("title", "")) + (f", PMID {meta['publication']}"
                                                if meta.get("publication", "").isdigit() else "")
            ruo = _is_ruo(r)
            used_ruo = used_ruo or ruo
            if crit is not None and not r["oddspath"]:
                used_label = True
            desc = (f"{urn} ({ref}) score={r['score']}"
                    + (f" in '{r['label']}'" if r["label"] else "")
                    + (f", OddsPath {r['oddspath']}" if r["oddspath"] else "")
                    + (f" [research use only: {r.get('cal_title') or 'RUO calibration'}]"
                       if ruo else ""))
            if crit is None:
                notes.append(f"{desc}: indeterminate")
                continue
            cap = ps3_max if crit == "PS3" else bs3_max
            capped = _cap(st, cap)
            if capped is None:
                notes.append(f"{desc}: {crit} not applicable for {pc.gene_symbol} ({spec})")
                continue
            if crit == "BS3" and level[urn] == "protein" and (splice is None or splice >= SPLICE_CUTOFF):
                why = "no splice prediction" if splice is None else f"splice score {splice:.2f} >= {SPLICE_CUTOFF}"
                notes.append(f"{desc}: BS3 withheld (protein-level assay; {why})")
                continue
            tag = "" if capped == st else f" (capped from {st.value}: {spec})"
            calls.append((crit, capped, f"{desc} -> {crit}_{capped.value}{tag}"))
        basis = ("calibration: MaveDB primary" if self._mode == "primary"
                 else "calibration: most conservative") + "; strength per Brnich et al. 2019 OddsPath"
        if used_label:
            basis += "; strength from the MaveDB calibration label where no OddsPath is given"
        if used_ruo:
            basis += ("; includes research-use-only calibrations (not reviewed for clinical "
                      "use; most conservative)")
        if beyond_vcep:
            basis += f"; outside the VCEP specification ({spec})"
        if not calls:
            return FunctionalCall(None, None, f"{SRC_MAVEDB} " + "; ".join(notes) + f" [{basis}]")
        best = {c: max((x for x in calls if x[0] == c), key=lambda x: _RANK[x[1]], default=None)
                for c in ("PS3", "BS3")}
        ps3, bs3 = best["PS3"], best["BS3"]
        listing = "; ".join(x[2] for x in calls) + ("; " + "; ".join(notes) if notes else "")
        if ps3 and bs3:
            if _RANK[ps3[1]] == _RANK[bs3[1]]:
                return FunctionalCall(None, None, f"{SRC_MAVEDB} conflicting assays of equal strength "
                                      f"(not used, Brnich 2019): {listing} [{basis}]")
            win = ps3 if _RANK[ps3[1]] > _RANK[bs3[1]] else bs3
        else:
            win = ps3 or bs3
        return FunctionalCall(ACMGCriterion(win[0]), win[1], f"{SRC_MAVEDB} {listing} [{basis}]")
