"""DuckDB query layer for gnomAD exome data (BA1/BS1/BS2/PM2)."""
from __future__ import annotations
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.models.annotation import GnomADData
from acmg_classifier.utils.chrom import chrom_candidates

log = structlog.get_logger()


def _pass_filter(filters) -> bool:
    """gnomAD FILTER conventions: None / "" / "PASS" / "." mean passed QC."""
    return filters is None or str(filters).strip().upper() in ("", "PASS", ".")


def _merge_rows(rows: list) -> GnomADData:
    """Merge multiple gnomAD rows for one variant by per-field MAX.

    Row layout: (af, an, ac, nhomalt, nhemi, popmax_af, popmax_pop,
    faf95_popmax, af_xy, filters). Only PASS rows are merged when any exist
    (a filtered record must not contribute a frequency); if every row is
    filtered the variant is reported as filter-failed."""
    pass_rows = [r for r in rows if _pass_filter(r[9])]
    use = pass_rows or rows

    def fmax(idx: int):
        vals = [r[idx] for r in use if r[idx] is not None]
        return max(vals) if vals else None

    # popmax_pop comes from whichever used row has the highest popmax AF.
    best = max(use, key=lambda r: (r[5] if r[5] is not None else -1.0))
    return GnomADData(
        af=fmax(0),
        an=fmax(1),
        ac=fmax(2),
        nhomalt=fmax(3),
        nhemi=fmax(4),
        popmax_af=fmax(5),
        popmax_pop=best[6],
        faf95_popmax=fmax(7),
        af_xy=fmax(8),
        filter_pass=bool(pass_rows),
    )


class GnomADDB:
    def __init__(self, db_path: Path, constraint_tsv: Path) -> None:
        self._db_path = db_path
        self._constraint_tsv = constraint_tsv
        self._constraint: dict[str, tuple[float | None, float | None, float | None]] = {}
        if constraint_tsv.exists():
            self._constraint = _load_constraint(constraint_tsv)
        # Whether the variants table carries the af_xy column (added for X-linked
        # "in males" BA1/BS1). Probed once and cached; a DB built before this
        # column lacks it, so we degrade gracefully to NULL rather than erroring.
        self._has_af_xy: bool | None = None

    def query(self, chrom: str, pos: int, ref: str, alt: str) -> Optional[GnomADData]:
        """Fetch population statistics for a specific variant.

        Returns a synthetic "absent" record (AF=0, AC=0, filter_pass=True)
        when the variant is missing from the database. This is intentional —
        downstream criteria distinguish "absent" (rare, supports PM2) from
        "filter-failed" (untrustworthy, supports neither side), so we must
        differentiate them at this layer."""
        if not self._db_path.exists():
            log.warning("gnomad_db_missing", path=str(self._db_path))
            return None
        # chrom_candidates handles "1" vs "chr1" — gnomAD raw files have
        # historically used both depending on version/source.
        c1, c2 = chrom_candidates(chrom)
        try:
            import duckdb  # lazy: keeps the merge helpers importable without duckdb
            con = duckdb.connect(str(self._db_path), read_only=True)
            # Select af_xy only when the schema has it; otherwise NULL keeps the
            # result tuple shape constant for older DBs.
            xy_expr = "af_xy" if self._af_xy_available(con) else "NULL AS af_xy"
            rows = con.execute(
                f"""
                SELECT af, an, ac, nhomalt, nhemi,
                       popmax_af, popmax_pop, faf95_popmax, {xy_expr},
                       filters
                FROM variants
                WHERE chrom IN (?, ?) AND pos = ? AND ref = ? AND alt = ?
                """,
                [c1, c2, pos, ref, alt],
            ).fetchall()
            con.close()
        except Exception as exc:
            log.error("gnomad_query_error", error=str(exc))
            return None

        # "Absent from gnomAD" is the most-informative PM2 signal — emit a
        # well-formed record rather than None so callers don't need a special
        # branch for the missing-record case.
        if not rows:
            return GnomADData(af=0.0, ac=0, an=0, filter_pass=True)

        # A variant may have several rows: the GRCh37 build loads gnomAD exomes
        # AND genomes (no joint release at v2.1.1), so a variant present in both
        # appears twice. Merge by per-field MAX over the PASS records — this is
        # the "either dataset meets the criterion" semantics (the higher AF wins
        # for BA1/BS1, and PM2's rarity check sees the higher AF too). GRCh38's
        # joint build has a single row, so the merge is a no-op there.
        return _merge_rows(rows)

    def _af_xy_available(self, con) -> bool:
        """True if the variants table has the af_xy column (cached per instance)."""
        if self._has_af_xy is None:
            cols = {r[1] for r in con.execute("PRAGMA table_info('variants')").fetchall()}
            self._has_af_xy = "af_xy" in cols
        return self._has_af_xy

    def get_constraint(
        self, gene_symbol: str
    ) -> tuple[float | None, float | None, float | None]:
        """Return (pLI, LOEUF, missense Z-score) for the gene, or all-None if absent."""
        return self._constraint.get(gene_symbol, (None, None, None))

    def enrich_with_constraint(self, data: GnomADData, gene_symbol: str) -> GnomADData:
        pli, loeuf, mis_z = self.get_constraint(gene_symbol)
        return data.model_copy(update={"pli": pli, "loeuf": loeuf, "mis_z": mis_z})


def _load_constraint(
    tsv: Path,
) -> dict[str, tuple[float | None, float | None, float | None]]:
    """Load per-gene (pLI, LOEUF, missense Z-score) from a gnomAD constraint table.

    The table has one row per *transcript*, so each gene appears many times. A
    single representative row is chosen per gene with this priority:
      1. value-bearing rows (a real LOEUF) before NA rows;
      2. MANE Select, else canonical, else any transcript;
      3. ties broken by the smallest (most constrained) LOEUF.

    Missense Z-score (mis.z_score in v4.1 / mis_z in v2.1.1) is read from the
    same chosen row and fed to PP2 as an alternative qualifier alongside the
    ClinVar benign-missense rate.
    """
    import csv
    import math

    # Column names differ between gnomAD v4.1 (newer dotted form) and
    # v2.1.1 (legacy short form). Candidate lists let one builder work
    # against either constraint release without manual translation.
    _PLI_COLS = ["pLI", "lof.pLI"]
    _LOEUF_COLS = ["oe_lof_upper", "lof.oe_ci.upper"]
    _MIS_Z_COLS = ["mis_z", "mis.z_score"]

    def _get_float(row: dict, cols: list[str]) -> float | None:
        for c in cols:
            v = row.get(c, "")
            if v and v not in ("NA", ".", "nan"):
                try:
                    return float(v)
                except ValueError:
                    pass
        return None

    def _is_true(row: dict, col: str) -> bool:
        return str(row.get(col, "")).strip().lower() == "true"

    result: dict[str, tuple[float | None, float | None, float | None]] = {}
    best_key: dict[str, tuple[int, int, float]] = {}
    with tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = row.get("gene", row.get("gene_id", "")).strip()
            if not gene:
                continue
            loeuf = _get_float(row, _LOEUF_COLS)
            rank = 0 if _is_true(row, "mane_select") else (1 if _is_true(row, "canonical") else 2)
            # Lower key wins: value-bearing first, then MANE/canonical, then min LOEUF.
            key = (0 if loeuf is not None else 1, rank, loeuf if loeuf is not None else math.inf)
            if gene not in best_key or key < best_key[gene]:
                best_key[gene] = key
                result[gene] = (
                    _get_float(row, _PLI_COLS),
                    loeuf,
                    _get_float(row, _MIS_Z_COLS),
                )
    return result
