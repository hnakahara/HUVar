"""DuckDB query layer for gnomAD exome data (BA1/BS1/BS2/PM2)."""
from __future__ import annotations
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import structlog

from acmg_classifier.models.annotation import GnomADData
from acmg_classifier.utils.chrom import chrom_candidates

if TYPE_CHECKING:
    from acmg_classifier.models.variant import VariantRecord

log = structlog.get_logger()


def _select_expr(cols: set[str]) -> str:
    """Build the gnomAD `variants` SELECT column list, degrading gracefully
    against DBs built before af_xy / ac_xx / grpmax / af_non_cancer existed.

    Shared verbatim by the single-variant query() and the batch precompute()
    so both produce the identical 15-field row layout consumed by _merge_rows
    (indices: 0 af … 11 filters, 12/13 grpmax, 14 af_non_cancer)."""
    xy_expr = "af_xy" if "af_xy" in cols else "NULL AS af_xy"
    has_xx = "ac_xx" in cols and "nhomalt_xx" in cols
    xx_expr = "ac_xx, nhomalt_xx" if has_xx else "NULL AS ac_xx, NULL AS nhomalt_xx"
    has_grpmax = "ac_grpmax" in cols and "an_grpmax" in cols
    grpmax_expr = (
        "ac_grpmax, an_grpmax" if has_grpmax
        else "NULL AS ac_grpmax, NULL AS an_grpmax"
    )
    nc_expr = "af_non_cancer" if "af_non_cancer" in cols else "NULL AS af_non_cancer"
    return (
        "af, an, ac, nhomalt, nhemi, popmax_af, popmax_pop, faf95_popmax, "
        f"{xy_expr}, {xx_expr}, filters, {grpmax_expr}, {nc_expr}"
    )


def _pass_filter(filters) -> bool:
    """gnomAD FILTER conventions: None / "" / "PASS" / "." mean passed QC."""
    return filters is None or str(filters).strip().upper() in ("", "PASS", ".")


def _merge_rows(rows: list) -> GnomADData:
    """Merge multiple gnomAD rows for one variant by per-field MAX.

    Row layout: (af, an, ac, nhomalt, nhemi, popmax_af, popmax_pop,
    faf95_popmax, af_xy, ac_xx, nhomalt_xx, filters[, ac_grpmax, an_grpmax]).
    Only PASS rows are merged when any exist (a filtered record must not
    contribute a frequency); if every row is filtered the variant is reported as
    filter-failed. The optional trailing ac_grpmax/an_grpmax are present only for
    DBs built after that schema addition (older DBs degrade to None)."""
    pass_rows = [r for r in rows if _pass_filter(r[11])]
    use = pass_rows or rows

    def fmax(idx: int):
        vals = [r[idx] for r in use if r[idx] is not None]
        return max(vals) if vals else None

    # popmax_pop (and the co-located grpmax AC/AN) come from whichever used row
    # has the highest popmax AF. The tiebreak is fully deterministic and
    # independent of row order: on equal popmax AF prefer the larger overall AN
    # (the dataset with more samples), then the population name as a final
    # total-order key. Without it, exomes+genomes rows sharing a popmax AF (the
    # GRCh37 case) would pick popmax_pop by list position — which differs between
    # the single-variant query() and the batch precompute() JOIN, and even
    # between runs of the same path (DuckDB row order is not guaranteed).
    best = max(
        use,
        key=lambda r: (
            r[5] if r[5] is not None else -1.0,
            r[1] if r[1] is not None else -1,
            r[6] or "",
        ),
    )
    # GrpMax AC and AN must come from the SAME row (the dataset where the
    # subpopulation frequency is highest) so the upper-CI is computed on a
    # consistent count/number pair — never a per-field max that could mix them.
    has_grpmax = len(best) > 13
    # Non-cancer subset AF trails at index 14 (present only for DBs built after
    # that schema addition; older DBs / shorter test tuples degrade to None).
    af_non_cancer = None
    if use and all(len(r) > 14 for r in use):
        nc_vals = [r[14] for r in use if r[14] is not None]
        af_non_cancer = max(nc_vals) if nc_vals else None
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
        ac_xx=fmax(9),
        nhomalt_xx=fmax(10),
        ac_grpmax=best[12] if has_grpmax else None,
        an_grpmax=best[13] if has_grpmax else None,
        af_non_cancer=af_non_cancer,
        filter_pass=bool(pass_rows),
    )


class GnomADDB:
    def __init__(
        self,
        db_path: Path,
        constraint_tsv: Path,
        noncancer_db_path: Optional[Path] = None,
    ) -> None:
        self._db_path = db_path
        self._constraint_tsv = constraint_tsv
        # Companion non-cancer-subset DB (GRCh38 gnomAD v3.1.2). The main v4.1
        # build has af_non_cancer = NULL (v4 dropped the subset), so PM2 for
        # ENIGMA BRCA1/2 consults this as a fallback. None when not configured
        # (GRCh37, whose v2.1.1 build carries the subset inline). See query().
        self._noncancer_db_path = noncancer_db_path
        self._constraint: dict[str, tuple[float | None, float | None, float | None]] = {}
        if constraint_tsv.exists():
            self._constraint = _load_constraint(constraint_tsv)
        # The variants table's columns, probed once and cached. A DB built before
        # af_xy (X-linked "in males" BA1/BS1) or ac_xx/nhomalt_xx (female-only
        # BS2, e.g. TP53) lacks them, so query() degrades gracefully to NULL
        # rather than erroring (a female-only BS2 gene then withholds BS2).
        self._cols: set[str] | None = None
        # Per-batch cache populated by precompute(): variant.key -> GnomADData.
        # cached() reads this instead of opening a connection per variant. Left
        # empty for the single-variant explain path, where cached() falls back
        # to query().
        self._cache: dict[str, GnomADData] = {}

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
            # Select af_xy / ac_xx / nhomalt_xx only when the schema has them;
            # otherwise NULL keeps the result tuple shape constant for older DBs.
            # Column list (with graceful NULLs for older schemas) is shared with
            # precompute() via _select_expr so both paths return the identical
            # 15-field layout that _merge_rows consumes.
            select_expr = _select_expr(self._columns(con))
            rows = con.execute(
                f"""
                SELECT {select_expr}
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
        data = _merge_rows(rows)

        # Backfill the non-cancer-subset AF from the companion v3.1.2 DB when the
        # main build doesn't carry it (GRCh38 v4.1 dropped the non-cancer subset,
        # so its af_non_cancer is always NULL). Only consulted for a variant
        # present in the main DB — one absent there is absent in every subset
        # too, so PM2's "absent" path already covers it and we skip the extra
        # lookup for the overwhelmingly common novel variant.
        if (self._noncancer_db_path is not None
                and self._noncancer_db_path.exists()):
            # The companion subset is consultable → record that fact so a None
            # af_non_cancer downstream reads as genuine absence in the non-cancer
            # subset (PM2's "present only in cancer cohorts" case), not as missing
            # data. Set regardless of hit/miss; a miss IS the absence signal.
            update = {"non_cancer_queried": True}
            if data.af_non_cancer is None or data.faf95_non_cancer is None:
                nc = self._query_noncancer(c1, c2, pos, ref, alt)
                if nc is not None:
                    af_nc, faf_nc = nc
                    if data.af_non_cancer is None and af_nc is not None:
                        update["af_non_cancer"] = af_nc
                    if data.faf95_non_cancer is None and faf_nc is not None:
                        update["faf95_non_cancer"] = faf_nc
            data = data.model_copy(update=update)
        return data

    def precompute(self, variants: "list[VariantRecord]") -> None:
        """Batch-fetch gnomAD stats for every variant in a single JOIN.

        Replaces query()'s connection-per-variant with one connection and one
        hash-join against a temp key table, caching the result by variant.key.
        The row-processing (_merge_rows + non-cancer backfill) is byte-for-byte
        identical to query(); only the access pattern changes. cached() then
        serves the batch without touching DuckDB again."""
        if not variants:
            return
        if not self._db_path.exists():
            # Leave the cache empty: cached() falls back to query(), which emits
            # the same gnomad_db_missing warning. This keeps the missing-DB
            # behaviour identical to the pre-batch path.
            log.warning("gnomad_db_missing", path=str(self._db_path))
            return
        try:
            import duckdb
            con = duckdb.connect(str(self._db_path), read_only=True)
            select_expr = _select_expr(self._columns(con))
            con.execute(
                "CREATE TEMP TABLE _keys "
                "(vkey TEXT, chrom TEXT, pos BIGINT, ref TEXT, alt TEXT)"
            )
            # One key row per (variant, chrom-spelling). Only one spelling exists
            # in the DB, so emitting both candidates can never double-count — it
            # just lets the equijoin match whichever convention the DB was built
            # with (the IN (?, ?) contract of query(), expressed as a join).
            key_rows = [
                (v.key, cand, v.pos, v.ref, v.alt)
                for v in variants
                for cand in chrom_candidates(v.chrom)
            ]
            con.executemany("INSERT INTO _keys VALUES (?, ?, ?, ?, ?)", key_rows)
            rows = con.execute(
                f"""
                SELECT k.vkey, {select_expr}
                FROM _keys k
                JOIN variants v
                  ON v.chrom = k.chrom AND v.pos = k.pos
                 AND v.ref = k.ref AND v.alt = k.alt
                """
            ).fetchall()
            con.close()
        except Exception as exc:
            log.error("gnomad_precompute_error", error=str(exc))
            return

        # Group matched DB rows by variant key; strip the leading vkey column so
        # each grouped row matches the exact tuple layout _merge_rows expects.
        by_key: dict[str, list] = {}
        for row in rows:
            by_key.setdefault(row[0], []).append(row[1:])

        companion_ok = (
            self._noncancer_db_path is not None
            and self._noncancer_db_path.exists()
        )
        cache: dict[str, GnomADData] = {}
        need_nc: dict[str, tuple[str, str, int, str, str]] = {}
        for v in variants:
            matched = by_key.get(v.key)
            if not matched:
                # "Absent from gnomAD": the same synthetic record query() emits.
                cache[v.key] = GnomADData(af=0.0, ac=0, an=0, filter_pass=True)
                continue
            data = _merge_rows(matched)
            if companion_ok:
                # Mark the subset consultable regardless of hit/miss — mirrors
                # query()'s unconditional non_cancer_queried=True for present rows.
                data = data.model_copy(update={"non_cancer_queried": True})
                if data.af_non_cancer is None or data.faf95_non_cancer is None:
                    c1, c2 = chrom_candidates(v.chrom)
                    need_nc[v.key] = (c1, c2, v.pos, v.ref, v.alt)
            cache[v.key] = data

        if need_nc:
            self._batch_noncancer(cache, need_nc)

        self._cache = cache

    def _batch_noncancer(
        self,
        cache: dict[str, GnomADData],
        need_nc: dict[str, tuple[str, str, int, str, str]],
    ) -> None:
        """Batch non-cancer-subset backfill: one JOIN mirroring _query_noncancer.

        Only variants whose af_non_cancer / faf95_non_cancer are still None are
        passed in (non_cancer_queried is already set by the caller). Per-field MAX
        over matching rows matches the single-variant path exactly."""
        try:
            import duckdb
            con = duckdb.connect(str(self._noncancer_db_path), read_only=True)
            cols = {
                r[1] for r in con.execute(
                    "PRAGMA table_info('non_cancer')"
                ).fetchall()
            }
            faf_expr = (
                "faf95_non_cancer" if "faf95_non_cancer" in cols
                else "NULL AS faf95_non_cancer"
            )
            con.execute(
                "CREATE TEMP TABLE _nckeys "
                "(vkey TEXT, chrom TEXT, pos BIGINT, ref TEXT, alt TEXT)"
            )
            key_rows: list[tuple] = []
            for vkey, (c1, c2, pos, ref, alt) in need_nc.items():
                key_rows.append((vkey, c1, pos, ref, alt))
                key_rows.append((vkey, c2, pos, ref, alt))
            con.executemany("INSERT INTO _nckeys VALUES (?, ?, ?, ?, ?)", key_rows)
            rows = con.execute(
                f"""
                SELECT k.vkey, n.af_non_cancer, {faf_expr}
                FROM _nckeys k
                JOIN non_cancer n
                  ON n.chrom = k.chrom AND n.pos = k.pos
                 AND n.ref = k.ref AND n.alt = k.alt
                """
            ).fetchall()
            con.close()
        except Exception as exc:
            log.error("gnomad_noncancer_query_error", error=str(exc))
            return

        nc_by_key: dict[str, list] = {}
        for r in rows:
            nc_by_key.setdefault(r[0], []).append((r[1], r[2]))
        for vkey, matched in nc_by_key.items():
            af_vals = [a for a, _ in matched if a is not None]
            faf_vals = [f for _, f in matched if f is not None]
            af_nc = max(af_vals) if af_vals else None
            faf_nc = max(faf_vals) if faf_vals else None
            data = cache[vkey]
            update: dict = {}
            if data.af_non_cancer is None and af_nc is not None:
                update["af_non_cancer"] = af_nc
            if data.faf95_non_cancer is None and faf_nc is not None:
                update["faf95_non_cancer"] = faf_nc
            if update:
                cache[vkey] = data.model_copy(update=update)

    def cached(self, variant: "VariantRecord") -> Optional[GnomADData]:
        """Return the precomputed record for *variant*.

        Falls back to a live query() when precompute() was not run or missed the
        key (the single-variant explain path never precomputes)."""
        hit = self._cache.get(variant.key)
        if hit is not None:
            return hit
        return self.query(variant.chrom, variant.pos, variant.ref, variant.alt)

    def _query_noncancer(
        self, c1: str, c2: str, pos: int, ref: str, alt: str
    ) -> Optional[tuple[Optional[float], Optional[float]]]:
        """Non-cancer-subset (AF, popmax FAF95) from the companion v3.1.2 DB
        (GRCh38), or None when the variant is absent / the DB is missing.

        Returns the per-field MAX over any matching rows (mirroring the main
        merge). ``faf95_non_cancer`` is absent in companion DBs built before that
        column was added — a schema probe degrades it to None there (BA1/BS1 then
        fall back to the overall FAF95). Connection-per-call mirrors query()'s
        pattern — DuckDB read-only opens are cheap and keep the object stateless
        across threads/processes."""
        if not self._noncancer_db_path.exists():
            log.warning("gnomad_noncancer_db_missing",
                        path=str(self._noncancer_db_path))
            return None
        try:
            import duckdb
            con = duckdb.connect(str(self._noncancer_db_path), read_only=True)
            cols = {
                r[1] for r in con.execute(
                    "PRAGMA table_info('non_cancer')"
                ).fetchall()
            }
            faf_expr = (
                "faf95_non_cancer" if "faf95_non_cancer" in cols
                else "NULL AS faf95_non_cancer"
            )
            rows = con.execute(
                f"""
                SELECT af_non_cancer, {faf_expr} FROM non_cancer
                WHERE chrom IN (?, ?) AND pos = ? AND ref = ? AND alt = ?
                """,
                [c1, c2, pos, ref, alt],
            ).fetchall()
            con.close()
        except Exception as exc:
            log.error("gnomad_noncancer_query_error", error=str(exc))
            return None
        if not rows:
            return None
        af_vals = [r[0] for r in rows if r[0] is not None]
        faf_vals = [r[1] for r in rows if r[1] is not None]
        af_nc = max(af_vals) if af_vals else None
        faf_nc = max(faf_vals) if faf_vals else None
        if af_nc is None and faf_nc is None:
            return None
        return (af_nc, faf_nc)

    def _columns(self, con) -> set[str]:
        """The variants table's column names (cached per instance). Used to
        degrade gracefully against DBs built before af_xy / ac_xx / nhomalt_xx."""
        if self._cols is None:
            self._cols = {
                r[1] for r in con.execute("PRAGMA table_info('variants')").fetchall()
            }
        return self._cols

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
