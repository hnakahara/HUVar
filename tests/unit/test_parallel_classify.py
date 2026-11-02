"""Speed-ups for large VCF batches: process-parallel classification, concurrent
VEP chunks, VEP / OpenSpliceAI overlap and the per-evaluator timing report.

The parallel path must give exactly the in-process result, in input order."""
import sys
import threading
import time
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from acmg_classifier.config import Config
from acmg_classifier.criteria.registry import CriteriaRegistry
from acmg_classifier.classification.classifier_2015 import Classifier2015
from acmg_classifier.classification.classifier_bayesian import ClassifierBayesian
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, ESM1bData, GnomADData,
)
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pipeline import pipeline as pl

_LINUX = sys.platform.startswith("linux")


def _variants_and_annotations(n=320):
    variants, anns = [], {}
    genes = ["BRCA1", "BRCA2", "MECP2", "TP53"]
    csqs = [ConsequenceType.MISSENSE, ConsequenceType.SYNONYMOUS,
            ConsequenceType.STOP_GAINED, ConsequenceType.INTRON]
    for i in range(n):
        v = VariantRecord(chrom="chr17", pos=43_000_000 + i, ref="G", alt="ACT"[i % 3],
                          assembly=Assembly.GRCH38)
        variants.append(v)
        if i % 37 == 5:  # a few records VEP "dropped"
            continue
        gene, csq = genes[i % 4], csqs[(i // 4) % 4]
        pos = 10 + i // 3
        anns[v.key] = AnnotationData(
            gnomad=GnomADData(af=0.0 if i % 5 else 0.02, ac=0, an=100_000),
            consequences=[ConsequenceInfo(
                transcript_id=f"NM_{gene}.1", gene_id="ENSG", gene_symbol=gene,
                consequence=csq, biotype="protein_coding", is_mane_select=True,
                protein_position=pos, amino_acids="R/C",
                hgvs_c=f"NM_{gene}.1:c.{pos * 3}G>A",
                hgvs_p=f"NP_{gene}.1:p.Arg{pos}Cys")],
            esm1b=ESM1bData(llr=-12.0 + (i % 9) * 1.5),
        )
    return variants, anns


def _ctx(cfg, anns):
    return pl._ClassifyContext(cfg, CriteriaRegistry(cfg), Classifier2015(),
                               ClassifierBayesian(), anns, {})


def _dump(r):
    return r.model_dump(exclude={"annotation"}, mode="json")


@pytest.mark.skipif(not _LINUX, reason="parallel classification needs fork (Linux)")
def test_parallel_classification_matches_in_process():
    cfg = Config(data_dir=Path("./data"), workers=3)
    variants, anns = _variants_and_annotations()
    assert pl.classify_process_count(cfg, len(variants)) == 3

    serial = pl._classify_all(variants, _ctx(cfg, anns), serial=True)
    parallel = pl._classify_all(variants, _ctx(cfg, anns))

    assert [r.variant_id for r in parallel] == [v.key for v in variants]
    assert [_dump(r) for r in parallel] == [_dump(r) for r in serial]
    # The parent's AnnotationData is re-attached (not a pickled copy).
    for v, r in zip(variants, parallel):
        if v.key in anns:
            assert r.annotation is anns[v.key]
        else:
            assert isinstance(r.annotation, AnnotationData)
            assert "VEP returned no annotation for this variant" in r.warnings


def test_process_count_rules():
    cfg = Config(data_dir=Path("./data"), workers=8)
    assert pl.classify_process_count(cfg, 10) == 1          # too small to pay off
    one = Config(data_dir=Path("./data"), workers=8, classify_workers=1)
    assert pl.classify_process_count(one, 50_000) == 1      # explicitly in-process
    if _LINUX:
        assert pl.classify_process_count(cfg, 50_000) == 8
        capped = Config(data_dir=Path("./data"), workers=36, classify_workers=4)
        assert pl.classify_process_count(capped, 50_000) == 4
        # never more workers than ~50 variants each
        assert pl.classify_process_count(cfg, 250) == 5


def test_warmup_picks_first_record_of_each_gene():
    variants, anns = _variants_and_annotations(40)
    idx = pl._warmup_indices(variants, anns)
    genes = [anns[variants[i].key].primary_consequence.gene_symbol
             if variants[i].key in anns else None for i in idx]
    assert len(genes) == len(set(genes))
    assert {"BRCA1", "BRCA2", "MECP2", "TP53"} <= set(genes)


def test_registry_timing_report():
    cfg = Config(data_dir=Path("./data"))
    variants, anns = _variants_and_annotations(12)
    reg = CriteriaRegistry(cfg)
    assert reg.timings() == {}
    reg.enable_timing()
    for v in variants:
        reg.evaluate_all(v, anns.get(v.key, AnnotationData()))
    t = reg.timings()
    assert "PVS1Evaluator" in t and "combination rules" in t
    assert list(t.values()) == sorted(t.values(), reverse=True)
    text = pl._format_timings(t, len(variants))
    assert "ms/variant" in text and "PVS1Evaluator" in text


# --- VEP chunking -------------------------------------------------------------

def _runner(workers):
    from acmg_classifier.local_db.vep_runner import LocalVEPRunner
    r = LocalVEPRunner(vep_cmd="vep", cache_dir=Path("."), fasta=Path("x.fa"),
                       assembly="GRCh38", workers=workers)
    calls, lock = [], threading.Lock()

    def fake_run(chunk, forks=None):
        with lock:
            calls.append((len(chunk), forks))
        time.sleep(0.01)
        return {v.key: [] for v in chunk}

    r._run_vep = fake_run
    return r, calls


def _vs(n):
    return [VariantRecord(chrom="chr1", pos=1000 + i, ref="A", alt="G",
                          assembly=Assembly.GRCH38) for i in range(n)]


def test_vep_process_count():
    from acmg_classifier.local_db.vep_runner import vep_process_count
    assert vep_process_count(1, 36) == 1
    assert vep_process_count(11, 36) == 4      # workers // 8
    assert vep_process_count(11, 4) == 1
    assert vep_process_count(2, 36) == 2       # never more than the chunks
    assert vep_process_count(11, 36, processes=6) == 6


def test_vep_chunks_run_concurrently_with_even_sizes_and_split_forks():
    r, calls = _runner(36)
    out = r.annotate_batch(_vs(53_343), batch_size=5000)
    assert len(out) == 53_343
    assert len(calls) == 12                     # 11 rounded up to a multiple of 4
    assert {f for _, f in calls} == {9}         # 36 forks shared by 4 processes
    sizes = sorted(n for n, _ in calls)
    assert sizes[-1] - sizes[0] <= sizes[-1] // 2 and sum(sizes) == 53_343


def test_vep_single_chunk_keeps_all_forks():
    r, calls = _runner(8)
    r.annotate_batch(_vs(300), batch_size=5000)
    assert calls == [(300, 8)]
    r2, calls2 = _runner(4)
    r2.annotate_batch(_vs(1200), batch_size=500)
    assert [f for _, f in calls2] == [4, 4, 4] and sum(n for n, _ in calls2) == 1200


# --- VEP / precompute overlap ---------------------------------------------------

def _orchestrator(overlap):
    from acmg_classifier.annotation.orchestrator import AnnotationOrchestrator
    o = object.__new__(AnnotationOrchestrator)
    o._cfg = Config(data_dir=Path("./data"), workers=2, overlap_precompute=overlap)
    events, lock = [], threading.Lock()

    def note(name):
        with lock:
            events.append((name, threading.current_thread().name))

    def vep(variants, batch_size, processes=0):
        note("vep_start")
        time.sleep(0.2)
        note("vep_end")
        return {v.key: [] for v in variants}

    o._vep = MagicMock()
    o._vep.annotate_batch.side_effect = vep
    o._gnomad = MagicMock()
    o._gnomad.precompute.side_effect = lambda vs: note("gnomad")
    o._splice = MagicMock()
    o._splice.precompute.side_effect = lambda vs: note("splice")
    o._annotate_one = lambda v, csq: AnnotationData()
    return o, events


def test_precompute_overlaps_vep():
    o, events = _orchestrator(True)
    out = o.annotate_batch(_vs(3))
    assert len(out) == 3
    names = [e for e, _ in events]
    # splice/gnomAD ran while VEP was still running, on a helper thread
    assert names.index("splice") < names.index("vep_end")
    assert dict(events)["splice"] != dict(events)["vep_start"]
    kwargs = o._vep.annotate_batch.call_args.kwargs
    assert kwargs["batch_size"] == 5000 and kwargs["processes"] == 0


def test_precompute_sequential_when_disabled():
    o, events = _orchestrator(False)
    o.annotate_batch(_vs(3))
    assert [e for e, _ in events] == ["vep_start", "vep_end", "gnomad", "splice"]


def test_precompute_error_is_raised():
    o, _ = _orchestrator(True)
    o._splice.precompute.side_effect = RuntimeError("osai broke")
    with pytest.raises(RuntimeError, match="osai broke"):
        o.annotate_batch(_vs(3))


# --- database handles ------------------------------------------------------------

def test_clinvar_exists_cache_is_positive_only(tmp_path):
    from acmg_classifier.local_db import clinvar_sqlite as cs
    p = tmp_path / "clinvar.db"
    assert cs._db_exists(p) is False
    p.write_bytes(b"")
    assert cs._db_exists(p) is True            # picked up once it exists
    assert str(p) in cs._EXISTING


def test_clinvar_connections_dropped_after_fork():
    from acmg_classifier.local_db import clinvar_sqlite as cs
    sentinel = object()
    cs._CONN_CACHE["x"] = sentinel
    try:
        cs._reset_after_fork()
        assert "x" not in cs._CONN_CACHE and sentinel in cs._FORK_ORPHANS
    finally:
        cs._FORK_ORPHANS.remove(sentinel)


def test_coverage_connection_reused_and_reopened_after_close_all(tmp_path):
    duckdb = pytest.importorskip("duckdb")
    from acmg_classifier.local_db import coverage_db
    db = tmp_path / "cov.duckdb"
    con = duckdb.connect(str(db))
    con.execute("CREATE TABLE coverage (chrom TEXT, pos BIGINT, mean_dp DOUBLE)")
    con.execute("INSERT INTO coverage VALUES ('chr17', 100, 30.0), ('chr17', 101, 10.0)")
    con.close()
    cov = coverage_db.CoverageDB(db)
    assert cov.mean_depth("17", 100) == 30.0
    first = cov._connection()
    assert cov._connection() is first           # one connection, not one per call
    coverage_db.close_all()                     # what happens right before fork()
    assert cov.mean_depth("17", 100, 101) == 20.0
    assert cov._connection() is not first


# --- run_pipeline wiring -----------------------------------------------------------

def _run(tmp_path, monkeypatch, name, **cfg_kw):
    variants, anns = _variants_and_annotations(260)
    vcf = tmp_path / "in.vcf"
    with vcf.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in variants:
            fh.write(f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\tPASS\t.\n")

    class _Orch:
        def __init__(self, cfg):
            pass

        def annotate_batch(self, vs):
            return {v.key: anns[v.key] for v in vs if v.key in anns}

    import acmg_classifier.annotation.orchestrator as orch
    monkeypatch.setattr(orch, "AnnotationOrchestrator", _Orch)
    out = tmp_path / f"{name}.tsv"
    cfg = Config(data_dir=Path("./data"), assembly=Assembly.GRCH38, **cfg_kw)
    res = pl.run_pipeline(vcf, cfg, output_path=out,
                          profile_path=cfg_kw.pop("_profile", None))
    return res, out.read_text()


@pytest.mark.skipif(not _LINUX, reason="parallel classification needs fork (Linux)")
def test_run_pipeline_parallel_output_identical(tmp_path, monkeypatch):
    r1, tsv1 = _run(tmp_path, monkeypatch, "serial", workers=3, classify_workers=1)
    r3, tsv3 = _run(tmp_path, monkeypatch, "parallel", workers=3)
    assert len(r1) == len(r3) == 260
    assert tsv1 == tsv3


def test_run_pipeline_profile_prints_evaluator_times(tmp_path, monkeypatch, capsys):
    variants, anns = _variants_and_annotations(30)
    vcf = tmp_path / "in.vcf"
    with vcf.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in variants:
            fh.write(f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\tPASS\t.\n")

    class _Orch:
        def __init__(self, cfg):
            pass

        def annotate_batch(self, vs):
            return {v.key: anns[v.key] for v in vs if v.key in anns}

    import acmg_classifier.annotation.orchestrator as orch
    monkeypatch.setattr(orch, "AnnotationOrchestrator", _Orch)
    pl.run_pipeline(vcf, Config(data_dir=Path("./data"), assembly=Assembly.GRCH38),
                    output_path=tmp_path / "o.tsv", profile_path=tmp_path / "p.prof")
    err = capsys.readouterr().err
    assert "criteria evaluation time per evaluator" in err and "PM2Evaluator" in err
