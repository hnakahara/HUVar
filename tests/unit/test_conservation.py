"""PhyloPReader lazy-open behaviour (the BP7 conservation speedup)."""
import sys
import types

from acmg_classifier.local_db.conservation import PhyloPReader


def _install_fake_pybigwig(monkeypatch, score=1.5):
    """Inject a fake `pyBigWig` module that counts open() calls."""
    calls = {"open": 0}

    class _FakeBW:
        def values(self, chrom, start, end):
            return [score]

    mod = types.ModuleType("pyBigWig")

    def _open(path):
        calls["open"] += 1
        return _FakeBW()

    mod.open = _open
    monkeypatch.setitem(sys.modules, "pyBigWig", mod)
    return calls


def test_construction_does_not_open(monkeypatch, tmp_path):
    calls = _install_fake_pybigwig(monkeypatch)
    bw = tmp_path / "phylop.bw"
    bw.write_bytes(b"")
    reader = PhyloPReader(bw)
    # The ~9 GB file must NOT be opened just by constructing the reader.
    assert calls["open"] == 0
    # ...yet the gate reports available (openable on demand).
    assert reader.is_available()


def test_value_triggers_single_lazy_open(monkeypatch, tmp_path):
    calls = _install_fake_pybigwig(monkeypatch, score=2.7)
    bw = tmp_path / "phylop.bw"
    bw.write_bytes(b"")
    reader = PhyloPReader(bw)

    assert reader.value("chr1", 100) == 2.7
    assert calls["open"] == 1          # opened on first lookup
    assert reader.value("chr1", 200) == 2.7
    assert calls["open"] == 1          # not re-opened on subsequent lookups


def test_none_path_unavailable(monkeypatch):
    _install_fake_pybigwig(monkeypatch)
    reader = PhyloPReader(None)
    assert not reader.is_available()
    assert reader.value("chr1", 100) is None


def test_missing_pybigwig_unavailable(monkeypatch, tmp_path):
    # Simulate pyBigWig not importable: removing it makes `import pyBigWig` fail.
    monkeypatch.setitem(sys.modules, "pyBigWig", None)  # forces ImportError
    bw = tmp_path / "phylop.bw"
    bw.write_bytes(b"")
    reader = PhyloPReader(bw)
    assert not reader.is_available()
