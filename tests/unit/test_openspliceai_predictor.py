"""OpenSpliceAI predictor wiring: reference-genome requirement + command build."""
from pathlib import Path
from unittest.mock import patch

from acmg_classifier.local_db.splice.openspliceai_predictor import OpenSpliceAIPredictor
from acmg_classifier.models.enums import Assembly
from acmg_classifier.models.variant import VariantRecord


def _files(tmp_path: Path) -> tuple[Path, Path, Path]:
    model = tmp_path / "model"
    model.mkdir()
    ref = tmp_path / "GRCh38.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    ann = tmp_path / "grch38.txt"
    ann.write_text("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n",
                   encoding="utf-8")
    return model, ref, ann


class TestAvailability:
    def test_unavailable_without_ref_genome(self, tmp_path):
        # `openspliceai variant` REQUIRES -R; without a reference FASTA the tool
        # must report unavailable rather than failing silently per-variant.
        model, _, ann = _files(tmp_path)
        with patch("shutil.which", return_value="/usr/bin/openspliceai"):
            pred = OpenSpliceAIPredictor(model, Assembly.GRCH38, 2000,
                                         ref_genome=None, annotation_file=ann)
        assert not pred.is_available()

    def test_unavailable_without_annotation(self, tmp_path):
        # The -A keyword is broken from an installed package, so an explicit
        # annotation file is required too.
        model, ref, _ = _files(tmp_path)
        with patch("shutil.which", return_value="/usr/bin/openspliceai"):
            pred = OpenSpliceAIPredictor(model, Assembly.GRCH38, 2000,
                                         ref_genome=ref, annotation_file=None)
        assert not pred.is_available()

    def test_available_with_all_inputs(self, tmp_path):
        model, ref, ann = _files(tmp_path)
        with patch("shutil.which", return_value="/usr/bin/openspliceai"):
            pred = OpenSpliceAIPredictor(model, Assembly.GRCH38, 2000,
                                         ref_genome=ref, annotation_file=ann)
        assert pred.is_available()


class TestCommand:
    def test_precompute_passes_required_flags(self, tmp_path):
        from acmg_classifier.utils import progress
        model, ref, ann = _files(tmp_path)
        with patch("shutil.which", return_value="/usr/bin/openspliceai"):
            pred = OpenSpliceAIPredictor(model, Assembly.GRCH38, 2000,
                                         ref_genome=ref, annotation_file=ann)

        captured = {}

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            raise RuntimeError("stop after capture")  # avoid real subprocess

        v = VariantRecord(chrom="chr17", pos=43115780, ref="C", alt="T",
                          assembly=Assembly.GRCH38)
        progress.set_enabled(False)  # force the simple (no-polling) subprocess path
        try:
            with patch("subprocess.run", side_effect=fake_run):
                try:
                    pred.precompute([v])
                except RuntimeError:
                    pass
        finally:
            progress.set_enabled(None)
        cmd = captured["cmd"]
        assert "-R" in cmd and str(ref) in cmd
        # -A must be the explicit annotation FILE path, not the bare keyword.
        a_idx = cmd.index("-A")
        assert cmd[a_idx + 1] == str(ann)
        assert "grch38" != cmd[a_idx + 1]
        # model-type flag must be the long form accepted by the CLI.
        assert "--model-type" in cmd


class TestProgressPath:
    def test_progress_branch_polls_output_and_caches(self, tmp_path):
        # With progress enabled, precompute runs via Popen and drives the bar by
        # polling the output VCF. A fake Popen writes one record so we verify the
        # branch runs end-to-end and the score is cached.
        from acmg_classifier.utils import progress
        model, ref, ann = _files(tmp_path)
        with patch("shutil.which", return_value="/usr/bin/openspliceai"):
            pred = OpenSpliceAIPredictor(model, Assembly.GRCH38, 80,
                                         ref_genome=ref, annotation_file=ann)

        class FakeProc:
            returncode = 0

            def __init__(self, cmd, **kwargs):
                out = cmd[cmd.index("-O") + 1]
                with open(out, "w") as f:
                    f.write("##fileformat=VCFv4.2\n")
                    f.write("chr17\t43115780\t.\tC\tT\t.\t.\t"
                            "OpenSpliceAI=T|BRCA1|0.91|0.0|0.0|0.0\n")

            def wait(self, timeout=None):
                return 0  # finishes immediately

        v = VariantRecord(chrom="chr17", pos=43115780, ref="C", alt="T",
                          assembly=Assembly.GRCH38)
        progress.set_enabled(True)
        try:
            with patch("subprocess.Popen", FakeProc):
                pred.precompute([v])
        finally:
            progress.set_enabled(None)

        score = pred.predict(v)
        assert score.is_available
        assert score.max_delta == 0.91
