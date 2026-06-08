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
        with patch("subprocess.run", side_effect=fake_run):
            try:
                pred.precompute([v])
            except RuntimeError:
                pass
        cmd = captured["cmd"]
        assert "-R" in cmd and str(ref) in cmd
        # -A must be the explicit annotation FILE path, not the bare keyword.
        a_idx = cmd.index("-A")
        assert cmd[a_idx + 1] == str(ann)
        assert "grch38" != cmd[a_idx + 1]
        # model-type flag must be the long form accepted by the CLI.
        assert "--model-type" in cmd
