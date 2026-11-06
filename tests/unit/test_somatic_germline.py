"""`somatic` reports the germline ACMG/AMP classification of `classify` next to the
oncogenicity (germline_* columns), from the same annotation, and accepts the
classify options that shape it."""
import csv
from pathlib import Path

from click.testing import CliRunner

from acmg_classifier.config import Config
from acmg_classifier.io.tsv_writer import GERMLINE_HEADER, _HEADER
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import Assembly, ConsequenceType as CT
from acmg_classifier.models.variant import SampleCall, VariantRecord
from acmg_classifier.pipeline.pipeline import classify_annotated
from acmg_classifier.somatic.pipeline import SomaticClassifier
from acmg_classifier.somatic.report import HEADER, write_somatic_tsv

_RES = Path(__file__).resolve().parents[2] / "resources" / "shared"


def _cfg(**kw):
    return Config(data_dir=Path("./nonexistent-data"), **kw)


def _ann():
    c = ConsequenceInfo(transcript_id="NM_004333.6", gene_id="g", gene_symbol="BRAF",
                        consequence=CT.MISSENSE, biotype="protein_coding", is_mane_select=True,
                        protein_position=600, amino_acids="V/E",
                        hgvs_p="NP_004324.2:p.Val600Glu")
    return AnnotationData(consequences=[c], gnomad=GnomADData(af=0.0, ac=0, an=0))


def _v():
    return VariantRecord(chrom="chr7", pos=140753336, ref="A", alt="T",
                         assembly=Assembly.GRCH38)


def _rows(path):
    with Path(path).open(encoding="utf-8") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def test_germline_block_follows_somatic_columns(tmp_path):
    cfg, v, ann = _cfg(), _v(), _ann()
    g = classify_annotated(v, ann, cfg)
    r = SomaticClassifier(cfg).classify(v, ann, germline=g.criteria_results)
    r.germline = g
    r.sample_calls = [SampleCall(sample="T1", gt="0/1", dp=100, vaf=0.25),
                      SampleCall(sample="T2", gt="0/1", dp=80, vaf=0.1)]
    out = tmp_path / "s.tsv"
    write_somatic_tsv([r], out)
    head, *rows = _rows(out)
    assert head == HEADER + GERMLINE_HEADER
    assert GERMLINE_HEADER[:len(_HEADER) - 15] == ["germline_" + h for h in _HEADER[15:]]
    assert "germline_classification_2015" in head and "germline_PS3_evidence" in head
    assert len(rows) == 2 and all(len(x) == len(head) for x in rows)
    by = dict(zip(head, rows[0]))
    assert by["germline_classification_2015"] == g.classification_2015.value
    assert by["germline_classification_bayesian"] == g.classification_bayesian.value
    assert by["sample"] == "T1" and by["oncogenicity"] == r.classification.value


def test_without_germline_header_unchanged(tmp_path):
    cfg = _cfg()
    r = SomaticClassifier(cfg).classify(_v(), _ann())
    out = tmp_path / "s.tsv"
    write_somatic_tsv([r], out)
    assert _rows(out)[0] == HEADER
    # Requested explicitly but missing for a row: empty germline cells.
    write_somatic_tsv([r], out, germline=True)
    head, row = _rows(out)
    assert head == HEADER + GERMLINE_HEADER and row[len(HEADER):] == [""] * len(GERMLINE_HEADER)


def test_germline_criteria_reused(monkeypatch):
    cfg, v, ann = _cfg(), _v(), _ann()
    clf = SomaticClassifier(cfg)
    fresh = clf.classify(v, ann)
    g = classify_annotated(v, ann, cfg)

    def boom(*a, **k):
        raise AssertionError("germline criteria evaluated twice")
    monkeypatch.setattr(clf._registry, "evaluate_all", boom)
    reused = clf.classify(v, ann, germline=g.criteria_results)
    assert reused.met_criteria == fresh.met_criteria and reused.score == fresh.score


def _vcf(tmp_path):
    p = tmp_path / "in.vcf"
    p.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr7,length=159345973>\n"
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                 '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT1\n"
                 "chr7\t140753336\t.\tA\tT\t.\tPASS\t.\tGT:AD\t0/1:75,25\n"
                 "chr7\t140753337\t.\tC\tG\t.\tPASS\t.\tGT:AD\t0/1:90,10\n")
    return p


def _fake_annotation(monkeypatch):
    from acmg_classifier.annotation import orchestrator
    monkeypatch.setattr(orchestrator.AnnotationOrchestrator, "__init__",
                        lambda self, cfg: None)
    monkeypatch.setattr(orchestrator.AnnotationOrchestrator, "annotate_batch",
                        lambda self, vs: {x.key: _ann() for x in vs})


def test_pipeline_writes_germline_and_limit(tmp_path, monkeypatch):
    from acmg_classifier.somatic.pipeline import run_somatic_pipeline
    _fake_annotation(monkeypatch)
    out = tmp_path / "s.tsv"
    res = run_somatic_pipeline(_vcf(tmp_path), _cfg(), output_path=out, limit=1)
    assert len(res) == 1 and res[0].germline is not None
    head, row = _rows(out)
    by = dict(zip(head, row))
    assert by["germline_classification_2015"] == res[0].germline.classification_2015.value
    assert by["vaf"] == "0.25" and by["sample"] == "T1"

    run_somatic_pipeline(_vcf(tmp_path), _cfg(), output_path=out, germline=False)
    head, *rows = _rows(out)
    assert head == HEADER and len(rows) == 2


def test_germline_supplement_only_affects_germline(tmp_path, monkeypatch):
    from acmg_classifier.pipeline import pipeline as pl
    from acmg_classifier.somatic.pipeline import run_somatic_pipeline
    _fake_annotation(monkeypatch)
    seen = {}
    real = pl.classify_variants

    def spy(variants, annotations, cfg, supplement=None):
        seen["curated"], seen["sup"] = cfg.curated_evidence, supplement
        return real(variants, annotations, cfg, supplement)
    monkeypatch.setattr(pl, "classify_variants", spy)
    sup = tmp_path / "sup.tsv"
    sup.write_text("variant_id\tcriterion\tstrength\tevidence\n"
                   "chr7:140753336:A:T\tPS3\tStrong\tPMID 1\n")
    out = tmp_path / "s.tsv"
    res = run_somatic_pipeline(_vcf(tmp_path), _cfg(), output_path=out,
                               supplement_path=sup)
    assert seen["curated"] is True and "chr7:140753336:A:T" in seen["sup"]
    g = res[0].germline
    ps3 = next(c for c in g.criteria_results if c.criterion.value == "PS3")
    assert ps3.triggered and "PMID 1" in ps3.evidence
    # The oncogenicity is evaluated without the curator entries.
    plain = SomaticClassifier(_cfg()).classify(_v(), _ann())
    assert res[0].met_criteria == plain.met_criteria


def test_cli_passes_classify_options(tmp_path, monkeypatch):
    from acmg_classifier import cli as cli_mod
    from acmg_classifier.somatic import pipeline as sp
    got = {}
    monkeypatch.setattr(sp, "run_somatic_pipeline",
                        lambda vcf, cfg, **kw: got.update(cfg=cfg, **kw))
    data = tmp_path / "data"
    (data / "shared").mkdir(parents=True)
    (data / "shared" / "disease_prevalence.tsv").write_bytes(
        (_RES / "disease_prevalence.tsv").read_bytes())
    sup = tmp_path / "sup.tsv"
    sup.write_text("variant_id\tcriterion\tstrength\tevidence\n")
    osup = tmp_path / "onc.tsv"
    osup.write_text("variant_id\tcriterion\tevidence\nchr7:1:A:T\tOP2\tx\n")
    r = CliRunner().invoke(cli_mod.cli, [
        "somatic", str(_vcf(tmp_path)), "--data-dir", str(data), "--assembly", "GRCh38",
        "--insilico-tool", "revel", "--with-bayesdel", "--with-cadd",
        "--cspec", "RYR1=congenital_myopathy_ar", "--exclude-self-expert-panel",
        "--supplement", str(sup), "--supplement-mode", "manual-only",
        "--somatic-supplement", str(osup),
        "--workers", "6", "--classify-workers", "2", "--limit", "5", "--no-progress"])
    assert r.exit_code == 0, r.output
    cfg = got["cfg"]
    assert cfg.use_bayesdel and cfg.use_cadd and cfg.exclude_self_expert_panel
    assert cfg.cspec_selection == {"RYR1": "congenital_myopathy_ar"}
    assert cfg.supplement_mode.value == "manual-only"
    assert cfg.workers == 6 and cfg.classify_workers == 2
    assert got["limit"] == 5 and got["germline"] is True
    assert got["supplement_path"] == sup and got["somatic_supplement_path"] == osup

    r = CliRunner().invoke(cli_mod.cli, [
        "somatic", str(_vcf(tmp_path)), "--data-dir", str(data), "--no-germline",
        "--supplement", str(sup)])
    assert r.exit_code != 0 and "no-germline" in r.output

    # A TSV given to the wrong option points at the right one.
    r = CliRunner().invoke(cli_mod.cli, [
        "somatic", str(_vcf(tmp_path)), "--data-dir", str(data), "--supplement", str(osup)])
    assert r.exit_code != 0 and "--somatic-supplement" in r.output
    acmg = tmp_path / "acmg.tsv"
    acmg.write_text("variant_id\tcriterion\tstrength\tevidence\nchr7:1:A:T\tPS3\tStrong\tx\n")
    r = CliRunner().invoke(cli_mod.cli, [
        "somatic", str(_vcf(tmp_path)), "--data-dir", str(data),
        "--somatic-supplement", str(acmg)])
    assert r.exit_code != 0 and "pass it with --supplement" in r.output


def test_cli_somatic_explain_evidence_options(tmp_path, monkeypatch):
    from acmg_classifier import cli as cli_mod
    from acmg_classifier.somatic import pipeline as sp
    got = {}
    monkeypatch.setattr(sp, "run_somatic_single",
                        lambda *a, **kw: got.update(kw))
    key = "chr7:140753336:A:T"
    osup = tmp_path / "onc.tsv"
    osup.write_text(f"variant_id\tcriterion\tevidence\n{key}\tOS2\tPMID 2\n")
    sup = tmp_path / "acmg.tsv"
    sup.write_text(f"variant_id\tcriterion\tstrength\tevidence\n{key}\tPM1\tModerate\tx\n")
    base = ["somatic-explain", "chr7", "140753336", "A", "T", "--data-dir", str(tmp_path)]
    r = CliRunner().invoke(cli_mod.cli, base + [
        "--evidence", "PS3:strong", "--supplement", str(sup),
        "--somatic-evidence", "OP2:retinoblastoma", "--somatic-supplement", str(osup)])
    assert r.exit_code == 0, r.output
    assert {e.criterion.value for e in got["supplement"]} == {"PS3", "PM1"}
    assert {e.criterion.value for e in got["manual"]} == {"OS2", "OP2"}
    assert got["germline"] is True

    r = CliRunner().invoke(cli_mod.cli, base + ["--evidence", "OP2:retinoblastoma"])
    assert r.exit_code != 0 and "--somatic-evidence" in r.output
    r = CliRunner().invoke(cli_mod.cli, base + ["--no-germline", "--evidence", "PS3:strong"])
    assert r.exit_code != 0 and "no-germline" in r.output


def test_somatic_single_prints_germline_report(monkeypatch, capsys):
    from acmg_classifier.pipeline import pipeline as pl
    from acmg_classifier.somatic.pipeline import run_somatic_single
    monkeypatch.setattr(pl, "annotate_one", lambda *a, **k: (_v(), _ann()))
    r = run_somatic_single("chr7", 140753336, "A", "T", _cfg(), germline=True)
    out = capsys.readouterr().out
    assert r.germline is not None
    assert "Somatic Oncogenicity Report" in out and "ACMG Variant Classification Report" in out
    r = run_somatic_single("chr7", 140753336, "A", "T", _cfg(), print_result=False)
    assert r.germline is None
