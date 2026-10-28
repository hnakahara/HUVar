"""OS2 / SBS2 from the laboratory-curated functional table."""
import re
from pathlib import Path

from acmg_classifier.config import Config
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.somatic.evaluator import SomaticEvaluator
from acmg_classifier.somatic.lab_curated import load_curated
from acmg_classifier.somatic.models import OncCriterion as C, Oncogenicity, classify_points
from acmg_classifier.somatic.resources import load_gene_roles, resolve

A, S, CT = ACMGCriterion, CriterionStrength, ConsequenceType
_PACKAGED = Path(__file__).resolve().parents[2] / "resources" / "shared" / "somatic_functional_curated.tsv"
_HDR = "gene\tprotein_change\tcriterion\tpmids\tevidence\tverification\tstatus\tcurator\treviewed\tnote\n"


def _braf(cfg, germ=()):
    c = ConsequenceInfo(transcript_id="NM_004333.6", gene_id="g", gene_symbol="BRAF",
                        consequence=CT.MISSENSE, biotype="protein_coding", is_mane_select=True,
                        protein_position=600, amino_acids="V/E", hgvs_p="NP_004324.2:p.Val600Glu")
    ann = AnnotationData(consequences=[c], gnomad=GnomADData(af=1.24e-6, ac=2, an=1600000,
                                                             popmax_af=1.1e-5, popmax_pop="sas"))
    g = {x: CriteriaResult.not_met(x) for x in (A.PP3, A.BP4, A.BP7, A.PM1, A.PS3, A.BS3)}
    for x, st, ev in germ:
        g[x] = CriteriaResult.met(x, st, ev)
    v = VariantRecord(chrom="chr7", pos=140753336, ref="A", alt="T", assembly=Assembly.GRCH38)
    crit, _ = SomaticEvaluator(cfg).evaluate(v, ann, list(g.values()))
    return {r.criterion: r for r in crit}, sum(r.points for r in crit)


def _cfg_with(tmp_path, rows):
    cfg = Config(data_dir=tmp_path)
    cfg.somatic_functional_curated_tsv.parent.mkdir(parents=True, exist_ok=True)
    cfg.somatic_functional_curated_tsv.write_text(_HDR + "".join(r + "\n" for r in rows), encoding="utf-8")
    load_curated.cache_clear()
    return cfg


def test_packaged_table_is_well_formed():
    load_curated.cache_clear()
    db = load_curated(_PACKAGED)
    assert len(db) >= 60 and ("BRAF", "V600E") in db
    roles = load_gene_roles(resolve(None, "somatic_gene_roles.tsv"))
    key = re.compile(r"^[A-Z]\d+(?:_[A-Z]\d+)?(?:[A-Z]|del|dup|ins[A-Z]+)$")
    for (gene, change), row in db.items():
        assert key.match(change), change
        assert row.criterion == "OS2" and row.pmids and row.status == "draft"
        assert all(p.strip().isdigit() for p in row.pmids.split(","))
        assert gene in roles


def test_draft_row_is_shown_not_counted(tmp_path):
    load_curated.cache_clear()
    by, _ = _braf(Config(data_dir=tmp_path))          # packaged table, all draft
    assert not by[C.OS2].met and "not counted until status=approved" in by[C.OS2].evidence


def test_approved_row_gives_os2(tmp_path):
    cfg = _cfg_with(tmp_path, ["BRAF\tV600E\tOS2\t12068308,15035987\tKinase activation\tpmid\t"
                               "approved\tHN\t2026-10-28\t"])
    by, score = _braf(cfg, [(A.PP3, S.STRONG, "ESM1b")])
    assert by[C.OS2].active and "PMID 12068308" in by[C.OS2].evidence and "HN" in by[C.OS2].evidence
    assert score == 10 and classify_points(score) == Oncogenicity.ONCOGENIC


def test_merge_rules():
    from acmg_classifier.somatic.models import OncCriterionResult
    res = {c: OncCriterionResult(criterion=c) for c in C}
    res[C.SBS2].met, res[C.SBS2].evidence = True, "assay normal"
    st: dict = {}
    SomaticEvaluator._merge_functional(res, "OS2", "CIViC GoF", "CIViC", st)
    assert not res[C.OS2].met and not res[C.SBS2].met and st["conflict"]   # automated: cancel
    SomaticEvaluator._merge_functional(res, "OS2", "lab row", "lab-curated table", st,
                                       authoritative=True)
    assert res[C.OS2].met and not res[C.SBS2].met and not st["conflict"]   # lab row wins
