"""BS2 VCEP-applicability: cspec extraction for genes where BS2 cannot be
derived from gnomAD population counts.

Our BS2 evaluator is gnomAD-based, so a VCEP must resolve to ``not_applicable``
when it bars population data or scores BS2 purely on clinical phenotype/points
with no gnomAD-countable rule (RPE65, the RASopathy point specs, the Fanconi
"points per proband" cancer specs, ...).
"""
import importlib.util
from pathlib import Path

# --- load the build script as a module (it lives under scripts/, not the pkg) -
_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


def _bs2_code(desc, applicability="Applicable"):
    return {"label": "BS2", "evidenceStrengths": [
        {"label": "Strong", "applicability": applicability, "description": desc},
    ]}


def _rs(gene, *codes):
    return {"genes": [{"label": gene}], "criteriaCodes": list(codes)}


class TestBS2NotApplicable:
    def test_bars_population_data(self):
        # RASopathy GN004 wording.
        rs = _rs("NRAS", _bs2_code(
            "general population data should not be used for this criterion."))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_gnomad_not_considered(self):
        # RPE65 GN120 wording.
        rs = _rs("RPE65", _bs2_code(
            "Variant is present in >= 3 homozygotes without any features of the "
            "phenotype. Presence in databases such as gnomAD are not considered."))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_not_applicable_due_to(self):
        # KCNQ1 GN112 wording.
        rs = _rs("KCNQ1", _bs2_code("Not applicable due to incomplete penetrance."))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_pure_point_scoring(self):
        # Single-gene RASopathy specs (NRAS/GN039 etc.) inherit point-based,
        # phenotype-driven scoring with no gnomAD-countable rule.
        rs = _rs("NRAS", _bs2_code("\\-4 Points."))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_points_per_proband_phenotype(self):
        # BRCA1/2, PALB2: Fanconi-Anemia "points per proband", no gnomAD count.
        rs = _rs("BRCA1", _bs2_code(
            "Applied in absence of features of recessive disease (Fanconi Anemia "
            "phenotype). Approach to assign points per proband. BS2 = >= 4 points"))
        assert bdt._bs2_applicability(rs) == "not_applicable"

    def test_no_applicable_strength_declined(self):
        rs = _rs("FOO", _bs2_code("", applicability="Not Applicable"))
        assert bdt._bs2_applicability(rs) == "not_applicable"


class TestBS2StaysApplicable:
    def test_points_with_homozygous_gnomad_path(self):
        # APC: points OR a homozygous-state count — gnomAD-derivable, keep on.
        rs = _rs("APC", _bs2_code(
            ">= 10 points for healthy individuals OR >= 2 times in homozygous state."))
        assert bdt._bs2_applicability(rs) == "applicable"

    def test_plain_count_based(self):
        # GP1BA-style: plain homozygote/heterozygote count, no clinical qualifier.
        rs = _rs("GP1BA", _bs2_code(
            "Variant is identified in >= 1 homozygous or >= 2 heterozygous "
            "healthy adults (unrelated)."))
        assert bdt._bs2_applicability(rs) == "applicable"

    def test_hemizygotes_in_gnomad(self):
        # PDHA1: ">= 16 hemizygotes in gnomAD" — explicitly gnomAD-countable.
        rs = _rs("PDHA1", _bs2_code(
            "Observed in at least two healthy male adults ... AND/OR >= 16 "
            "hemizygotes in gnomAD"))
        assert bdt._bs2_applicability(rs) == "applicable"

    def test_no_bs2_code(self):
        rs = _rs("FOO")
        assert bdt._bs2_applicability(rs) == ""


class TestBS2HomOnly:
    """Incomplete-penetrance dominant genes whose BS2 counts homozygotes only."""

    def test_homozygotes_only_flagged(self):
        # BMPR2/GN125 wording: BS2 defined purely on gnomAD homozygotes.
        rs = _rs("BMPR2", _bs2_code(
            "Observed in >=3 homozygotes in gnomAD controls or reported in the "
            "literature (healthy adult individuals)."))
        assert bdt._bs2_hom_only(rs) == "1"

    def test_heterozygote_mention_not_flagged(self):
        # PIK3R2/GN018: also allows heterozygous family members → not hom-only.
        rs = _rs("PIK3R2", _bs2_code(
            "Award BS2 if >=3 homozygotes present in gnomAD or >=3 heterozygous "
            "in well phenotyped family members."))
        assert bdt._bs2_hom_only(rs) == ""

    def test_carrier_mention_not_flagged(self):
        rs = _rs("GENE", _bs2_code(
            "Observed in >=2 homozygous or >=4 carrier healthy adults."))
        assert bdt._bs2_hom_only(rs) == ""

    def test_no_bs2_code(self):
        assert bdt._bs2_hom_only(_rs("FOO")) == ""


class TestBS2Count:
    """The per-gene BS2 count is the LOWEST (Supporting) threshold across all
    applicable strengths, so a binary BS2 fires at the gene's minimum bar."""

    @staticmethod
    def _tiered(*descs):
        return {"label": "BS2", "evidenceStrengths": [
            {"label": "Strong", "applicability": "Applicable", "description": d}
            for d in descs
        ]}

    def test_min_across_strengths(self):
        # GUCY2D: Strong ">=6 ... OR ... >=3 homozytes (literature)", Supporting
        # ">=3 homozygotes in gnomAD" — the operative bar is 3, not 6.
        rs = _rs("GUCY2D", self._tiered(
            "Variant is present in >= 3 homozygotes without features ... "
            "Alternatively >= 6 homozygotes in gnomAD v.4.1.0 or later.",
            "Variant is present in >= 3 homozygotes in gnomAD v.4.1.0 or later.",
        ))
        assert bdt._bs2_count(rs) == "3"

    def test_spelled_out_count_ignored(self):
        # PDHA1: "at least two healthy male adults" ("two" is a word, not a
        # gnomAD digit count) must not lower the ">=16 hemizygotes" bar.
        rs = _rs("PDHA1", self._tiered(
            "Observed in at least two healthy male adults ... AND/OR >= 16 "
            "hemizygotes in gnomAD.",
        ))
        assert bdt._bs2_count(rs) == "16"

    def test_no_count(self):
        rs = _rs("PAH", _bs2_code(
            "Only to be used when variant is observed in the homozygous state "
            "in a healthy adult."))
        assert bdt._bs2_count(rs) == ""


class TestResolvedTSV:
    """Spot-check the committed TSV reflects the extraction."""

    def test_committed_tsv_marks_target_genes_not_applicable(self):
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        for gene in ("RPE65", "NRAS", "BRAF", "KCNQ1", "BRCA1", "BRCA2", "PALB2"):
            assert rows[gene]["bs2"] == "not_applicable", gene

    def test_committed_tsv_withholds_bs2_needing_clinical_confirmation(self):
        # Genes whose VCEP BS2 demands phenotype/lab/functional confirmation or
        # internal-only cohort data that gnomAD cannot supply — BS2 is withheld so
        # a population-count BS2 never falsely fires (and never double-counts the
        # gnomAD individuals that already drove BS1/BA1). See
        # build_disease_thresholds._BS2_CLINICAL_CONFIRMATION.
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        clinical_confirmation = (
            # batch 1
            "HNF4A", "RYR1", "LDLR", "GAA", "ITGA2B", "ITGB3",
            "CDH1", "UBE3A", "DICER1", "SERPINC1", "TP53",
            # batch 2 — phase / lab assay / specific phenotype
            "HNF1A", "GCK", "IDUA", "GP1BA", "GP1BB", "GP9", "VHL", "PTEN",
            "MLH1", "MSH2", "MSH6", "PMS2",
            "CDH23", "GJB2", "MYO6", "MYO7A", "SLC26A4", "TECTA", "USH2A",
            "MYO15A", "OTOF",
            # batch 2 — "healthy/unaffected adult" not gnomAD-confirmable
            "ACTA1", "DNM2", "NEB", "SCN1A", "SCN2A", "SCN3A", "SCN8A",
            "DCLRE1C", "GATM", "HBB", "HBA2",
            "FOXG1", "TCF4", "APC",
        )
        for gene in clinical_confirmation:
            assert rows[gene]["bs2"] == "not_applicable", gene
            # Dead count / female-only flags must be cleared once BS2 is off.
            assert rows[gene]["bs2_count"] == "", gene
            assert rows[gene]["bs2_female_only"] == "", gene

    def test_committed_tsv_keeps_explicit_gnomad_bs2_applicable(self):
        # BMPR2 / PIK3R2 sanction gnomAD homozygote counting directly
        # ("≥3 homozygotes in gnomAD"), so BS2 stays gnomAD-automatable.
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        for gene in ("BMPR2", "PIK3R2"):
            assert rows[gene]["bs2"] == "applicable", gene
        # BMPR2 is dominant with incomplete penetrance: BS2 counts homozygotes
        # only, so healthy heterozygous carriers must not fire it (was the
        # reported false positive). PIK3R2 also sanctions heterozygous family
        # members, so it is NOT hom-only.
        assert rows["BMPR2"]["bs2_hom_only"] == "1"
        assert rows["PIK3R2"]["bs2_hom_only"] == ""

    def test_committed_tsv_reclassifies_gnomad_homozygote_genes(self):
        # Genes whose VCEP BS2 has an explicit gnomAD homozygote-count path
        # ("≥N homozygotes" / "homozygous state in a healthy adult"). These were
        # previously forced not_applicable and produced FALSE NEGATIVES; they are
        # recessive/mode-agnostic so the evaluator counts gnomAD homozygotes.
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        for gene in ("IL7R", "RAG1", "RAG2", "PAH", "POLG", "GAMT", "ETHE1", "ADA"):
            assert rows[gene]["bs2"] == "applicable", gene

    def test_committed_tsv_tiered_count_uses_supporting_threshold(self):
        # Tiered specs (Supporting ">=3 homozygotes", Strong ">=6"): the binary
        # evaluator must fire at the LOWEST (Supporting) bar. GUCY2D nhomalt=3 was
        # the reported false negative when the count resolved to the Strong 6.
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        assert rows["GUCY2D"]["bs2_count"] == "3"
        assert rows["AIPL1"]["bs2_count"] == "3"
        # (PDHA1's ">=16 hemizygotes" min-extraction is covered at the function
        # level in TestBS2Count; the gene itself is forced not_applicable as an
        # X-linked internal-cohort spec, see TestBS2XLinkedInternal.)

    def test_committed_tsv_xlinked_internal_not_applicable(self):
        # X-linked genes whose VCEP BS2 requires a phenotyped male/hemizygote
        # cohort gnomAD cannot supply (functional/lab confirmation or phenotyped
        # unaffected relatives). Forced not_applicable so a gnomAD-count BS2 never
        # falsely fires on a pathogenic X-linked variant.
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
        for gene in ("CDKL5", "RS1", "PDHA1", "RPGR", "IL2RG",
                     "SLC9A6", "F9", "SLC6A8", "MECP2", "F8"):
            assert rows[gene]["bs2"] == "not_applicable", gene
            assert rows[gene]["bs2_count"] == "", gene
