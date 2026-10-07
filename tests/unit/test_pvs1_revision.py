"""Tests for the Genome Medicine revision: PTC-based NMD, ClinGen HI/GDV, null counting."""
from pathlib import Path

import pytest

from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType, Pathogenicity
from acmg_classifier.pvs1.nmd_predictor import (
    nmd_status, ptc_codon, transcript_structure,
)


def _pc(ctype, hgvs_p, exon, tx="NM_000314.8", pos=None):
    return ConsequenceInfo(
        transcript_id=tx, gene_id="G", gene_symbol="PTEN", consequence=ctype,
        biotype="protein_coding", is_mane_select=True, exon=exon,
        hgvs_p=hgvs_p, protein_position=pos,
    )


class TestPTC:
    def test_stop_gained_codon(self):
        assert ptc_codon(_pc(ConsequenceType.STOP_GAINED, "NP_000305.3:p.Arg130Ter", "5/9")) == 130

    def test_frameshift_new_stop(self):
        pc = _pc(ConsequenceType.FRAMESHIFT, "NP_000305.3:p.(Lys267ArgfsTer9)", "7/9")
        assert ptc_codon(pc) == 275

    def test_frameshift_without_stop(self):
        assert ptc_codon(_pc(ConsequenceType.FRAMESHIFT, "p.Lys267ArgfsTer?", "7/9")) is None


class TestNMDStructure:
    def test_pten_structure_loaded(self):
        st = transcript_structure("NM_000314.8")
        assert st is not None and st.n_exons == 9
        # PTEN CDS = 403 aa + stop = 1212 nt
        assert st.cds_end - st.cds_start + 1 == 1212

    def test_early_ptc_nmd(self):
        nmd, note = nmd_status(_pc(ConsequenceType.STOP_GAINED, "p.Arg130Ter", "5/9"))
        assert nmd is True and "NMD predicted" in note

    def test_last_exon_ptc_escapes(self):
        # PTEN codon 380 lies in exon 9 (last).
        nmd, note = nmd_status(_pc(ConsequenceType.STOP_GAINED, "p.Glu380Ter", "9/9"))
        assert nmd is False and "last exon" in note

    def test_frameshift_stop_in_last_exon_escapes(self):
        # Variant in exon 8 (penultimate) but the new stop falls in exon 9.
        st = transcript_structure("NM_000314.8")
        last_exon_first_codon = (st.last_junction - st.cds_start + 1) // 3 + 2
        pc = _pc(ConsequenceType.FRAMESHIFT,
                 f"p.Lys{last_exon_first_codon - 5}ArgfsTer20", "8/9")
        nmd, note = nmd_status(pc)
        assert nmd is False

    def test_penultimate_last_50nt_escapes(self):
        st = transcript_structure("NM_000314.8")
        # codon whose last nt is 30 nt before the last junction
        codon = (st.last_junction - 30 - st.cds_start + 1) // 3
        nmd, note = nmd_status(_pc(ConsequenceType.STOP_GAINED, f"p.Arg{codon}Ter", "8/9"))
        assert nmd is False and "penultimate" in note

    def test_unknown_transcript_falls_back_to_exon_rule(self):
        nmd, note = nmd_status(_pc(ConsequenceType.STOP_GAINED, "p.Arg10Ter", "3/9", tx="NM_999999.1"))
        assert nmd is True and "exon-number rule" in note


class TestNullCount:
    def test_is_null_change(self):
        from acmg_classifier.local_db.clinvar_sqlite import _is_null_change
        assert _is_null_change(None, "NP_1:p.Arg10Ter")
        assert _is_null_change(None, "p.Lys5fs")
        assert _is_null_change("NM_1:c.100+1G>A", None)
        assert _is_null_change("NM_1:c.101-2A>G", None)
        assert not _is_null_change("NM_1:c.100+5G>A", None)
        assert not _is_null_change("NM_1:c.100+12G>A", None)
        assert not _is_null_change(None, "p.Ter100ArgextTer5")
        assert not _is_null_change("c.10A>G", "p.Lys4Arg")


class TestClinGenGeneDB:
    def test_hi_scores(self, tmp_path):
        from acmg_classifier.local_db.clingen_gene_db import hi_score, load_hi_scores
        p = tmp_path / "dosage.tsv"
        p.write_text("#ClinGen Gene Curation Results\n#Genome Build GRCh38\n"
                     "#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\t"
                     "Haploinsufficiency Score\tHaploinsufficiency Description\n"
                     "PTEN\t5728\t10q23.31\tchr10:1-2\t3\tSufficient Evidence\n"
                     "ABC1\t1\t1p\tchr1:1-2\t1\tLittle Evidence\n", encoding="utf-8")
        load_hi_scores.cache_clear()
        assert hi_score(p, "PTEN") == "3"
        assert hi_score(p, "ABC1") == "1"
        assert hi_score(p, "NONE") is None

    def _gdv(self, tmp_path):
        p = tmp_path / "gdv.csv"
        p.write_text('"CLINGEN GENE DISEASE VALIDITY CURATIONS"\n"FILE CREATED: 2026-06-01"\n'
                     '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP",'
                     '"CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n'
                     '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++","+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","+++++++++++++++++"\n'
                     '"GOOD1","HGNC:1","d","MONDO:1","AD","SOP9","Definitive","x","2020","g"\n'
                     '"WEAK1","HGNC:2","d","MONDO:2","AD","SOP9","Limited","x","2020","g"\n'
                     '"WEAK1","HGNC:2","d2","MONDO:3","AR","SOP9","Disputed","x","2020","g"\n'
                     '"MIX1","HGNC:3","d","MONDO:4","AD","SOP9","Limited","x","2020","g"\n'
                     '"MIX1","HGNC:3","d2","MONDO:5","AR","SOP9","Moderate","x","2020","g"\n',
                     encoding="utf-8")
        return p

    def test_validity_status(self, tmp_path):
        from acmg_classifier.local_db.clingen_gene_db import gene_validity_status, load_gene_validity
        p = self._gdv(tmp_path)
        load_gene_validity.cache_clear()
        assert gene_validity_status(p, "GOOD1")[0] == "supported"
        assert gene_validity_status(p, "WEAK1")[0] == "unsupported"
        assert gene_validity_status(p, "MIX1")[0] == "supported"
        assert gene_validity_status(p, "OTHER")[0] == "not_curated"

    def test_cap_at_vus(self, tmp_path):
        from types import SimpleNamespace
        from acmg_classifier.classification.gene_validity import apply_gene_validity
        from acmg_classifier.local_db.clingen_gene_db import load_gene_validity
        load_gene_validity.cache_clear()
        cfg = SimpleNamespace(clingen_gene_validity_csv=self._gdv(tmp_path))
        a, b, w = apply_gene_validity("WEAK1", Pathogenicity.PATHOGENIC,
                                      Pathogenicity.LIKELY_PATHOGENIC, cfg)
        assert a == b == Pathogenicity.VUS and w and "capped at VUS" in w[0]
        a, b, w = apply_gene_validity("GOOD1", Pathogenicity.PATHOGENIC,
                                      Pathogenicity.PATHOGENIC, cfg)
        assert a == Pathogenicity.PATHOGENIC and not w
        a, b, w = apply_gene_validity("WEAK1", Pathogenicity.LIKELY_BENIGN,
                                      Pathogenicity.BENIGN, cfg)
        assert b == Pathogenicity.BENIGN  # benign side untouched
        _, _, w = apply_gene_validity("OTHER", Pathogenicity.PATHOGENIC,
                                      Pathogenicity.PATHOGENIC, cfg)
        assert w and w[0].startswith("GENE_VALIDITY_NOT_CURATED")
