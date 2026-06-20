"""Gene-specific VCEP PVS1 trees added for GATM (GN025) and PIK3R1 (GN160).

GATM mirrors the sibling GAMT NMD-based rule (NMD → PVS1; escape → 10% rule;
initiation codon → Moderate). PIK3R1 uses codon bands: NMD region (codons
306-630) → PVS1; 631-718 → Strong; 719-724 → Moderate; the N-terminal
agammaglobulinaemia-only region (codons 1-305) → N/A; initiation codon → N/A.
"""
from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType, CriterionStrength
from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1

S = CriterionStrength


def _pc(gene, consequence, protein_position=None, exon=None, intron=None):
    return ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=consequence, biotype="protein_coding",
        protein_position=protein_position, exon=exon, intron=intron,
    )


class TestGATM:
    def test_nmd_predicted_very_strong(self):
        # Mid-transcript PTC → NMD predicted → PVS1.
        assert evaluate_vcep_pvs1(
            _pc("GATM", ConsequenceType.FRAMESHIFT, 100, exon="3/9"))[0] == S.VERY_STRONG

    def test_nmd_escape_small_loss_moderate(self):
        # Last-exon PTC near the C-terminus removes <10% (423 aa) → Moderate.
        assert evaluate_vcep_pvs1(
            _pc("GATM", ConsequenceType.STOP_GAINED, 420, exon="9/9"))[0] == S.MODERATE

    def test_start_lost_moderate(self):
        assert evaluate_vcep_pvs1(_pc("GATM", ConsequenceType.START_LOST))[0] == S.MODERATE

    def test_splice_very_strong(self):
        assert evaluate_vcep_pvs1(
            _pc("GATM", ConsequenceType.SPLICE_DONOR, intron="3/8"))[0] == S.VERY_STRONG

    def test_deletion_very_strong(self):
        assert evaluate_vcep_pvs1(
            _pc("GATM", ConsequenceType.TRANSCRIPT_ABLATION))[0] == S.VERY_STRONG


class TestPIK3R1:
    def test_nmd_band_very_strong(self):
        # Codon 400 (in c.917-1890 NMD region, codons 306-630) → PVS1.
        assert evaluate_vcep_pvs1(
            _pc("PIK3R1", ConsequenceType.STOP_GAINED, 400))[0] == S.VERY_STRONG

    def test_csh2_band_strong(self):
        # Codon 700 (cSH2 631-718) → Strong.
        assert evaluate_vcep_pvs1(
            _pc("PIK3R1", ConsequenceType.FRAMESHIFT, 700))[0] == S.STRONG

    def test_cterminal_band_moderate(self):
        assert evaluate_vcep_pvs1(
            _pc("PIK3R1", ConsequenceType.STOP_GAINED, 720))[0] == S.MODERATE

    def test_nterminal_region_na(self):
        # Codon 100 (c.4-916 agammaglobulinaemia-only region) → outside bands → N/A.
        s, ev = evaluate_vcep_pvs1(_pc("PIK3R1", ConsequenceType.STOP_GAINED, 100))
        assert s == S.NOT_MET and "outside" in ev

    def test_start_lost_na(self):
        assert evaluate_vcep_pvs1(_pc("PIK3R1", ConsequenceType.START_LOST))[0] == S.NOT_MET

    def test_splice_very_strong(self):
        assert evaluate_vcep_pvs1(
            _pc("PIK3R1", ConsequenceType.SPLICE_ACCEPTOR, intron="9/15"))[0] == S.VERY_STRONG

    def test_deletion_very_strong(self):
        assert evaluate_vcep_pvs1(
            _pc("PIK3R1", ConsequenceType.TRANSCRIPT_ABLATION))[0] == S.VERY_STRONG


class TestOTOFcarveout:
    # OTOF exon 46 = codons 1947-1997: PVS1-exception region (withheld).
    def test_truncation_in_exon46_withheld(self):
        s, ev = evaluate_vcep_pvs1(_pc("OTOF", ConsequenceType.STOP_GAINED, 1960))
        assert s == S.NOT_MET and "exception region" in ev

    def test_truncation_outside_defers(self):
        # Codon 1000 is outside the carve-out → no per-gene rule → defer (None).
        assert evaluate_vcep_pvs1(_pc("OTOF", ConsequenceType.STOP_GAINED, 1000)) is None

    def test_splice_and_start_defer(self):
        assert evaluate_vcep_pvs1(_pc("OTOF", ConsequenceType.SPLICE_DONOR, intron="5/46")) is None
        assert evaluate_vcep_pvs1(_pc("OTOF", ConsequenceType.START_LOST)) is None
        assert evaluate_vcep_pvs1(_pc("OTOF", ConsequenceType.TRANSCRIPT_ABLATION)) is None


class TestMYO15Acarveout:
    # Exon 8 = codons 1345-1346; exon 26 = codons 1971-1988.
    def test_exon8_withheld(self):
        assert evaluate_vcep_pvs1(_pc("MYO15A", ConsequenceType.FRAMESHIFT, 1345))[0] == S.NOT_MET

    def test_exon26_withheld(self):
        assert evaluate_vcep_pvs1(_pc("MYO15A", ConsequenceType.STOP_GAINED, 1980))[0] == S.NOT_MET

    def test_between_carveouts_defers(self):
        assert evaluate_vcep_pvs1(_pc("MYO15A", ConsequenceType.STOP_GAINED, 1500)) is None

    def test_early_truncation_defers(self):
        assert evaluate_vcep_pvs1(_pc("MYO15A", ConsequenceType.FRAMESHIFT, 500)) is None


def test_specs_registered():
    from acmg_classifier.pvs1.vcep_pvs1 import _SPECS
    assert _SPECS["GATM"].transcript == "NM_001482.3"
    assert _SPECS["GATM"].aa_len == 423
    assert _SPECS["PIK3R1"].transcript == "NM_181523.3"
    assert _SPECS["PIK3R1"].aa_len == 724
    assert _SPECS["OTOF"].trunc_exception_bands == ((1947, 1997),)
    assert _SPECS["MYO15A"].trunc_exception_bands == ((1345, 1346), (1971, 1988))
