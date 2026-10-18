"""REVEL per-transcript row selection (PP3 discordance audit, Genome Medicine
revision): REVEL v1.3 has one row per transcript-specific amino-acid change, and
the MANE row must be used (IDUA c.530T>G: MANE F>C 0.961, not S>A 0.400; FBN1
c.7016G>A: MANE C>Y 0.973, not A>T 0.209)."""
from acmg_classifier.annotation.orchestrator import _revel_keys
from acmg_classifier.local_db.revel_db import _match_revel
from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType

IDUA = [
    "4\t1001504\tT\tG\t0.400\tS\tA\tENST00000504568",
    "4\t1001504\tT\tG\t0.961\tF\tC\tENST00000247933;ENST00000453894;ENST00000514224",
]
FBN1 = [
    "15\t48427755\tC\tT\t0.209\tA\tT\tENST00000389087",
    "15\t48427755\tC\tT\t0.973\tC\tY\tENST00000316623",
]


def test_transcript_match_wins():
    r = _match_revel(IDUA, 1001504, "T", "G", aa=("S", "A"), transcripts={"ENST00000514224"})
    assert r.score == 0.961 and r.match == "transcript"


def test_amino_acid_match():
    r = _match_revel(FBN1, 48427755, "C", "T", aa=("C", "Y"), transcripts={"ENST00000999999"})
    assert r.score == 0.973 and r.match == "amino_acid"


def test_fallback_max():
    r = _match_revel(FBN1, 48427755, "C", "T")
    assert r.score == 0.973 and r.match == "max"


def test_order_independent():
    a = _match_revel(IDUA, 1001504, "T", "G", transcripts={"ENST00000514224.2"})
    b = _match_revel(list(reversed(IDUA)), 1001504, "T", "G", transcripts={"ENST00000514224"})
    assert a == b


def test_legacy_five_column():
    r = _match_revel(["4\t1001504\tT\tG\t0.400", "4\t1001504\tT\tG\t0.961"], 1001504, "T", "G")
    assert r.score == 0.400 and r.match == "legacy"


def test_no_match():
    assert _match_revel(IDUA, 1001504, "T", "C") is None


def _c(tid, mane=False, cons=ConsequenceType.MISSENSE, aa="F/C"):
    return ConsequenceInfo(transcript_id=tid, gene_id="ENSG", gene_symbol="IDUA",
                           consequence=cons, biotype="protein_coding",
                           is_mane_select=mane, amino_acids=aa)


def test_revel_keys():
    p = _c("NM_000203.5", mane=True)
    aa, tx = _revel_keys(p, [p, _c("ENST00000514224.2", mane=True), _c("ENST00000504568.1", aa="S/A")],
                         {"IDUA": ("NM_000203", "ENST00000514224")})
    assert aa == ("F", "C") and tx == {"ENST00000514224"}


def test_revel_keys_non_missense():
    p = _c("NM_000203.5", mane=True, cons=ConsequenceType.SYNONYMOUS, aa="F")
    aa, tx = _revel_keys(p, [p], {})
    assert aa is None and tx is None
