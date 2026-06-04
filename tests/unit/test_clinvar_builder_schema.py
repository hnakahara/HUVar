"""ClinVar XML schema parsing — new RCV_release vs legacy RCV_xml_old_format.

ClinVar's current ``RCV_release/ClinVarRCVRelease`` (ClinVar_RCV_2.3.xsd) moved
the aggregate germline call out of the flat ``<ClinicalSignificance>`` element
into ``<Classifications>/<GermlineClassification>`` (RCV level) and into a
``<Classification>`` with a ``<GermlineClassification>`` text child (SCV level).
The builder must read both so the stale legacy file and the current file both
build correctly.
"""
import xml.etree.ElementTree as ET

from acmg_classifier.setup.clinvar_builder import (
    _gene_symbol, _parse_clinvarset, _rcv_classification, _scv_significance,
    _star_rating,
)

# New RCV_release format (ClinVar_RCV_2.3): expert-panel-style record with a
# coding+protein HGVS and an affected, P/LP SCV.
_NEW = """
<ClinVarSet ID="1">
  <ReferenceClinVarAssertion ID="9001">
    <ClinVarAccession Acc="RCV000000001" Type="RCV"/>
    <Classifications>
      <GermlineClassification>
        <ReviewStatus>reviewed by expert panel</ReviewStatus>
        <Description DateLastEvaluated="2026-01-17" SubmissionCount="2">Likely pathogenic</Description>
      </GermlineClassification>
    </Classifications>
    <MeasureSet Type="Variant" ID="18397">
      <Measure Type="single nucleotide variant" ID="33436">
        <AttributeSet>
          <Attribute Type="HGVS, coding, RefSeq">NM_005026.5:c.1570T>A</Attribute>
        </AttributeSet>
        <AttributeSet>
          <Attribute Type="HGVS, protein, RefSeq">NP_005017.3:p.Tyr524Asn</Attribute>
        </AttributeSet>
        <SequenceLocation Assembly="GRCh38" Chr="1" positionVCF="9720791"
                          referenceAlleleVCF="A" alternateAlleleVCF="T"/>
        <MeasureRelationship>
          <Symbol><ElementValue Type="Preferred">PIK3CD</ElementValue></Symbol>
        </MeasureRelationship>
      </Measure>
    </MeasureSet>
  </ReferenceClinVarAssertion>
  <ClinVarAssertion>
    <Classification DateLastEvaluated="2026-01-17">
      <ReviewStatus>criteria provided, single submitter</ReviewStatus>
      <GermlineClassification>Likely pathogenic</GermlineClassification>
    </Classification>
    <ObservedIn>
      <Sample><AffectedStatus>yes</AffectedStatus></Sample>
    </ObservedIn>
  </ClinVarAssertion>
</ClinVarSet>
"""

# Legacy RCV_xml_old_format equivalent (flat ClinicalSignificance).
_OLD = """
<ClinVarSet ID="1">
  <ReferenceClinVarAssertion ID="9002">
    <ClinVarAccession Acc="RCV000000002" Type="RCV"/>
    <ClinicalSignificance DateLastEvaluated="2025-07-01">
      <ReviewStatus>criteria provided, multiple submitters, no conflicts</ReviewStatus>
      <Description>Pathogenic</Description>
    </ClinicalSignificance>
    <MeasureSet Type="Variant" ID="18398">
      <Measure Type="single nucleotide variant" ID="33437">
        <AttributeSet>
          <Attribute Type="HGVS, protein, RefSeq">NP_005017.3:p.Trp275Cys</Attribute>
        </AttributeSet>
        <SequenceLocation Assembly="GRCh38" Chr="1" positionVCF="40819463"
                          referenceAlleleVCF="G" alternateAlleleVCF="C"/>
        <MeasureRelationship>
          <Symbol><ElementValue Type="Preferred">KCNQ4</ElementValue></Symbol>
        </MeasureRelationship>
      </Measure>
    </MeasureSet>
  </ReferenceClinVarAssertion>
  <ClinVarAssertion>
    <ClinicalSignificance>
      <Description>Pathogenic</Description>
    </ClinicalSignificance>
    <ObservedIn>
      <Sample><AffectedStatus>yes</AffectedStatus></Sample>
    </ObservedIn>
  </ClinVarAssertion>
</ClinVarSet>
"""


def _field(row, name):
    cols = ("variation_id", "chrom", "pos", "ref", "alt", "gene_symbol",
            "hgvs_c", "hgvs_p", "amino_acid_change", "codon_position",
            "clinical_significance", "review_status", "star_rating",
            "last_evaluated", "affected_cases", "functional_evidence",
            "segregation_evidence")
    return row[cols.index(name)]


class TestRcvClassification:
    def test_new_schema(self):
        ra = ET.fromstring(_NEW).find(".//ReferenceClinVarAssertion")
        sig, rev, dle = _rcv_classification(ra)
        assert sig == "Likely pathogenic"
        assert rev == "reviewed by expert panel"
        assert dle == "2026-01-17"
        assert _star_rating(rev) == 3

    def test_legacy_schema(self):
        ra = ET.fromstring(_OLD).find(".//ReferenceClinVarAssertion")
        sig, rev, dle = _rcv_classification(ra)
        assert sig == "Pathogenic"
        assert rev == "criteria provided, multiple submitters, no conflicts"
        assert dle == "2025-07-01"
        assert _star_rating(rev) == 2


_OVERLAP_GENE = """
<ClinVarSet ID="3">
  <ReferenceClinVarAssertion ID="9003">
    <MeasureSet Type="Variant" ID="1">
      <Measure Type="single nucleotide variant" ID="1">
        <MeasureRelationship Type="within multiple genes by overlap">
          <Symbol><ElementValue Type="Preferred">LOC126861615</ElementValue></Symbol>
        </MeasureRelationship>
        <MeasureRelationship Type="variant in gene">
          <Symbol><ElementValue Type="Preferred">PAH</ElementValue></Symbol>
        </MeasureRelationship>
      </Measure>
    </MeasureSet>
  </ReferenceClinVarAssertion>
</ClinVarSet>
"""

# Both relationships are 'overlap' — fall back to the non-LOC symbol.
_OVERLAP_BOTH = """
<ClinVarSet ID="4">
  <ReferenceClinVarAssertion ID="9004">
    <MeasureSet Type="Variant" ID="1">
      <Measure Type="single nucleotide variant" ID="1">
        <MeasureRelationship Type="within multiple genes by overlap">
          <Symbol><ElementValue Type="Preferred">LOC100</ElementValue></Symbol>
        </MeasureRelationship>
        <MeasureRelationship Type="within multiple genes by overlap">
          <Symbol><ElementValue Type="Preferred">BRCA1</ElementValue></Symbol>
        </MeasureRelationship>
      </Measure>
    </MeasureSet>
  </ReferenceClinVarAssertion>
</ClinVarSet>
"""


class TestGeneSymbol:
    def test_variant_in_gene_beats_overlap_locus(self):
        ra = ET.fromstring(_OVERLAP_GENE).find(".//ReferenceClinVarAssertion")
        assert _gene_symbol(ra) == "PAH"   # not LOC126861615 (the first listed)

    def test_non_loc_preferred_on_type_tie(self):
        ra = ET.fromstring(_OVERLAP_BOTH).find(".//ReferenceClinVarAssertion")
        assert _gene_symbol(ra) == "BRCA1"

    def test_single_gene(self):
        ra = ET.fromstring(_NEW).find(".//ReferenceClinVarAssertion")
        assert _gene_symbol(ra) == "PIK3CD"


class TestScvSignificance:
    def test_new_schema(self):
        scv = ET.fromstring(_NEW).find(".//ClinVarAssertion")
        assert _scv_significance(scv) == "likely pathogenic"

    def test_legacy_schema(self):
        scv = ET.fromstring(_OLD).find(".//ClinVarAssertion")
        assert _scv_significance(scv) == "pathogenic"


class TestParseClinvarsetEndToEnd:
    def test_new_format_row(self):
        row = _parse_clinvarset(ET.fromstring(_NEW), assembly="GRCh38")
        assert row is not None
        assert _field(row, "gene_symbol") == "PIK3CD"
        assert _field(row, "hgvs_p") == "NP_005017.3:p.Tyr524Asn"
        assert _field(row, "codon_position") == 524
        assert _field(row, "clinical_significance") == "Likely pathogenic"
        assert _field(row, "star_rating") == 3            # expert panel
        assert _field(row, "last_evaluated") == "2026-01-17"
        assert _field(row, "affected_cases") == 1         # P/LP SCV, affected yes
        assert _field(row, "chrom") == "1"
        assert _field(row, "pos") == 9720791

    def test_legacy_format_row(self):
        row = _parse_clinvarset(ET.fromstring(_OLD), assembly="GRCh38")
        assert row is not None
        assert _field(row, "gene_symbol") == "KCNQ4"
        assert _field(row, "codon_position") == 275
        assert _field(row, "clinical_significance") == "Pathogenic"
        assert _field(row, "star_rating") == 2
        assert _field(row, "affected_cases") == 1
