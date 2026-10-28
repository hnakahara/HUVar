"""Somatic oncogenicity criteria (ClinGen/CGC/VICC SOP, Horak et al. 2022).

Automated criteria and their data (no licensed resources):

  OVS1   null variant in a bona fide tumour suppressor (NCG canonical TSG-only
         role, curated overrides; NMD / 3'-end check)
  OS1    same amino-acid change as a ClinVar P/LP variant (>= 2 stars), different
         nucleotide change
  OS2    calibrated functional assay (MaveDB / ENIGMA Table 9 / TP53 VCEP data)
         showing loss of function at Strong or above, tumour suppressors only
  OS3    cancerhotspots.org residue with >= 50 samples and >= 10 with this change
  OM1    critical region without benign variation (HUVar PM1: VCEP hotspot /
         domain tables or the ClinVar P/LP cluster)
  OM2    in-frame indel in an oncogene / TSG, or stop-loss in a TSG
  OM3    cancerhotspots.org residue with < 50 samples and >= 10 with this change
  OM4    different missense at a residue with a ClinVar P/LP missense (>= 2 stars),
         Grantham distance at least ~ that of the comparator
  OP1    in-silico evidence for an effect (HUVar PP3)
  OP3    cancerhotspots.org residue / indel region, this change in < 10 samples
  OP4    absent from gnomAD or grpmax AF <= somatic_op4_max_af
  SBVS1  gnomAD grpmax AF > 5 %
  SBS1   gnomAD grpmax AF > 1 %
  SBS2   calibrated functional assay showing normal function at Strong or above
  SBP1   in-silico evidence for no effect (HUVar BP4)
  SBP2   synonymous, no predicted splicing effect, not highly conserved (HUVar BP7)

OP2 (malignancy with a single genetic aetiology) needs the tumour type and is
applied only from manual evidence; any criterion can be set or overridden that way.
"""
from __future__ import annotations

from typing import Optional

from acmg_classifier.config import Config
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.somatic.models import (
    MANUAL_ONLY, ManualOncEvidence, OncCriterion, OncCriterionResult,
)
from acmg_classifier.somatic.resources import (
    GeneRole, HotspotSite, indel_site, load_gene_roles, load_hotspots, resolve,
    residue_site, short_protein_change,
)

C = OncCriterion
CT = ConsequenceType
_NULL = {CT.STOP_GAINED, CT.FRAMESHIFT, CT.SPLICE_ACCEPTOR, CT.SPLICE_DONOR, CT.START_LOST,
         CT.TRANSCRIPT_ABLATION}
_INFRAME = {CT.INFRAME_DELETION, CT.INFRAME_INSERTION}
_STRONG_UP = {CriterionStrength.STRONG, CriterionStrength.VERY_STRONG}
# Functional evidence sources that come from calibrated assay data (not text mining).
_FUNCTIONAL_SOURCES = ("[MaveDB]", "[ENIGMA Table 9]", "[TP53 VCEP functional data]")
# A truncation in the last 10 % of the protein that escapes NMD is not counted
# as a null allele (SOP: caution at the extreme 3' end).
_NMD_ESCAPE_MAX_FRACTION = 0.9
# Splice disruption makes amino-acid level comparisons (OS1 / OM4) unreliable.
_SPLICE_CAUTION = 0.5


def _fmt_af(v: Optional[float]) -> str:
    return "n/a" if v is None else f"{v:.3g}"


class SomaticEvaluator:
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._roles = load_gene_roles(resolve(cfg.somatic_gene_roles_tsv, "somatic_gene_roles.tsv"))
        self._hotspots = load_hotspots(resolve(cfg.cancer_hotspots_tsv, "cancer_hotspots_v2.tsv"))

    # ------------------------------------------------------------------ helpers
    def gene_role(self, gene: Optional[str]) -> Optional[GeneRole]:
        return self._roles.get(gene or "")

    def hotspot(self, ann: AnnotationData) -> tuple[Optional[HotspotSite], Optional[str], str]:
        """(site, change key, note) for a missense or in-frame indel."""
        pc = ann.primary_consequence
        if pc is None:
            return None, None, ""
        gene = pc.gene_symbol
        if pc.consequence == CT.MISSENSE:
            site = residue_site(self._hotspots, gene, pc.protein_position)
            if site is None:
                return None, None, ""
            ref, _, alt = (pc.amino_acids or "/").partition("/")
            if ref and site.ref_aa and ref != site.ref_aa:
                return None, None, (f"cancerhotspots {gene} residue {pc.protein_position} is "
                                    f"{site.ref_aa}, the variant reference is {ref} (different "
                                    "transcript numbering) - hotspot not used")
            return site, alt, ""
        if pc.consequence in _INFRAME:
            change = short_protein_change(pc.hgvs_p)
            if not change:
                return None, None, ""
            import re
            nums = [int(n) for n in re.findall(r"[A-Z*](\d+)", change)]
            site = indel_site(self._hotspots, gene, min(nums, default=None),
                              max(nums, default=None))
            return site, change, ""
        return None, None, ""

    # --------------------------------------------------------------- evaluation
    def evaluate(
        self,
        variant: VariantRecord,
        ann: AnnotationData,
        germline: list[CriteriaResult],
        manual: Optional[list[ManualOncEvidence]] = None,
    ) -> tuple[list[OncCriterionResult], dict]:
        """Evaluate all criteria. ``germline`` is the HUVar ACMG/AMP evaluation of
        the same variant (PP3/BP4/BP7/PM1/PS3/BS3 are reused). Returns the results
        (one per criterion, SOP order) and context (gene role, hotspot text)."""
        pc = ann.primary_consequence
        gene = pc.gene_symbol if pc else None
        role = self.gene_role(gene)
        g = {r.criterion: r for r in germline}
        res: dict[OncCriterion, OncCriterionResult] = {c: OncCriterionResult(criterion=c) for c in C}

        def met(c: OncCriterion, ev: str) -> None:
            res[c].met, res[c].evidence = True, ev

        def note(c: OncCriterion, ev: str) -> None:
            res[c].evidence = ev

        csq = pc.consequence if pc else None
        role_txt = role.describe() if role else f"{gene or 'gene'}: not a canonical cancer gene (NCG)"

        # OVS1 -------------------------------------------------------------
        if csq in _NULL:
            if role is None or not role.is_tsg:
                note(C.OVS1, f"Null variant ({csq.value}) but {role_txt} - OVS1 needs a "
                             "bona fide tumour suppressor")
            else:
                ok, why = self._null_ok(pc)
                (met if ok else note)(C.OVS1, f"{csq.value} in tumour suppressor {role_txt}; {why}")
        else:
            note(C.OVS1, "Not a null variant")

        # OM2 --------------------------------------------------------------
        if csq in _INFRAME:
            if role is not None and role.is_cancer_gene:
                met(C.OM2, f"In-frame {csq.value.split('_')[1]} changes protein length; {role_txt}")
            else:
                note(C.OM2, f"In-frame indel but {role_txt}")
        elif csq == CT.STOP_LOST:
            if role is not None and role.role in ("TSG", "dual"):
                met(C.OM2, f"Stop-loss extends the protein; {role_txt}")
            else:
                note(C.OM2, f"Stop-loss but {role_txt} (OM2 needs a tumour suppressor)")
        else:
            note(C.OM2, "Not an in-frame indel or stop-loss")

        # OS1 / OM4 (ClinVar comparators) ------------------------------------
        self._clinvar_comparators(variant, ann, res)

        # Hotspots OS3 / OM3 / OP3 ------------------------------------------
        site, change, hs_note = self.hotspot(ann)
        hotspot_txt = ""
        if site is not None:
            n_pos, n_var = site.position_count, site.variants.get(change, 0)
            hotspot_txt = site.describe()
            ev = f"{hotspot_txt}; this change ({change}) in {n_var} samples"
            if n_var >= 10 and n_pos >= 50:
                met(C.OS3, ev)
            elif n_var >= 10:
                met(C.OM3, ev)
            else:
                met(C.OP3, ev)
        else:
            for c in (C.OS3, C.OM3, C.OP3):
                note(c, hs_note or "Not at a cancerhotspots.org hotspot")

        # OM1 (HUVar PM1) ----------------------------------------------------
        pm1 = g.get(ACMGCriterion.PM1)
        if pm1 is not None and pm1.triggered and not pm1.suppressed:
            met(C.OM1, f"Critical region (HUVar PM1): {pm1.evidence}")
        else:
            note(C.OM1, "Not in a critical region without benign variation (HUVar PM1 not met)")

        # OS2 / SBS2 (calibrated functional data) ----------------------------
        self._functional(g, role, role_txt, res)

        # OP1 / SBP1 / SBP2 (in silico) --------------------------------------
        pp3, bp4, bp7 = g.get(ACMGCriterion.PP3), g.get(ACMGCriterion.BP4), g.get(ACMGCriterion.BP7)
        if pp3 is not None and pp3.triggered:
            met(C.OP1, f"In-silico evidence for an effect (HUVar PP3 {pp3.strength.value}): {pp3.evidence}")
        else:
            note(C.OP1, (pp3.evidence if pp3 else "") or "HUVar PP3 not met")
        if bp4 is not None and bp4.triggered:
            met(C.SBP1, f"In-silico evidence for no effect (HUVar BP4 {bp4.strength.value}): {bp4.evidence}")
        else:
            note(C.SBP1, (bp4.evidence if bp4 else "") or "HUVar BP4 not met")
        if csq == CT.SYNONYMOUS and bp7 is not None and bp7.triggered:
            met(C.SBP2, f"Synonymous, no predicted splice effect, not highly conserved (HUVar BP7): {bp7.evidence}")
        else:
            note(C.SBP2, "Not synonymous" if csq != CT.SYNONYMOUS else
                 ((bp7.evidence if bp7 else "") or "HUVar BP7 not met"))

        # Population frequency OP4 / SBS1 / SBVS1 ---------------------------
        self._frequency(ann, res)

        # OP2 / manual ---------------------------------------------------------
        note(C.OP2, "Needs the tumour type (malignancy with a single genetic aetiology) - "
                    "apply manually")
        self._apply_manual(res, manual)
        self._combination_rules(res)
        return [res[c] for c in C], {"role": role, "role_text": role_txt, "hotspot": hotspot_txt}

    # ------------------------------------------------------------- components
    @staticmethod
    def _null_ok(pc) -> tuple[bool, str]:
        if pc.consequence not in (CT.STOP_GAINED, CT.FRAMESHIFT):
            return True, "canonical splice site / start loss / ablation"
        from acmg_classifier.pvs1.nmd_predictor import nmd_status, ptc_codon, transcript_structure
        nmd, why = nmd_status(pc)
        if nmd:
            return True, why
        st = transcript_structure(pc.transcript_id)
        codon = ptc_codon(pc) or pc.protein_position
        if st is not None and codon:
            aa_len = (st.cds_end - st.cds_start + 1) // 3 - 1
            if aa_len > 0:
                frac = codon / aa_len
                if frac <= _NMD_ESCAPE_MAX_FRACTION:
                    return True, (f"{why}; truncates at residue {codon}/{aa_len} "
                                  f"({100 * (1 - frac):.0f}% of the protein removed)")
                return False, (f"{why}; truncation in the last {100 * (1 - frac):.0f}% of the "
                               f"protein (residue {codon}/{aa_len}) - review manually")
        return False, f"{why}; protein length unknown - review manually"

    def _clinvar_comparators(self, variant: VariantRecord, ann: AnnotationData, res) -> None:
        pc = ann.rule_consequence or ann.primary_consequence
        if pc is None or pc.consequence != CT.MISSENSE:
            res[C.OS1].evidence = "Not a missense variant"
            res[C.OM4].evidence = "Not a missense variant"
            return
        sp = ann.splice
        if sp is not None and sp.is_available and (sp.max_delta or 0) >= _SPLICE_CAUTION:
            msg = (f"Predicted to affect splicing ({sp.tool} {sp.max_delta:.2f}) - amino-acid "
                   "comparison not used")
            res[C.OS1].evidence = res[C.OM4].evidence = msg
            return
        db = self._cfg.clinvar_sqlite
        stars = int(getattr(self._cfg, "somatic_clinvar_min_stars", 2))
        from acmg_classifier.local_db.clinvar_sqlite import (
            query_same_aa_change, query_same_codon_different_aa,
        )
        same = query_same_aa_change(db, pc.gene_symbol, pc.hgvs_p, variant.chrom, variant.pos,
                                    variant.ref, variant.alt, min_stars=stars)
        if same:
            r = same[0]
            res[C.OS1].met = True
            res[C.OS1].evidence = (f"Same amino-acid change as ClinVar {r.clinical_significance} "
                                   f"{r.hgvs_c or ''} (VCV {r.variation_id}, {r.star_rating} stars)")
        else:
            res[C.OS1].evidence = (f"No ClinVar P/LP (>= {stars} stars) variant with the same "
                                   "amino-acid change")
        others = query_same_codon_different_aa(db, pc.gene_symbol, pc.protein_position, pc.hgvs_p,
                                               min_stars=stars, query_chrom=variant.chrom,
                                               query_pos=variant.pos)
        if not others:
            res[C.OM4].evidence = (f"No ClinVar P/LP (>= {stars} stars) missense at residue "
                                   f"{pc.protein_position}")
            return
        from acmg_classifier.criteria.grantham import grantham_distance
        import re
        ref, _, alt = (pc.amino_acids or "/").partition("/")
        cand = grantham_distance(ref, alt)
        best, seen = None, []
        for r in others:
            m = re.search(r"p\.\(?([A-Z][a-z]{2})\d+([A-Z][a-z]{2})", r.hgvs_p or "")
            d = grantham_distance(m.group(1), m.group(2)) if m else None
            seen.append(d)
            if cand is None or d is None or cand >= 0.9 * d:
                best = (r, d)
                break
        if best is None:
            res[C.OM4].evidence = (f"ClinVar P/LP missense at residue {pc.protein_position} have a "
                                   f"larger amino-acid difference (Grantham {cand} vs "
                                   f"{', '.join(str(x) for x in seen[:5])})")
            return
        r, d = best
        res[C.OM4].met = True
        res[C.OM4].evidence = (f"Different missense at residue {pc.protein_position}: ClinVar "
                               f"{r.clinical_significance} {(r.hgvs_p or '').split(':')[-1]} "
                               f"(VCV {r.variation_id}, {r.star_rating} stars); Grantham {cand} "
                               f"vs {d}")

    @staticmethod
    def _functional(g, role, role_txt, res) -> None:
        for crit, onc, label in ((ACMGCriterion.PS3, C.OS2, "abnormal (loss of) function"),
                                 (ACMGCriterion.BS3, C.SBS2, "normal function")):
            r = g.get(crit)
            if r is None or not (r.evidence or "").startswith(_FUNCTIONAL_SOURCES):
                res[onc].evidence = "No calibrated functional assay data (MaveDB / ENIGMA / TP53 VCEP)"
                continue
            if not r.triggered:
                res[onc].evidence = r.evidence
                continue
            if role is None or not role.is_tsg:
                res[onc].evidence = (f"{r.evidence} - not used: the assays measure loss of function, "
                                     f"which supports oncogenicity only in a tumour suppressor "
                                     f"({role_txt})")
                continue
            if r.strength not in _STRONG_UP:
                res[onc].evidence = (f"{r.evidence} - {crit.value}_{r.strength.value} only; "
                                     f"{onc.value} needs well-established (Strong) evidence")
                continue
            res[onc].met = True
            res[onc].evidence = f"Calibrated assay, {label} ({crit.value}_{r.strength.value}): {r.evidence}"

    def _frequency(self, ann: AnnotationData, res) -> None:
        gd = ann.gnomad
        if gd is None:
            for c in (C.OP4, C.SBS1, C.SBVS1):
                res[c].evidence = "gnomAD data unavailable"
            return
        af = gd.popmax_af if gd.popmax_af is not None else gd.af
        absent = (gd.ac or 0) == 0 and not af
        cut = float(getattr(self._cfg, "somatic_op4_max_af", 1e-5))
        detail = (f"gnomAD grpmax AF {_fmt_af(gd.popmax_af)} ({gd.popmax_pop or '-'}), "
                  f"AF {_fmt_af(gd.af)}, AC {gd.ac}")
        if af is not None and af > 0.05:
            res[C.SBVS1].met, res[C.SBVS1].evidence = True, f"{detail} > 5%"
        elif af is not None and af > 0.01:
            res[C.SBS1].met, res[C.SBS1].evidence = True, f"{detail} > 1%"
        else:
            res[C.SBVS1].evidence = res[C.SBS1].evidence = f"{detail} <= 1%"
        if absent:
            res[C.OP4].met, res[C.OP4].evidence = True, "Absent from gnomAD"
        elif af is not None and af <= cut:
            res[C.OP4].met, res[C.OP4].evidence = True, f"{detail} <= {cut:g} (extremely low)"
        else:
            res[C.OP4].evidence = f"{detail} > {cut:g}"
        if not gd.filter_pass:
            for c in (C.OP4, C.SBS1, C.SBVS1):
                res[c].evidence += " (gnomAD QC filter not passed)"

    @staticmethod
    def _apply_manual(res, manual: Optional[list[ManualOncEvidence]]) -> None:
        for e in manual or []:
            r = res[e.criterion]
            r.met, r.manual = e.met, True
            r.evidence = f"[Manual] {e.evidence}".strip() if e.met else f"[Manual] not met. {e.evidence}".strip()
        for c in MANUAL_ONLY:
            if not res[c].manual:
                res[c].met = False

    @staticmethod
    def _combination_rules(res) -> None:
        """SOP Table 2 exclusions, applied as suppression (kept in the audit trail)."""
        def active(c):
            return res[c].met and not res[c].suppressed

        def suppress(c, why):
            if res[c].met and not res[c].suppressed and not res[c].manual:
                res[c].suppressed = True
                res[c].evidence += f" [not counted: {why}]"

        # Hotspot counts are amino-acid level: with OS1 (same amino-acid change
        # already established) they would count the same observation twice.
        if active(C.OS1):
            for c in (C.OS3, C.OM3, C.OP3):
                suppress(c, "OS1 met (same amino-acid change)")
        if active(C.OS1) or active(C.OS3):
            suppress(C.OM1, "OS1/OS3 met")
        if active(C.OS1) or active(C.OS3) or active(C.OM1):
            suppress(C.OM4, "OS1/OS3/OM1 met")
        if active(C.OM1) or active(C.OM4):
            suppress(C.OM3, "OM1/OM4 met")
        if active(C.OVS1):
            suppress(C.OM2, "OVS1 met")
            suppress(C.OP1, "OVS1 met (null variant)")
        # Frequency criteria are exclusive: SBVS1 > SBS1 > OP4.
        if active(C.SBVS1):
            suppress(C.SBS1, "SBVS1 met")
            suppress(C.OP4, "SBVS1 met")
        if active(C.SBS1):
            suppress(C.OP4, "SBS1 met")
        # Opposite in-silico calls cannot both count.
        if active(C.OP1) and active(C.SBP1):
            suppress(C.SBP1, "OP1 met")
