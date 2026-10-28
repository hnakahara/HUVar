"""Somatic oncogenicity criteria (ClinGen/CGC/VICC SOP, Horak et al. 2022).

Automated criteria and their data (no licensed resources):

  OVS1   null variant in a bona fide tumour suppressor (NCG canonical TSG-only
         role, curated overrides; NMD / 3'-end check)
  OS1    same amino-acid change, different nucleotide change, as a ClinVar
         oncogenic / likely oncogenic variant (somatic oncogenicity, >= 1 star);
         the variant's own ClinVar classification is shown but not counted
  OS2    calibrated functional assay (MaveDB / ENIGMA Table 9 / TP53 VCEP data)
         showing loss of function at Strong or above, tumour suppressors only; or
         CIViC: an Oncogenic assertion coded OS2, or Functional evidence items
         (Gain of Function / Neomorphic in oncogenes, Loss of Function / Dominant
         Negative in TSGs) rated >= 3 from >= 2 publications
  OS3    cancerhotspots.org residue with >= 50 samples and >= 10 with this change
  OM1    critical region without benign variation (HUVar PM1: VCEP hotspot /
         domain tables or the ClinVar P/LP cluster)
  OM2    in-frame indel in an oncogene / TSG, or stop-loss in a TSG
  OM3    cancerhotspots.org residue with < 50 samples and >= 10 with this change
  OM4    different missense at a residue with a ClinVar oncogenic missense (>= 2 stars),
         Grantham distance at least ~ that of the comparator
  OP1    in-silico evidence for an effect (HUVar PP3)
  OP3    cancerhotspots.org residue / indel region, this change in < 10 samples
  OP4    absent from gnomAD or grpmax AF <= somatic_op4_max_af (1e-4; 5e-4 for known
         oncogenic variants: ClinVar oncogenic itself, or OS1 / OS3 / OM3 met)
  SBVS1  gnomAD grpmax AF > 5 %
  SBS1   gnomAD grpmax AF > 1 %
  SBS2   calibrated functional assay showing normal function at Strong or above, or
         CIViC Unaltered Function evidence (same thresholds) / assertion coded SBS2
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
        from acmg_classifier.somatic.civic import load_civic
        self._civic_path = getattr(cfg, "civic_somatic_tsv", None)
        self._civic_db = (load_civic(self._civic_path)
                          if getattr(cfg, "somatic_civic", True) and self._civic_path else {})

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
        known_self = self._clinvar_comparators(variant, ann, res)

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
        self._civic(ann, role, res)

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
        self._frequency(ann, res, known_self)

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

    def _clinvar_comparators(self, variant: VariantRecord, ann: AnnotationData, res) -> bool:
        """OS1 / OM4 from previously established variants in ClinVar: somatic
        oncogenicity classifications (Oncogenic / Likely oncogenic) only, always a
        different nucleotide change than the variant. Germline P/LP is not used."""
        pc = ann.rule_consequence or ann.primary_consequence
        if pc is None or pc.consequence != CT.MISSENSE:
            res[C.OS1].evidence = "Not a missense variant"
            res[C.OM4].evidence = "Not a missense variant"
            return False
        sp = ann.splice
        if sp is not None and sp.is_available and (sp.max_delta or 0) >= _SPLICE_CAUTION:
            msg = (f"Predicted to affect splicing ({sp.tool} {sp.max_delta:.2f}) - amino-acid "
                   "comparison not used")
            res[C.OS1].evidence = res[C.OM4].evidence = msg
            return False
        import re
        from acmg_classifier.criteria.grantham import grantham_distance
        from acmg_classifier.local_db.clinvar_sqlite import _aa_change_key

        stars = int(getattr(self._cfg, "somatic_clinvar_min_stars", 2))
        os1_onc = int(getattr(self._cfg, "somatic_os1_onc_min_stars", 1))
        key = _aa_change_key(pc.hgvs_p)
        comps, _, self_note = self._comparators(variant, pc, key, stars)

        def ok(c, need):
            return c["onc_stars"] >= need

        same_all = [c for c in comps if c["aa"] == key]
        same = [c for c in same_all if ok(c, os1_onc)]
        if same:
            c = same[0]
            res[C.OS1].met = True
            res[C.OS1].evidence = f"Same amino-acid change as {c['text']}"
        else:
            res[C.OS1].evidence = (f"No other ClinVar oncogenic / likely oncogenic variant "
                                   f"(>= {os1_onc} star) with the same amino-acid change")
        ref, _, alt = (pc.amino_acids or "/").partition("/")
        cand = grantham_distance(ref, alt)
        pos_txt = str(pc.protein_position)
        at_residue = [c for c in comps if c["aa"] and c["aa"] != key and c["aa"][-1] != "*"
                      and re.sub(r"\D", "", c["aa"]) == pos_txt]
        others = [c for c in at_residue if ok(c, stars)]
        best, seen = None, []
        for c in others:
            d = grantham_distance(c["aa"][0], c["aa"][-1])
            seen.append(d)
            if cand is None or d is None or cand >= 0.9 * d:
                best = (c, d)
                break
        if best is not None:
            c, d = best
            res[C.OM4].met = True
            res[C.OM4].evidence = (f"Different missense at residue {pos_txt}: {c['text']}; "
                                   f"Grantham {cand} vs {d}")
        elif others:
            res[C.OM4].evidence = (f"ClinVar established missense at residue {pos_txt} have a "
                                   f"larger amino-acid difference (Grantham {cand} vs "
                                   f"{', '.join(str(x) for x in seen[:5])})")
        else:
            res[C.OM4].evidence = (f"No ClinVar oncogenic / likely oncogenic missense "
                                   f"(>= {stars} stars) at residue {pos_txt}")
        weak = {C.OS1: [c["text"] for c in same_all if not ok(c, os1_onc)],
                C.OM4: [c["text"] for c in at_residue if not ok(c, stars)]}
        for crit, texts in weak.items():
            if texts and not res[crit].met:
                res[crit].evidence += f" (below the star threshold: {'; '.join(texts[:3])})"
        if self_note:
            res[C.OS1].evidence += f"; {self_note}"
        return bool(self_note)

    def _comparators(self, variant, pc, key, stars) -> tuple[list[dict], list[str], str]:
        """Previously established variants in the codon: [{aa, text}], notes on
        ClinVar oncogenic records below the star threshold, and a note on the
        variant's own ClinVar oncogenicity (reported, never counted: OS1 needs a
        different nucleotide change, as PS1)."""
        from acmg_classifier.local_db.clinvar_sqlite import (
            query_aa_change_at, query_oncogenic_at_codon,
        )
        from acmg_classifier.local_db.clinvar_vcf import ONCOGENIC_LABELS, query_clinvar_oncogenicity
        from acmg_classifier.utils.chrom import strip_chr

        db = self._cfg.clinvar_sqlite
        comps: list[dict] = []
        weak: list[str] = []
        self_note = ""

        def is_self(r) -> bool:
            return (strip_chr(str(r["chrom"])) == strip_chr(variant.chrom) and r["pos"] == variant.pos
                    and r["ref"] == variant.ref and r["alt"] == variant.alt)

        onc = query_oncogenic_at_codon(db, pc.gene_symbol, pc.protein_position,
                                       variant.chrom, variant.pos)
        if onc is None:  # DB predates the oncogenicity columns: ClinVar VCF if present
            onc = []
            for r in query_clinvar_oncogenicity(self._cfg.clinvar_vcf, variant.chrom,
                                                variant.pos - 2, variant.pos + 2):
                if r["oncogenicity"].lower() not in ONCOGENIC_LABELS:
                    continue
                if is_self(r):
                    r["aa"], r["hgvs_p"] = key, pc.hgvs_p
                else:
                    hit = query_aa_change_at(db, r["chrom"], r["pos"], r["ref"], r["alt"])
                    if hit is None or hit[0] != pc.gene_symbol or not hit[1]:
                        continue
                    r["aa"], r["hgvs_p"] = hit[1], hit[2]
                onc.append(r)
        for r in onc:
            label = (r.get("hgvs_p") or "").split(":")[-1] or r.get("aa") or ""
            text = (f"ClinVar {r['oncogenicity']} (oncogenicity) {label} (VCV {r['variation_id']}, "
                    f"{r['stars']} stars: {r['review_status']})")
            if is_self(r):
                self_note = (f"this variant is {text.replace(' ' + label, '', 1)} - not counted "
                             "(OS1 needs a different nucleotide change)")
                continue
            comps.append({"aa": r.get("aa"), "text": text, "onc_stars": r["stars"]})
        return comps, weak, self_note

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

    def _civic(self, ann: AnnotationData, role: Optional[GeneRole], res) -> None:
        """OS2 / SBS2 from CIViC (functional evidence items, SOP-based oncogenic
        assertions), combined with the calibrated-assay result above."""
        if not getattr(self._cfg, "somatic_civic", True):
            return
        if not self._civic_db:
            for c in (C.OS2, C.SBS2):
                if not res[c].met:
                    res[c].evidence += ("; CIViC table not installed (setup_data.py --only civic)"
                                        if self._civic_path and not self._civic_path.exists() else "")
            return
        pc = ann.primary_consequence
        if pc is None:
            return
        key = None
        if pc.consequence == CT.MISSENSE and pc.protein_position and pc.amino_acids:
            ref, _, alt = pc.amino_acids.partition("/")
            key = f"{ref}{pc.protein_position}{alt}" if ref and alt else None
        elif pc.consequence in _INFRAME:
            key = short_protein_change(pc.hgvs_p)
        if not key:
            return
        from acmg_classifier.somatic.civic import civic_call
        call = civic_call(self._civic_db.get((pc.gene_symbol, key), []),
                          role.role if role else None,
                          int(getattr(self._cfg, "somatic_civic_min_rating", 3)),
                          int(getattr(self._cfg, "somatic_civic_min_sources", 2)))
        if call is None:
            for c in (C.OS2, C.SBS2):
                if not res[c].met:
                    res[c].evidence += f"; no CIViC functional evidence for {pc.gene_symbol} {key}"
            return
        if call.criterion is None:
            for c in (C.OS2, C.SBS2):
                if not res[c].met:
                    res[c].evidence += f"; {call.evidence}"
            return
        crit = C.OS2 if call.criterion == "OS2" else C.SBS2
        other = C.SBS2 if crit == C.OS2 else C.OS2
        if res[other].met:  # calibrated assay says the opposite
            res[other].met = False
            res[other].evidence += f" - conflicts with {call.evidence}; neither counted"
            res[crit].evidence = f"{call.evidence} - conflicts with the calibrated assay; not counted"
            return
        if res[crit].met:
            res[crit].evidence += f"; also {call.evidence}"
        else:
            res[crit].met, res[crit].evidence = True, call.evidence
        res[other].evidence += f"; CIViC supports {crit.value} instead"

    def _frequency(self, ann: AnnotationData, res, known_self: bool = False) -> None:
        gd = ann.gnomad
        if gd is None:
            for c in (C.OP4, C.SBS1, C.SBVS1):
                res[c].evidence = "gnomAD data unavailable"
            return
        af = gd.popmax_af if gd.popmax_af is not None else gd.af
        absent = (gd.ac or 0) == 0 and not af
        cut = float(getattr(self._cfg, "somatic_op4_max_af", 1e-4))
        # Known oncogenic variants (ClinVar-established or cancerhotspots >= 10
        # samples) are seen at low frequency in gnomAD through clonal haematopoiesis
        # / mosaicism; they get the higher "pathogenic outlier" ceiling.
        outlier = known_self or res[C.OS1].met or res[C.OS3].met or res[C.OM3].met
        if outlier:
            cut = max(cut, float(getattr(self._cfg, "somatic_op4_outlier_max_af", 5e-4)))
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
            res[C.OP4].met, res[C.OP4].evidence = True, (
                f"{detail} <= {cut:g} (extremely low"
                + ("; known oncogenic variant, outlier threshold" if outlier else "") + ")")
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
