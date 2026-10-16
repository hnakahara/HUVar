#!/usr/bin/env python3
"""Build a per-gene PM1 hotspot table from the ClinGen cspec_summary.json export.

PM1 ("located in a mutational hotspot / critical functional domain") is defined
per gene in free text by each VCEP. This script mines those PM1 descriptions for
machine-readable hotspot evidence — residue RANGES (e.g. "codons 167-931",
"amino acids 271-292") and explicit RESIDUE positions (e.g. "Arg158", "R107",
"residues 175, 245, 248") — and writes one row per (gene, strength) to
``pm1_hotspots.tsv`` (columns: gene_symbol, strength, regions, residues).

Coverage is best-effort: rules expressed only as exon numbers, external tables
("Supp. Table 4"), or cancerhotspots.org occurrence counts cannot be resolved to
residue positions and are skipped. Multi-gene PM1 entries are skipped because a
region quoted in the prose usually applies to one gene, not all of them.

Usage:
    python scripts/build_pm1_hotspots.py \
        --summary resources/clingen/cspec_json/cspec_summary.json \
        --out resources/shared/pm1_hotspots.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import re

_AA3 = "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val"
_AA1 = "ACDEFGHIKLMNPQRSTVWY"

# Tokens that carry stray numbers (transcripts, PMIDs, citations, HTML, coding
# HGVS) — stripped before parsing so they are not mistaken for residues.
# NOTE: the bracket rule strips ONLY numeric-citation brackets ("[12]",
# "[1-3]", "[10, 11]"). VCEPs (RASopathy etc.) write the actual hotspot ranges
# inside square brackets as "[AA 10-17]" — those contain letters and MUST be
# preserved, or every bracketed-domain PM1 (KRAS/HRAS/PTPN11/RAF1/…) loses its
# regions and then trips the not-applicable fallback. So the bracket class is
# restricted to digits/whitespace/separators only.
_NOISE = re.compile(
    r"ENS[TPG]\d+(?:\.\d+)?|N[MPR]_?\d+(?:\.\d+)?|c\.-?\d+|PMID:?\s*\d+|pmid_\d+"
    r"|\[[\s\d,;.–-]+\]|<[^>]+>|&[a-z]+;",
    re.IGNORECASE,
)
# An "AA<n> - AA<m>" span (e.g. "Ser151 - Pro153") → residue range.
_AA_RANGE = re.compile(
    rf"(?:{_AA3}|[{_AA1}])\s?(\d{{1,4}})\s*[-–]\s*(?:{_AA3}|[{_AA1}])\s?(\d{{1,4}})"
)
# A bare numeric range "167-931". Occurrence-count ranges (cancerhotspots) are
# excluded by a trailing-context check below.
_RANGE = re.compile(r"(?<![\d.])(\d{1,4})\s*[-–]\s*(\d{1,4})(?!\d)")
_AA3_RES = re.compile(rf"\b(?:{_AA3})\s?(\d{{1,4}})\b")
_AA1_RES = re.compile(rf"\b([{_AA1}])(\d{{2,4}})\b")
# "AA <n>"-prefixed single residues (e.g. PTPN11 "AA 247, AA 251, AA 256").
# Ranges in this notation ("AA 10-17") are already captured by _RANGE; this
# picks up the isolated residues a VCEP lists alongside them.
_AA_PREFIX_RES = re.compile(r"\bAA\s?(\d{1,4})\b", re.IGNORECASE)
# A comma-separated bare-number list introduced by codon/residue wording.
_CODON_LIST = re.compile(
    r"(?:codons?|residues?|amino acids?)[^:.\n]*?[:\s]\s*((?:\d{1,4}\s*,\s*)+\d{1,4})",
    re.IGNORECASE,
)
# Occurrence/instance wording that follows a count range, not a residue range.
_OCCURRENCE = re.compile(r"\b(occurrence|instances?|somatic occurrence)", re.IGNORECASE)
_NOT_APPLICABLE = re.compile(
    r"does not apply|not applicable|highly polymorphic", re.IGNORECASE
)
# A positive applicability clause ("Applicable only to ... domains/residues").
# When present, a trailing "Not applicable to specific amino acid residues
# (see PM5)" is a PM5 redirect for the *rest* of the gene — NOT a gene-level
# negation — so the gene must NOT be marked not_applicable.
_POSITIVE_APPLICABLE = re.compile(r"\bapplicable\b\s+(?:only\s+)?to\b", re.IGNORECASE)

_VALID_STRENGTH = {"Supporting", "Moderate", "Strong"}


def _normalise(text: str) -> str:
    # Remove markdown escape backslashes ("NM\_00546.4", "PM1\_strong") so the
    # noise regex can match transcript/identifier tokens; drop thousands
    # separators ("2,101" -> "2101"); then strip noise tokens carrying unrelated
    # numbers (transcripts, PMIDs, citations).
    text = text.replace("\\", "")
    text = re.sub(r"(\d),(\d{3})\b", r"\1\2", text)
    return _NOISE.sub(" ", text)


def parse_regions(text: str) -> tuple[list[tuple[int, int]], list[int]]:
    t = _normalise(text)
    ranges: set[tuple[int, int]] = set()
    residues: set[int] = set()

    # b > 4 rejects enumeration / domain-numbering noise ("Exons 1-3",
    # "Cys2-Cys3", "1. ... 4.") — no real protein hotspot ends within the first
    # four residues.
    for m in _AA_RANGE.finditer(t):
        a, b = int(m.group(1)), int(m.group(2))
        if 0 < a < b <= 9999 and b > 4 and b - a < 2000:
            ranges.add((a, b))

    for m in _RANGE.finditer(t):
        a, b = int(m.group(1)), int(m.group(2))
        if not (0 < a < b <= 9999 and b > 4 and b - a < 2000):
            continue
        # Skip "2-9 occurrences" style count ranges (cancerhotspots wording).
        if _OCCURRENCE.search(t[m.end():m.end() + 30]):
            continue
        ranges.add((a, b))

    for m in _AA3_RES.finditer(t):
        residues.add(int(m.group(1)))
    for m in _AA1_RES.finditer(t):
        residues.add(int(m.group(2)))
    for m in _AA_PREFIX_RES.finditer(t):
        residues.add(int(m.group(1)))
    for m in _CODON_LIST.finditer(t):
        for n in re.findall(r"\d{1,4}", m.group(1)):
            residues.add(int(n))

    # Drop residues <5 (disulfide-bond numbering like "Cys2-Cys3", start codon)
    # and those already covered by a range, to keep the table compact and clean.
    residues = {
        r for r in residues
        if r >= 5 and not any(a <= r <= b for a, b in ranges)
    }
    return sorted(ranges), sorted(residues)


# Manually curated PM1 hotspots, applied AFTER mining (replacing the mined
# value for that (gene, strength) key). Used for specs whose residue/region list
# the free-text miner cannot reliably extract: long explicit residue tables, a
# range written as two endpoints, or a hotspot defined in a multi-gene panel
# (the miner only mines single-gene specs). Every value is transcribed verbatim
# from the cited ClinGen spec. (regions, residues).
_CURATED: dict[tuple[str, str], tuple[list[tuple[int, int]], list[int]]] = {
    # SCN1A/2A/3A/8A GN067-070 (Epilepsy Sodium Channel VCEP) — PM1_Moderate
    # applies ONLY to the Pathogenic Enriched Region residues listed in the VCEP's
    # external "PM1 Table" (PM1_Table_v2.20241205.xlsx), transcribed here as 16
    # per-gene residue ranges. A variant outside these ranges does NOT get PM1.
    ("SCN1A", "Moderate"): ([(212, 230), (247, 255), (411, 424), (859, 867), (879, 887), (889, 902), (904, 912), (931, 939), (979, 997), (1321, 1364), (1468, 1476), (1478, 1491), (1493, 1511), (1631, 1649), (1656, 1674), (1771, 1784)], []),
    ("SCN2A", "Moderate"): ([(213, 231), (248, 256), (413, 426), (850, 858), (870, 878), (880, 893), (895, 903), (922, 930), (970, 988), (1311, 1354), (1458, 1466), (1468, 1481), (1483, 1501), (1621, 1639), (1646, 1664), (1761, 1774)], []),
    ("SCN3A", "Moderate"): ([(212, 230), (247, 255), (412, 425), (851, 859), (871, 879), (881, 894), (896, 904), (923, 931), (971, 989), (1309, 1352), (1453, 1461), (1463, 1476), (1478, 1496), (1616, 1634), (1641, 1659), (1756, 1769)], []),
    ("SCN8A", "Moderate"): ([(216, 234), (251, 259), (399, 412), (844, 852), (864, 872), (874, 887), (889, 897), (916, 924), (964, 982), (1301, 1344), (1449, 1457), (1459, 1472), (1474, 1492), (1612, 1630), (1637, 1655), (1751, 1764)], []),
    # RYR1 — two VCEPs define PM1 differently and BOTH apply:
    #  * Congenital Myopathies VCEP (GN150/GN179): PM1_Moderate in the
    #    pore/transmembrane region (aa 4800-4950).
    #  * Malignant Hyperthermia Susceptibility VCEP (GN012): PM1_Supporting in the
    #    N-terminal (1-552), central (2101-2458) and C-terminal (4631-4991)
    #    regions.
    # The mined table merged these into a single (wrong) Moderate row, so RYR1 is
    # pinned exclusively (see _CURATED_EXCLUSIVE) to these two curated rows.
    # GN012 (MHS VCEP) text: N-terminal 1-552 and central 2101-2458 are
    # PM1_Moderate ("if PS1/PM5 applicable then PM1 should be used at
    # supporting" — downgrade applied in the registry, PM1_DOWNGRADE_WITH);
    # C-terminal 4631-4991 is PM1_Supporting. GN150/GN179 (Congenital
    # Myopathies VCEP): 4800-4950 Moderate (highest strength wins).
    ("RYR1", "Moderate"): ([(1, 552), (2101, 2458), (4800, 4950)], []),
    ("RYR1", "Supporting"): ([(4631, 4991)], []),
    # OTC GN156 v1.0.0 — CP-binding, ornithine, catalytic + conserved residues
    # (the miner captured only Met-268). 21 critical residues.
    ("OTC", "Moderate"): ([], [
        90, 91, 92, 93, 117, 141, 163, 168, 171, 198, 199, 263, 267, 268, 269,
        277, 302, 303, 304, 305, 330,
    ]),
    # GALT GN158 v1.0.0 — active site Phe171–Gln188 (contiguous range; the miner
    # stored only the two endpoints as residues).
    ("GALT", "Moderate"): ([(171, 188)], []),
    # KCNQ4 GN005 v2.0.0 — pore-forming region aa 271-292 (NM_004700.4). In the
    # multi-gene Hearing Loss panel, so not mined.
    ("KCNQ4", "Moderate"): ([(271, 292)], []),
    # KCNQ1 GN112 v1.0.0 — pore helix aa 300-320. (Spec also requires
    # PM2_Supporting, a co-criterion the PM1 engine does not model.)
    ("KCNQ1", "Moderate"): ([(300, 320)], []),
    # JAK3 GN121 v2.3.0 — JH2 domain residues R651, C759.
    ("JAK3", "Moderate"): ([], [651, 759]),
    # PDHA1 GN014 v1.0.0 — TPP-binding + αβ / α2β2 interface + phospho-loop
    # residues (multi-gene panel, not mined). 59 critical residues.
    ("PDHA1", "Moderate"): ([], [
        88, 118, 119, 140, 160, 162, 164, 165, 166, 167, 169, 172, 173, 176,
        177, 179, 180, 183, 195, 196, 197, 198, 199, 200, 201, 202, 203, 205,
        209, 210, 213, 225, 227, 228, 229, 230, 231, 245, 287, 288, 289, 290,
        291, 292, 293, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305,
        314, 315, 316,
    ]),
    # VHL GN078 v1.1.0 — germline + somatic hotspots (the residue list lives in
    # the VCEP's external "Germline and Somatic Hotspots" table, not the JSON).
    # PM1_Moderate for the curated hotspot residues. (The somatic <10-instances →
    # PM1_Supporting split needs per-residue cancerhotspots counts, not encoded.)
    ("VHL", "Moderate"): ([], [
        # Germline hotspots, VCEP "Germline and Somatic Hotspots" table:
        # Stebbins et al. (R167, C162, L178, Y98, N78, P86) and Chiorean et al.
        # 2022 (65, 76, 78, 80, 86, 88, 96, 98, 112, 117, 161, 162, 167, 170,
        # 176). The struck-through somatic entries (S65, N78, S80, P86, W88,
        # Y112, W117, R161, C162) are used as germline hotspots.
        65, 76, 78, 80, 86, 88, 96, 98, 112, 117, 161, 162, 167, 170, 176, 178,
    ]),
    # LDLR GN013 v1.2.0 — missense in exon 4 (MANE NM_000527.5 codons 105-232) OR
    # one of the 60 conserved disulfide-bond cysteine residues (Supp. Table S4 of
    # the LDLR VCEP). PM1_Moderate. (Spec also requires PM2 — a co-criterion the
    # PM1 engine does not model.)
    ("LDLR", "Moderate"): ([(105, 232)], [
        27, 34, 39, 46, 52, 63, 68, 75, 82, 89, 95, 104, 109, 116, 121, 128,
        134, 143, 148, 155, 160, 167, 173, 184, 197, 204, 209, 216, 222, 231,
        236, 243, 248, 255, 261, 270, 276, 284, 289, 296, 302, 313, 318, 325,
        329, 338, 340, 352, 358, 364, 368, 377, 379, 392, 667, 677, 681, 696,
        698, 711,
    ]),
    # FBN1 GN022 v1.0.0 — Marfan VCEP. Cys residues in cbEGF-like (EGF-like
    # calcium-binding) domains -> PM1_Strong; Cys in EGF-like / TB / hybrid
    # domains AND the critical Gly residues (between Cys2-Cys3 of every cbEGF, and
    # between Cys3-Cys4 of cbEGF domains that have an upstream cbEGF) -> PM1_
    # Moderate; the cbEGF calcium-binding / hydroxylation / consensus-substitution
    # residues are also Moderate. Residues derived from UniProt P35555
    # (= MANE NM_000138.5, 2871 aa) domain + cysteine positions. The Ca-binding
    # consensus (Marfan review, Robinson & Booms; X-D-X-(N/D)-E-C…C-X-N*-X2-G-X-
    # (Y/F)…) maps relative to the domain cysteines: D=C1-4, N/D=C1-2, E=C1-1,
    # N*=C3+2 (hydroxylated), Y/F=C3+7 (validated 43/43 cbEGF). Cys-creating
    # variants are handled separately (alt=Cys). The N*→Ser / Gly→Ala tolerated-
    # exception caveat is not modelled.
    ("FBN1", "Strong"): ([], [
        250, 257, 262, 271, 273, 286, 292, 299, 304, 313, 315, 328, 494, 499,
        504, 513, 515, 528, 534, 541, 546, 555, 557, 570, 576, 582, 587, 596,
        598, 611, 617, 623, 628, 637, 639, 652, 727, 734, 739, 748, 750, 763,
        769, 776, 781, 790, 792, 805, 811, 816, 821, 830, 832, 845, 914, 921,
        926, 935, 937, 950, 1032, 1039, 1044, 1053, 1055, 1068, 1074, 1081,
        1086, 1095, 1097, 1111, 1117, 1124, 1129, 1138, 1140, 1153, 1159, 1166,
        1171, 1180, 1182, 1195, 1201, 1208, 1212, 1221, 1223, 1236, 1242, 1249,
        1254, 1263, 1265, 1278, 1284, 1291, 1296, 1305, 1307, 1320, 1326, 1333,
        1339, 1348, 1350, 1361, 1367, 1374, 1380, 1389, 1391, 1402, 1408, 1415,
        1420, 1429, 1431, 1444, 1450, 1456, 1461, 1470, 1472, 1485, 1491, 1497,
        1502, 1511, 1513, 1526, 1610, 1617, 1622, 1631, 1633, 1646, 1652, 1658,
        1663, 1672, 1674, 1687, 1770, 1777, 1782, 1791, 1793, 1806, 1812, 1818,
        1824, 1833, 1835, 1847, 1853, 1860, 1865, 1874, 1876, 1889, 1895, 1900,
        1905, 1914, 1916, 1928, 1934, 1942, 1947, 1956, 1958, 1971, 1977, 1984,
        1989, 1998, 2000, 2011, 2017, 2024, 2029, 2038, 2040, 2053, 2131, 2137,
        2142, 2151, 2153, 2164, 2170, 2176, 2181, 2190, 2192, 2204, 2210, 2217,
        2221, 2230, 2232, 2245, 2251, 2258, 2265, 2274, 2276, 2289, 2295, 2302,
        2307, 2316, 2318, 2331, 2406, 2413, 2418, 2427, 2429, 2442, 2448, 2455,
        2459, 2468, 2470, 2483, 2489, 2496, 2500, 2509, 2511, 2522, 2528, 2535,
        2541, 2550, 2552, 2565, 2571, 2577, 2581, 2590, 2592, 2605, 2611, 2617,
        2622, 2631, 2633, 2646, 2652, 2659, 2663, 2672, 2674, 2686,
    ]),
    ("FBN1", "Moderate"): ([], [
        85, 89, 94, 100, 102, 111, 119, 123, 129, 134, 136, 145, 150, 154, 160,
        166, 168, 177, 186, 195, 204, 209, 210, 221, 224, 231, 246, 248, 249,
        259, 260, 264, 269, 288, 290, 291, 301, 302, 306, 311, 336, 345, 358,
        359, 360, 365, 377, 389, 453, 460, 465, 474, 476, 488, 490, 492, 493,
        501, 502, 506, 509, 511, 530, 532, 533, 544, 548, 551, 553, 572, 574,
        575, 585, 589, 592, 594, 613, 615, 616, 626, 630, 633, 635, 661, 670,
        683, 684, 685, 696, 699, 711, 723, 725, 726, 737, 741, 744, 746, 765,
        767, 768, 779, 783, 786, 788, 807, 809, 810, 819, 823, 826, 828, 853,
        862, 875, 876, 887, 890, 896, 910, 912, 913, 924, 928, 931, 933, 958,
        967, 980, 981, 982, 993, 996, 1008, 1028, 1030, 1031, 1042, 1046, 1049,
        1051, 1070, 1072, 1073, 1082, 1084, 1088, 1091, 1093, 1113, 1115, 1116,
        1126, 1127, 1131, 1134, 1136, 1155, 1157, 1158, 1169, 1173, 1176, 1178,
        1197, 1199, 1200, 1214, 1217, 1219, 1238, 1240, 1241, 1251, 1252, 1256,
        1259, 1261, 1280, 1282, 1283, 1294, 1298, 1301, 1303, 1322, 1324, 1325,
        1334, 1341, 1344, 1346, 1363, 1365, 1366, 1382, 1385, 1387, 1404, 1406,
        1407, 1416, 1418, 1422, 1425, 1426, 1427, 1446, 1448, 1449, 1459, 1463,
        1466, 1468, 1487, 1489, 1490, 1500, 1504, 1507, 1509, 1534, 1549, 1562,
        1563, 1564, 1574, 1577, 1589, 1606, 1608, 1609, 1619, 1620, 1624, 1627,
        1629, 1648, 1650, 1651, 1659, 1661, 1665, 1668, 1670, 1695, 1706, 1719,
        1720, 1721, 1733, 1736, 1748, 1766, 1768, 1769, 1780, 1784, 1787, 1789,
        1808, 1810, 1811, 1826, 1829, 1831, 1849, 1851, 1852, 1863, 1867, 1870,
        1872, 1891, 1893, 1894, 1901, 1903, 1907, 1910, 1912, 1930, 1932, 1933,
        1945, 1949, 1952, 1954, 1973, 1975, 1976, 1987, 1991, 1994, 1996, 2013,
        2015, 2016, 2027, 2031, 2034, 2036, 2061, 2070, 2083, 2084, 2085, 2096,
        2099, 2111, 2127, 2129, 2130, 2140, 2144, 2147, 2149, 2166, 2168, 2169,
        2177, 2179, 2183, 2186, 2187, 2188, 2206, 2208, 2209, 2223, 2226, 2228,
        2247, 2249, 2250, 2267, 2270, 2272, 2291, 2293, 2294, 2305, 2309, 2312,
        2314, 2339, 2348, 2363, 2364, 2365, 2375, 2378, 2390, 2402, 2404, 2405,
        2416, 2420, 2423, 2425, 2444, 2446, 2447, 2461, 2464, 2466, 2485, 2487,
        2488, 2502, 2505, 2506, 2507, 2524, 2526, 2527, 2536, 2539, 2543, 2546,
        2548, 2567, 2569, 2570, 2580, 2583, 2586, 2587, 2588, 2607, 2609, 2610,
        2618, 2619, 2624, 2627, 2629, 2648, 2650, 2651, 2662, 2665, 2668, 2669,
        2670,
    ]),
    # Brain Malformations VCEP (MTOR/AKT3/PIK3CA/PIK3R2; multi-gene spec, so the
    # miner skips it) — PM1_Supporting only; domains from the VCEP gene table (MTOR NM_004958.3,
    # AKT3 NM_005465.4, PIK3CA NM_006218.3, PIK3R2 NM_005027.3).
    ("MTOR", "Supporting"): ([(1382, 1982), (2015, 2114)], []),
    ("AKT3", "Supporting"): ([(5, 109), (151, 388), (425, 475)], []),
    ("PIK3CA", "Supporting"): ([(31, 108), (173, 292), (322, 483), (797, 1068)], []),
    ("PIK3R2", "Supporting"): ([(328, 716)], []),
    # IL2RG (SCID VCEP) — PM1_Strong: conserved Cys62/72/102/115; CpG Arg224,
    # Arg226, Arg285 (the spec's c. numbers use a legacy cDNA numbering, +14 nt);
    # WSxWS motif Trp237-Ser241. Transmembrane residues 263-283 only when a
    # charged/polar residue is introduced: see _CURATED_ALT. Co-requirements
    # (PM2 met; BA1/BS1/BS2 not met): PM1_CO_REQUIREMENTS in pm1_hotspots.py.
    ("IL2RG", "Strong"): ([(237, 241)], [62, 72, 102, 115, 224, 226, 285]),
    # HBA2 (Hemoglobinopathy VCEP) — AHSP-binding residues α32R, α104H, α118F,
    # α120P (HGVS numbering). The poly(A) signal is non-coding (not a missense
    # PM1 target).
    ("HBA2", "Moderate"): ([], [32, 104, 118, 120]),
    # ---- Corrections after re-checking every PM1 rule against the CSpec text
    # (2026-10, Genome Medicine revision). Positions are on the MANE Select
    # transcript HUVar annotates; exon-defined regions were converted with the
    # MANE v1.5 exon structure (resources/shared/mane_exons.tsv).
    # MECP2 GN143 (Rett/AS VCEP) defines MBD aa 90-162 and TRD aa 302-306 on
    # NM_004992.4 (e2, 486 aa); MANE Select is NM_001110792.2 (e1, 498 aa),
    # +12 residues (offset confirmed on eRepo records) -> 102-174, 314-318.
    ("MECP2", "Moderate"): ([(102, 174), (314, 318)], []),
    # UBE3A GN144: "3' cysteine binding site aa 820" is the catalytic Cys in
    # the 852-aa isoform (NM_130838; C843 in the 875-aa isoform). MANE
    # NM_130839.5 (872 aa) is +20 (offset confirmed on eRepo records) -> 840.
    ("UBE3A", "Moderate"): ([], [840]),
    # GP1BB: the miner dropped the last residue of "93, 95, 118, 141, and 147".
    ("GP1BB", "Moderate"): ([], [93, 95, 118, 141, 147]),
    # F9 Moderate: "variant is within exons 3, 4 or 5" (NM_000133.4: p.85-174).
    ("F9", "Moderate"): ([(85, 174)], []),
    # HNF4A Supporting DNA binding is codons 37-113 ("amino acids 33-36 do not
    # bind DNA"); the miner also kept 33-36 / 33-113.
    ("HNF4A", "Supporting"): ([(37, 113), (180, 220), (300, 350)], []),
    # TNNT2: codons 89-189 on NM_001276345.2 (= MANE); 79-179 is the same
    # region on ENST00000367318 and must not be added on MANE coordinates.
    ("TNNT2", "Supporting"): ([(89, 189)], []),
    # ENG: individual residues only; "C350-C382 disulfide" is not a range.
    ("ENG", "Moderate"): ([], [207, 278, 282, 350, 363, 382, 394, 412, 516, 549, 582]),
    # RPE65: active-site His/Glu residues plus the Ala107-Gly125 range.
    ("RPE65", "Moderate"): ([(107, 125)], [180, 182, 241, 313, 417, 527]),
    # GUCY2D Supporting: other missense in the catalytic domain p.873-951.
    ("GUCY2D", "Supporting"): ([(873, 951)], []),
    # DICER1 Supporting: RNase IIIb domain p.Y1682-p.S1846 (range, not 2 residues).
    ("DICER1", "Supporting"): ([(1682, 1846)], []),
    # ACADVL GN021 "table PM1" (HGVS numbering; Moderate): p.1-40, 214-223,
    # 249-251, R326, 381-382, R429, E441, R459, 481-516, 562.
    ("ACADVL", "Moderate"): ([(1, 40), (214, 223), (249, 251), (381, 382), (481, 516)],
                             [326, 429, 441, 459, 562]),
    # PTPN11: AA 4 was dropped by the miner's residue>=5 filter.
    ("PTPN11", "Moderate"): ([(7, 9), (58, 63), (69, 77), (278, 281)],
                             [4, 247, 251, 255, 256, 258, 261, 265, 284]),
    # F8 Moderate: FXa-binding Gly2267-Gly2304 "with the exception of Ser2283".
    ("F8", "Moderate"): ([(2267, 2282), (2284, 2304)], [1332, 1667]),
    # BRAF: exon 6 (p.238-287), exon 11 (p.439-478, incl. P-loop 459-474) and
    # CR3 activation segment 594-627 (NM_004333.6).
    ("BRAF", "Moderate"): ([(238, 287), (439, 478), (594, 627)], []),
    # RAF1: CR2 251-266 (exon 7), exon 14 (p.473-512) and exon 17 (p.602-648)
    # (NM_002880.4).
    ("RAF1", "Moderate"): ([(251, 266), (473, 512), (602, 648)], []),
    # FBN1 Cys-creating variants: a missense introducing a new cysteine anywhere
    # in a disulfide-bonded domain (EGF-like / cbEGF / TB / hybrid) → PM1_Moderate
    # (handled by the evaluator as alt=Cys within these ranges). UniProt P35555
    # domain ranges (merged).
    ("FBN1", "cys_creating"): ([
        (81, 112), (115, 178), (184, 236), (246, 329), (334, 389), (449, 653),
        (659, 711), (723, 846), (851, 902), (910, 951), (956, 1008),
        (1028, 1527), (1532, 1589), (1606, 1688), (1693, 1748), (1766, 2054),
        (2059, 2111), (2127, 2332), (2337, 2390), (2402, 2687),
    ], []),
}


# Curated rows restricted to substitutions introducing particular residues
# (written with an alt_aa column). gene -> (strength, regions, residues, alt_aa).
_CURATED_ALT: dict[str, list[tuple[str, list[tuple[int, int]], list[int], str]]] = {
    # IL2RG transmembrane domain: charged or polar residue introduced
    # (Asn, Asp, Arg, Cys, His, Glu, Gln, Lys, Ser, Thr, Tyr).
    "IL2RG": [("Strong", [(263, 283)], [], "NDRCHEQKSTY")],
}


# TP53 (GN009): PM1_Moderate also "for germline missense variants seen in
# cancerhotspots.org with >= 10 somatic occurrences for the same amino acid
# change"; PM1_Supporting for 2-9 occurrences. Built from the cancerhotspots.org
# single-residue export (resources/clingen/cancerhotspots_TP53.tsv) as one row
# per residue and strength, restricted to the qualifying alternate residues.
#
# VHL (GN078): PM1_Moderate for somatic variants with >=10 instances of the same
# amino-acid change in cancerhotspots.org, PM1_Supporting for <10 instances —
# restricted to the somatic hotspot residues of the VCEP "Germline and Somatic
# Hotspots" table that are not struck through (the struck ones are germline
# hotspots, Moderate for any missense via _CURATED). In the cancerhotspots.org
# export no VHL missense change reaches 10 samples, so these are Supporting for
# the changes observed (>=1 sample).
_CANCERHOTSPOTS = {
    # gene: (export, min samples for Supporting, residues allowed or None)
    "TP53": ("resources/clingen/cancerhotspots_TP53.tsv", 2, None),
    "VHL": ("resources/clingen/cancerhotspots_VHL.tsv", 1,
            {68, 74, 89, 111, 114, 115, 121, 135, 151, 158, 169}),
}


def cancerhotspot_rows(path: str, sup_min: int = 2,
                       residues: set[int] | None = None,
                       ) -> list[tuple[str, str, list, list, str]]:
    """``(gene, strength, ranges, residues, alt_aa)`` rows from a cancerhotspots
    export: per residue, the missense alternates with >=10 samples (Moderate)
    and with ``sup_min``-9 samples (Supporting), optionally restricted to
    ``residues``. Nonsense (*), splice (sp), synonymous and in-frame indel rows
    are ignored."""
    import os
    if not os.path.exists(path):
        return []
    out = []
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("Type") or "") != "single residue":
                continue
            m = re.fullmatch(r"([A-Z])(\d+)", (r.get("Residue") or "").strip())
            if not m:
                continue
            ref, pos = m.group(1), int(m.group(2))
            if residues is not None and pos not in residues:
                continue
            mod, sup = [], []
            for tok in (r.get("Variants") or "").split("|"):
                alt, _, n = tok.partition(":")
                if len(alt) != 1 or alt not in _AA1 or alt == ref or not n.isdigit():
                    continue
                (mod if int(n) >= 10 else sup if int(n) >= sup_min else []).append(alt)
            if mod:
                out.append((r["Gene"], "Moderate", [], [pos], "".join(sorted(mod))))
            if sup:
                out.append((r["Gene"], "Supporting", [], [pos], "".join(sorted(sup))))
    return out

# Curated not-applicable genes from multi-gene specs (skipped by the miner):
#  * Hearing Loss VCEP: PM1 applies only to the KCNQ4 pore-forming region
#    (271-292); every other gene of the panel has no PM1 hotspot.
#  * HBB (Hemoglobinopathy VCEP): PM1 covers only the TATA box and the poly(A)
#    signal (non-coding), so no missense variant qualifies.
_CURATED_NOT_APPLICABLE = {
    "CDH23", "COCH", "GJB2", "MYO15A", "MYO6", "MYO7A", "OTOF", "SLC26A4",
    "TECTA", "USH2A", "HBB",
}


def build(summary_path: str) -> dict[tuple[str, str], tuple[set, set]]:
    with open(summary_path, encoding="utf-8") as fh:
        data = json.load(fh)["data"]

    # (gene, strength) -> (ranges set, residues set); plus not_applicable genes.
    table: dict[tuple[str, str], tuple[set, set]] = {}
    not_applicable: set[str] = set()

    for item in data:
        genes = [
            g.get("label") if isinstance(g, dict) else g
            for g in (item.get("genes") or [])
        ]
        genes = [g for g in genes if g]
        if len(genes) != 1:  # single-gene PM1 only (avoid mis-assigning regions)
            continue
        gene = genes[0]
        for code in item.get("codes", []):
            if code.get("label") != "PM1":
                continue
            strength = code.get("strengthDescriptor")
            if strength not in _VALID_STRENGTH:
                continue
            text = (code.get("text") or "").strip()
            if not text:
                continue
            ranges, residues = parse_regions(text)
            if not ranges and not residues:
                # No resolvable region. Record an explicit not-applicable only
                # when the VCEP negates the rule for the WHOLE gene — not when
                # "not applicable" is merely a PM5 redirect sub-clause sitting
                # next to a positive "Applicable only to … domains" statement.
                if _NOT_APPLICABLE.search(text) and not _POSITIVE_APPLICABLE.search(text):
                    not_applicable.add(gene)
                continue
            key = (gene, strength)
            r, res = table.setdefault(key, (set(), set()))
            r.update(ranges)
            res.update(residues)

    # Genes whose PM1 is defined SOLELY by the curated override: drop every mined
    # row first so a conflicting broader spec cannot leak extra regions/strengths.
    # RYR1 — the pore-region spec (GN150/GN179, 4800-4950) supersedes the broad
    # Malignant Hyperthermia VCEP regions (GN012).
    _CURATED_EXCLUSIVE = {"RYR1"}
    for g in _CURATED_EXCLUSIVE:
        for key in [k for k in table if k[0] == g]:
            del table[key]

    # Apply manual curation overrides (replace the mined value for that key).
    for key, (ranges, residues) in _CURATED.items():
        table[key] = (set(ranges), set(residues))

    # Curated not-applicable: ITGA2B/ITGB3 (GN011) declare PM1 "does not apply due
    # to genes being highly polymorphic". They sit in a MULTI-gene spec, which the
    # single-gene miner above skips, so force them here.
    not_applicable |= {"ITGA2B", "ITGB3"} | _CURATED_NOT_APPLICABLE

    # Materialise not_applicable as its own strength row (only if the gene has no
    # positive hotspot rows from any spec).
    positive_genes = {g for g, _ in table}
    for gene in not_applicable - positive_genes:
        table[(gene, "not_applicable")] = (set(), set())
    return table


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="resources/clingen/cspec_json/cspec_summary.json")
    ap.add_argument("--out", default="resources/shared/pm1_hotspots.tsv")
    ap.add_argument("--cancerhotspots", default=None,
                    help="override the TP53 cancerhotspots.org export")
    args = ap.parse_args()

    table = build(args.summary)
    rows = []
    for (gene, strength), (ranges, residues) in sorted(table.items()):
        rows.append({
            "gene_symbol": gene,
            "strength": strength,
            "regions": ";".join(f"{a}-{b}" for a, b in sorted(ranges)),
            "residues": ",".join(str(r) for r in sorted(residues)),
            "alt_aa": "",
        })
    alt_rows = [(g, *t) for g, lst in _CURATED_ALT.items() for t in lst]
    for gene, (path, sup_min, res) in _CANCERHOTSPOTS.items():
        if gene == "TP53" and args.cancerhotspots:
            path = args.cancerhotspots
        alt_rows += cancerhotspot_rows(path, sup_min, res)
    for gene, strength, ranges, residues, alt in sorted(alt_rows):
        rows.append({
            "gene_symbol": gene, "strength": strength,
            "regions": ";".join(f"{a}-{b}" for a, b in sorted(ranges)),
            "residues": ",".join(str(r) for r in sorted(residues)), "alt_aa": alt,
        })
    rows.sort(key=lambda r: (r["gene_symbol"], r["strength"], r["alt_aa"]))
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["gene_symbol", "strength", "regions", "residues", "alt_aa"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    n_genes = len({r["gene_symbol"] for r in rows})
    print(f"PM1 hotspot rows: {len(rows)} | genes: {n_genes} | written → {args.out}")


if __name__ == "__main__":
    main()
