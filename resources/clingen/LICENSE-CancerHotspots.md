# Cancer Hotspots data (ODbL v1.0)

The following files contain data from **Cancer Hotspots**, Memorial Sloan Kettering
Cancer Center (https://www.cancerhotspots.org/#/home):

| File | Content |
|---|---|
| `resources/clingen/cancerhotspots_TP53.tsv` | Unmodified export for TP53 (retrieved 2026-10-17) |
| `resources/clingen/cancerhotspots_VHL.tsv` | Unmodified export for VHL (retrieved 2026-10-17) |
| `resources/shared/pm1_hotspots.tsv` (rows for TP53 and VHL) | Derived database: residues and amino-acid changes meeting the PM1 hotspot rules of the ClinGen TP53 and VHL VCEPs, built by `scripts/build_pm1_hotspots.py` |

Cancer Hotspots data are made available under the **ODC Open Database License
(ODbL) v1.0**: https://opendatacommons.org/licenses/odbl/1-0/ . The files above,
including the derived TP53/VHL rows, are redistributed under the same licence. You
are free to share and adapt them, provided that you attribute Cancer Hotspots, keep
any adapted database open and offer it under the ODbL.

The other rows of `pm1_hotspots.tsv` are derived from the ClinGen Criteria
Specification Registry (CC0 1.0). The HUVar program code is licensed separately
under the Apache License 2.0 (see `LICENSE`); the ODbL does not apply to it.
Variant classifications produced by HUVar are Produced Works under the ODbL and
require only attribution to Cancer Hotspots when they rely on these data.

Please cite:
- Chang MT, et al. Accelerating discovery of functional mutant alleles in cancer.
  Cancer Discov. 2018;8:174-183.
- Chang MT, et al. Identifying recurrent mutations in cancer reveals widespread
  lineage diversity and mutational specificity. Nat Biotechnol. 2016;34:155-163.
