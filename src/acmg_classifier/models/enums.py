from enum import Enum


class Assembly(str, Enum):
    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"


class ACMGCriterion(str, Enum):
    # Pathogenic very strong
    PVS1 = "PVS1"
    # Pathogenic strong
    PS1 = "PS1"
    PS2 = "PS2"
    PS3 = "PS3"
    PS4 = "PS4"
    # Pathogenic moderate
    PM1 = "PM1"
    PM2 = "PM2"
    PM3 = "PM3"
    PM4 = "PM4"
    PM5 = "PM5"
    PM6 = "PM6"
    # Pathogenic supporting
    PP1 = "PP1"
    PP2 = "PP2"
    PP3 = "PP3"
    PP4 = "PP4"
    PP5 = "PP5"
    # Benign stand-alone
    BA1 = "BA1"
    # Benign strong
    BS1 = "BS1"
    BS2 = "BS2"
    BS3 = "BS3"
    BS4 = "BS4"
    # Benign supporting
    BP1 = "BP1"
    BP2 = "BP2"
    BP3 = "BP3"
    BP4 = "BP4"
    BP5 = "BP5"
    BP6 = "BP6"
    BP7 = "BP7"


class CriterionStrength(str, Enum):
    """7-level Bayesian strength scale (Tavtigian 2020 + Bergquist 2024)."""
    VERY_STRONG = "VeryStrong"       # ±8 points
    STRONG = "Strong"                # ±4 points
    THREE_POINT = "ThreePoint"       # ±3 points (Bergquist 2024 extension)
    MODERATE = "Moderate"            # ±2 points
    SUPPORTING = "Supporting"        # ±1 point
    NOT_MET = "NotMet"               # 0 (criterion not triggered)
    INDETERMINATE = "Indeterminate"  # score exists but falls in grey zone


class CriterionDirection(str, Enum):
    PATHOGENIC = "Pathogenic"
    BENIGN = "Benign"


class Pathogenicity(str, Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "LikelyPathogenic"
    VUS = "VUS"
    LIKELY_BENIGN = "LikelyBenign"
    BENIGN = "Benign"


class Inheritance(str, Enum):
    AD = "AD"
    AR = "AR"
    XL = "XL"
    UNKNOWN = "Unknown"


class InSilicoTool(str, Enum):
    ALPHAMISSENSE = "alphamissense"
    ESM1B = "esm1b"


class SpliceTool(str, Enum):
    NONE = "none"
    OPENSPLICEAI = "openspliceai"  # GPL-3.0; runtime inference via openspliceai CLI
    SPLICEAI = "spliceai"          # requires Illumina commercial license (opt-in only)
    # SQUIRLS retained for when its precomputed DB is downloadable again; not
    # the default and not currently selectable from the CLI.
    SQUIRLS = "squirls"
    # MMSplice DISABLED: its dependency chain (numpy<2, cyvcf2<=0.30.x, pyranges
    # 0.0.x) conflicts with this project's cyvcf2/numpy. The integration code is
    # retained (commented out across the codebase) to re-enable later.
    # MMSPLICE = "mmsplice"


class VariantType(str, Enum):
    SNV = "SNV"
    INDEL = "INDEL"
    MNV = "MNV"


class ConsequenceType(str, Enum):
    FRAMESHIFT = "frameshift_variant"
    STOP_GAINED = "stop_gained"
    STOP_LOST = "stop_lost"
    START_LOST = "start_lost"
    SPLICE_ACCEPTOR = "splice_acceptor_variant"
    SPLICE_DONOR = "splice_donor_variant"
    MISSENSE = "missense_variant"
    SYNONYMOUS = "synonymous_variant"
    INFRAME_INSERTION = "inframe_insertion"
    INFRAME_DELETION = "inframe_deletion"
    SPLICE_REGION = "splice_region_variant"
    INTRON = "intron_variant"
    FIVE_PRIME_UTR = "5_prime_UTR_variant"
    THREE_PRIME_UTR = "3_prime_UTR_variant"
    UPSTREAM = "upstream_gene_variant"
    DOWNSTREAM = "downstream_gene_variant"
    INTERGENIC = "intergenic_variant"
    TRANSCRIPT_ABLATION = "transcript_ablation"
    OTHER = "other"
