class ACMGClassifierError(Exception):
    """Base exception for all classifier errors."""


class DataNotFoundError(ACMGClassifierError):
    """Required local data file or database is missing."""


class VEPRunError(ACMGClassifierError):
    """VEP subprocess failed."""


class AnnotationError(ACMGClassifierError):
    """Failed to annotate a variant."""


class ClassificationError(ACMGClassifierError):
    """Failed to classify a variant."""


class SetupError(ACMGClassifierError):
    """Setup or data download failed."""


class SupplementParseError(ACMGClassifierError):
    """Failed to parse the supplement TSV."""
