"""APPLE-II RADIA reference-model toolkit."""

from .analysis import (
    HarmonicResult,
    QualityAssessment,
    QualityCriteria,
    analyze_fundamental,
    analyze_harmonics,
    evaluate_quality,
    polarization_metrics,
)
from .geometry import (
    Apple2Device,
    Apple2Parameters,
    ArrayMotion,
    BlockRecord,
    RowHandle,
    build_device,
)

__all__ = [
    "Apple2Device",
    "Apple2Parameters",
    "ArrayMotion",
    "BlockRecord",
    "HarmonicResult",
    "QualityAssessment",
    "QualityCriteria",
    "RowHandle",
    "analyze_fundamental",
    "analyze_harmonics",
    "build_device",
    "evaluate_quality",
    "polarization_metrics",
]

__version__ = "0.3.0"
