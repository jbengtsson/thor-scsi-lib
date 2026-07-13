import math

import pytest

from apple_ii.analysis import HarmonicResult, QualityAssessment
from apple_ii.field import FieldIntegral
from apple_ii.optimize import PhaseEvaluation, optimize_phase, scan_phase


def evaluation(phase):
    circularity = max(0.0, 1.0 - ((phase - 2.345) / 4.0) ** 2)
    harmonic = HarmonicResult(1, 0.3, 0.0, 0.0, 0.001)
    quality = QualityAssessment(True, 0.424, 0.0014, 0.0033, ())
    return PhaseEvaluation(
        phase,
        harmonic,
        harmonic,
        {"circularity": circularity, "phase_deg": 90.0, "ratio": 1.0, "s1": 0.0, "s2": 0.0, "s3": circularity},
        quality,
        {1: (harmonic, harmonic)},
        FieldIntegral((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), -1.0, 1.0),
    )


def test_scan_requires_at_least_two_steps():
    with pytest.raises(ValueError, match="phase_steps"):
        scan_phase(evaluation, -1, 1, 1)


def test_bounded_refinement_finds_peak():
    result = optimize_phase(
        evaluation,
        -10.0,
        10.0,
        coarse_steps=21,
        refinement_levels=4,
        refinement_steps=11,
    )
    assert result.best.phase_mm == pytest.approx(2.345, abs=0.002)
    assert -10.0 <= result.best.phase_mm <= 10.0


def test_optimizer_rejects_when_no_point_passes_quality_gate():
    def invalid(phase):
        item = evaluation(phase)
        return PhaseEvaluation(
            item.phase_mm,
            item.bx1,
            item.by1,
            item.metrics,
            QualityAssessment(False, 0.0, 0.0, math.inf, ("too small",)),
            item.harmonics,
            item.integrals,
        )

    with pytest.raises(ValueError, match="quality gate"):
        optimize_phase(invalid, -1.0, 1.0, coarse_steps=3, refinement_levels=0)
