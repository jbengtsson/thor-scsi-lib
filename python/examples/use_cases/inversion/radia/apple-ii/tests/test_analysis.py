import math

import numpy as np
import pytest

from apple_ii.analysis import (
    HarmonicResult,
    QualityCriteria,
    analyze_fundamental,
    analyze_harmonics,
    evaluate_quality,
    polarization_metrics,
)


def synthetic_rows(ax=0.3, ay=0.3, phase_y=math.pi / 2, noise=0.0):
    period = 40.0
    z = np.linspace(-120, 120, 1001)
    k = 2 * np.pi / period
    rng = np.random.default_rng(7)
    bx = ax * np.cos(k * z)
    by = ay * np.cos(k * z - phase_y)
    if noise:
        bx = bx + rng.normal(0, noise, z.size)
        by = by + rng.normal(0, noise, z.size)
    return period, np.column_stack((z, bx, by, np.zeros_like(z)))


def test_fundamental_recovers_positive_circular_field():
    period, rows = synthetic_rows()
    bx, by = analyze_fundamental(rows, period)
    metrics = polarization_metrics(bx, by)
    assert bx.amplitude == pytest.approx(0.3, abs=1e-10)
    assert by.amplitude == pytest.approx(0.3, abs=1e-10)
    assert metrics["circularity"] == pytest.approx(1.0, abs=1e-12)
    assert metrics["s3"] == pytest.approx(1.0, abs=1e-12)


def test_negative_helicity_reverses_s3_sign():
    period, rows = synthetic_rows(phase_y=-math.pi / 2)
    bx, by = analyze_fundamental(rows, period)
    assert polarization_metrics(bx, by)["s3"] == pytest.approx(-1.0, abs=1e-12)


def test_linear_and_unequal_ellipse_metrics():
    period, rows = synthetic_rows(ax=0.4, ay=0.2, phase_y=0.0)
    bx, by = analyze_fundamental(rows, period)
    metrics = polarization_metrics(bx, by)
    assert metrics["ratio"] == pytest.approx(2.0)
    assert metrics["circularity"] == pytest.approx(0.0, abs=1e-12)


def test_metrics_are_clipped_to_physical_range():
    bx = HarmonicResult(1, 1.0, 0.0, 0.0, 0.0)
    by = HarmonicResult(1, 1.0, math.pi / 2, 0.0, 0.0)
    metrics = polarization_metrics(bx, by)
    assert 0.0 <= metrics["circularity"] <= 1.0
    assert -1.0 <= metrics["s3"] <= 1.0


def test_quality_gate_rejects_tiny_and_noisy_fields():
    tiny_bx = HarmonicResult(1, 1e-9, 0.0, 0.0, 0.0)
    tiny_by = HarmonicResult(1, 1e-9, math.pi / 2, 0.0, 0.0)
    assert not evaluate_quality(tiny_bx, tiny_by).passed

    noisy_bx = HarmonicResult(1, 0.3, 0.0, 0.0, 0.2)
    noisy_by = HarmonicResult(1, 0.3, math.pi / 2, 0.0, 0.2)
    assessment = evaluate_quality(noisy_bx, noisy_by, QualityCriteria(0.0, 0.1))
    assert not assessment.passed
    assert assessment.relative_residual > 0.1


def test_higher_harmonics_are_measured():
    period = 40.0
    z = np.linspace(-120, 120, 1201)
    k = 2 * np.pi / period
    rows = np.column_stack(
        (
            z,
            0.3 * np.cos(k * z) + 0.03 * np.cos(3 * k * z),
            0.2 * np.sin(k * z) + 0.02 * np.sin(3 * k * z),
            np.zeros_like(z),
        )
    )
    harmonics = analyze_harmonics(rows, period, (1, 3))
    assert harmonics[1][0].amplitude == pytest.approx(0.3, rel=1e-9)
    assert harmonics[3][0].amplitude == pytest.approx(0.03, rel=1e-9)
    assert harmonics[3][1].amplitude == pytest.approx(0.02, rel=1e-9)


@pytest.mark.parametrize(
    "rows, message",
    [
        ([[0, 0, 0]], "N x >=3"),
        (np.column_stack((np.arange(8), np.ones(8), np.full(8, np.nan))), "finite"),
    ],
)
def test_invalid_rows_are_rejected(rows, message):
    with pytest.raises(ValueError, match=message):
        analyze_fundamental(rows, 40.0)
