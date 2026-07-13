from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Callable, Iterable, Sequence

import numpy as np

from .analysis import (
    HarmonicResult,
    QualityAssessment,
    QualityCriteria,
    analyze_harmonics,
    evaluate_quality,
    polarization_metrics,
)
from .field import FieldIntegral, field_integrals, sample_axis
from .geometry import Apple2Device


@dataclass(frozen=True)
class PhaseEvaluation:
    phase_mm: float
    bx1: HarmonicResult
    by1: HarmonicResult
    metrics: dict[str, float]
    quality: QualityAssessment
    harmonics: dict[int, tuple[HarmonicResult, HarmonicResult]]
    integrals: FieldIntegral

    @property
    def score(self) -> float:
        return self.metrics["circularity"] if self.quality.passed else -math.inf

    def flat_record(self) -> dict[str, float | bool | str]:
        return {
            "phase_mm": self.phase_mm,
            "Bx1_T": self.bx1.amplitude,
            "By1_T": self.by1.amplitude,
            "phase_deg": self.metrics["phase_deg"],
            "ratio": self.metrics["ratio"],
            "s1": self.metrics["s1"],
            "s2": self.metrics["s2"],
            "s3": self.metrics["s3"],
            "circularity": self.metrics["circularity"],
            "quality_passed": self.quality.passed,
            "total_amplitude_T": self.quality.total_amplitude_t,
            "relative_residual": self.quality.relative_residual,
            "quality_reasons": "; ".join(self.quality.reasons),
        }


@dataclass(frozen=True)
class PhaseOptimizationResult:
    best: PhaseEvaluation
    coarse_scan: tuple[PhaseEvaluation, ...]
    refinement_evaluations: tuple[PhaseEvaluation, ...]
    bounds_mm: tuple[float, float]


def _validate_scan_bounds(lower: float, upper: float, steps: int) -> None:
    if not math.isfinite(lower) or not math.isfinite(upper):
        raise ValueError("phase bounds must be finite")
    if upper <= lower:
        raise ValueError("phase_max must be greater than phase_min")
    if not isinstance(steps, int) or isinstance(steps, bool) or steps < 2:
        raise ValueError("phase_steps must be an integer >= 2")


def evaluate_phase(
    device: Apple2Device,
    phase_mm: float,
    *,
    z_min_mm: float,
    z_max_mm: float,
    n_points: int,
    central_periods: int,
    harmonic_orders: Sequence[int],
    quality_criteria: QualityCriteria,
) -> PhaseEvaluation:
    device.set_phase(float(phase_mm))
    rows = sample_axis(device, z_min_mm, z_max_mm, n_points)
    harmonics = analyze_harmonics(
        rows,
        device.parameters.period_mm,
        harmonic_orders=harmonic_orders,
        central_periods=central_periods,
    )
    if 1 not in harmonics:
        raise ValueError("harmonic_orders must include the fundamental order 1")
    bx1, by1 = harmonics[1]
    metrics = polarization_metrics(bx1, by1)
    quality = evaluate_quality(bx1, by1, quality_criteria)
    return PhaseEvaluation(
        phase_mm=float(phase_mm),
        bx1=bx1,
        by1=by1,
        metrics=metrics,
        quality=quality,
        harmonics=harmonics,
        integrals=field_integrals(rows),
    )


def scan_phase(
    evaluator: Callable[[float], PhaseEvaluation],
    phase_min_mm: float,
    phase_max_mm: float,
    phase_steps: int,
) -> tuple[PhaseEvaluation, ...]:
    """Evaluate an already-built device through a caller-supplied phase evaluator."""
    _validate_scan_bounds(phase_min_mm, phase_max_mm, phase_steps)
    return tuple(
        evaluator(float(phase))
        for phase in np.linspace(phase_min_mm, phase_max_mm, phase_steps)
    )


def _best_valid(evaluations: Iterable[PhaseEvaluation]) -> PhaseEvaluation:
    valid = [evaluation for evaluation in evaluations if math.isfinite(evaluation.score)]
    if not valid:
        raise ValueError("no phase point passed the field-amplitude and residual quality gate")
    return max(
        valid,
        key=lambda evaluation: (
            evaluation.score,
            evaluation.quality.total_amplitude_t,
            -evaluation.quality.relative_residual,
        ),
    )


def optimize_phase(
    evaluator: Callable[[float], PhaseEvaluation],
    phase_min_mm: float,
    phase_max_mm: float,
    *,
    coarse_steps: int = 41,
    refinement_levels: int = 3,
    refinement_steps: int = 11,
) -> PhaseOptimizationResult:
    """Bounded coarse-to-fine phase optimization without a SciPy dependency."""
    _validate_scan_bounds(phase_min_mm, phase_max_mm, coarse_steps)
    if not isinstance(refinement_levels, int) or refinement_levels < 0:
        raise ValueError("refinement_levels must be an integer >= 0")
    if refinement_levels and (
        not isinstance(refinement_steps, int)
        or isinstance(refinement_steps, bool)
        or refinement_steps < 3
    ):
        raise ValueError("refinement_steps must be an integer >= 3")

    cache: dict[float, PhaseEvaluation] = {}

    def cached(phase: float) -> PhaseEvaluation:
        key = float(phase)
        if key not in cache:
            cache[key] = evaluator(key)
        return cache[key]

    coarse_phases = np.linspace(phase_min_mm, phase_max_mm, coarse_steps)
    coarse = tuple(cached(float(phase)) for phase in coarse_phases)
    best = _best_valid(coarse)
    refinement: list[PhaseEvaluation] = []
    step_width = (phase_max_mm - phase_min_mm) / (coarse_steps - 1)

    for _ in range(refinement_levels):
        lower = max(phase_min_mm, best.phase_mm - step_width)
        upper = min(phase_max_mm, best.phase_mm + step_width)
        level = [cached(float(phase)) for phase in np.linspace(lower, upper, refinement_steps)]
        refinement.extend(level)
        best = _best_valid((*coarse, *refinement))
        step_width = (upper - lower) / (refinement_steps - 1)

    return PhaseOptimizationResult(
        best=best,
        coarse_scan=coarse,
        refinement_evaluations=tuple(refinement),
        bounds_mm=(float(phase_min_mm), float(phase_max_mm)),
    )
