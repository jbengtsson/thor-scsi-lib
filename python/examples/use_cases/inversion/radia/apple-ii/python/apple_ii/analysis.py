from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable, Sequence

import numpy as np

from ._util import finite_number, positive_number


@dataclass(frozen=True)
class HarmonicResult:
    order: int
    amplitude: float
    phase_rad: float
    offset: float
    rms_residual: float


@dataclass(frozen=True)
class QualityCriteria:
    min_total_amplitude_t: float = 1.0e-4
    max_relative_residual: float = 0.10

    def validate(self) -> None:
        if self.min_total_amplitude_t < 0 or not math.isfinite(self.min_total_amplitude_t):
            raise ValueError("min_total_amplitude_t must be finite and >= 0")
        if self.max_relative_residual < 0 or not math.isfinite(self.max_relative_residual):
            raise ValueError("max_relative_residual must be finite and >= 0")


@dataclass(frozen=True)
class QualityAssessment:
    passed: bool
    total_amplitude_t: float
    combined_rms_residual_t: float
    relative_residual: float
    reasons: tuple[str, ...]


def _validated_array(rows: Iterable[Sequence[float]]) -> np.ndarray:
    array = np.asarray(list(rows), dtype=float)
    if array.ndim != 2 or array.shape[1] < 3 or array.shape[0] < 8:
        raise ValueError("rows must be an N x >=3 array with N >= 8")
    if not np.isfinite(array[:, :3]).all():
        raise ValueError("rows must contain only finite z, Bx, and By values")
    if not np.all(np.diff(array[:, 0]) > 0):
        raise ValueError("z coordinates must be strictly increasing")
    return array


def _fit_multiple(
    z: np.ndarray,
    values: np.ndarray,
    period_mm: float,
    orders: tuple[int, ...],
) -> dict[int, HarmonicResult]:
    columns: list[np.ndarray] = []
    for order in orders:
        if not isinstance(order, int) or isinstance(order, bool) or order < 1:
            raise ValueError("harmonic order must be an integer >= 1")
        wave_number = 2 * np.pi * order / period_mm
        columns.extend((np.cos(wave_number * z), np.sin(wave_number * z)))
    columns.append(np.ones_like(z))
    matrix = np.column_stack(columns)
    coefficients, *_ = np.linalg.lstsq(matrix, values, rcond=None)
    offset = float(coefficients[-1])
    results: dict[int, HarmonicResult] = {}
    for index, order in enumerate(orders):
        cosine = float(coefficients[2 * index])
        sine = float(coefficients[2 * index + 1])
        wave_number = 2 * np.pi * order / period_mm
        component = cosine * np.cos(wave_number * z) + sine * np.sin(wave_number * z)
        # Per-order residual measures everything not explained by this order plus the DC offset.
        residual = values - (component + offset)
        results[order] = HarmonicResult(
            order=order,
            amplitude=float(math.hypot(cosine, sine)),
            phase_rad=float(math.atan2(sine, cosine)),
            offset=offset,
            rms_residual=float(np.sqrt(np.mean(residual ** 2))),
        )
    return results

def _central_data(
    rows: Iterable[Sequence[float]], period_mm: float, central_periods: int
) -> np.ndarray:
    period_mm = positive_number("period_mm", period_mm)
    if not isinstance(central_periods, int) or isinstance(central_periods, bool) or central_periods < 1:
        raise ValueError("central_periods must be an integer >= 1")
    array = _validated_array(rows)
    mask = np.abs(array[:, 0]) <= 0.5 * central_periods * period_mm
    if int(mask.sum()) < 8:
        raise ValueError("not enough samples in the requested central interval")
    return array[mask]


def analyze_harmonics(
    rows: Iterable[Sequence[float]],
    period_mm: float,
    harmonic_orders: Sequence[int] = (1, 3, 5),
    central_periods: int = 6,
) -> dict[int, tuple[HarmonicResult, HarmonicResult]]:
    data = _central_data(rows, period_mm, central_periods)
    orders = tuple(harmonic_orders)
    if not orders:
        raise ValueError("harmonic_orders must not be empty")
    if len(set(orders)) != len(orders):
        raise ValueError("harmonic_orders must not contain duplicates")
    bx_results = _fit_multiple(data[:, 0], data[:, 1], float(period_mm), orders)
    by_results = _fit_multiple(data[:, 0], data[:, 2], float(period_mm), orders)
    return {order: (bx_results[order], by_results[order]) for order in orders}


def analyze_fundamental(
    rows: Iterable[Sequence[float]], period_mm: float, central_periods: int = 6
) -> tuple[HarmonicResult, HarmonicResult]:
    return analyze_harmonics(rows, period_mm, (1,), central_periods)[1]


def polarization_metrics(bx: HarmonicResult, by: HarmonicResult) -> dict[str, float]:
    ax = finite_number("Bx amplitude", bx.amplitude)
    ay = finite_number("By amplitude", by.amplitude)
    if ax < 0 or ay < 0:
        raise ValueError("harmonic amplitudes must be nonnegative")
    relative_phase = (by.phase_rad - bx.phase_rad + math.pi) % (2 * math.pi) - math.pi
    denominator = ax * ax + ay * ay
    if denominator == 0:
        return {
            "ratio": math.nan,
            "phase_deg": math.degrees(relative_phase),
            "s1": math.nan,
            "s2": math.nan,
            "s3": math.nan,
            "circularity": math.nan,
        }
    s1 = float(np.clip((ax * ax - ay * ay) / denominator, -1.0, 1.0))
    s2 = float(np.clip(2 * ax * ay * math.cos(relative_phase) / denominator, -1.0, 1.0))
    s3 = float(np.clip(2 * ax * ay * math.sin(relative_phase) / denominator, -1.0, 1.0))
    return {
        "ratio": ax / ay if ay else math.inf,
        "phase_deg": math.degrees(relative_phase),
        "s1": s1,
        "s2": s2,
        "s3": s3,
        "circularity": float(np.clip(abs(s3), 0.0, 1.0)),
    }


def evaluate_quality(
    bx: HarmonicResult,
    by: HarmonicResult,
    criteria: QualityCriteria = QualityCriteria(),
) -> QualityAssessment:
    criteria.validate()
    total_amplitude = math.hypot(bx.amplitude, by.amplitude)
    combined_residual = math.hypot(bx.rms_residual, by.rms_residual)
    relative_residual = combined_residual / total_amplitude if total_amplitude else math.inf
    reasons: list[str] = []
    if total_amplitude < criteria.min_total_amplitude_t:
        reasons.append(
            f"total fundamental amplitude {total_amplitude:.6g} T is below "
            f"{criteria.min_total_amplitude_t:.6g} T"
        )
    if relative_residual > criteria.max_relative_residual:
        reasons.append(
            f"relative harmonic residual {relative_residual:.6g} exceeds "
            f"{criteria.max_relative_residual:.6g}"
        )
    return QualityAssessment(
        passed=not reasons,
        total_amplitude_t=total_amplitude,
        combined_rms_residual_t=combined_residual,
        relative_residual=relative_residual,
        reasons=tuple(reasons),
    )


def harmonics_to_dict(
    harmonics: dict[int, tuple[HarmonicResult, HarmonicResult]]
) -> dict[str, dict[str, dict[str, float | int]]]:
    return {
        str(order): {"Bx": asdict(pair[0]), "By": asdict(pair[1])}
        for order, pair in harmonics.items()
    }
