from __future__ import annotations

from dataclasses import asdict
from types import ModuleType
from typing import Sequence

import numpy as np

from .analysis import analyze_harmonics, harmonics_to_dict, polarization_metrics
from .field import field_integrals, sample_axis
from .geometry import Apple2Device


def row_field_diagnostics(
    device: Apple2Device,
    z_min_mm: float,
    z_max_mm: float,
    n_points: int,
    *,
    period_mm: float,
    harmonic_orders: Sequence[int],
    central_periods: int,
    combined_rows: Sequence[Sequence[float]] | None = None,
    radia_module: ModuleType | object | None = None,
) -> dict[str, object]:
    """Analyze each magnetic row independently and verify field superposition.

    The result is intended for geometry debugging. It does not establish that a
    row construction or magnetization convention is physically correct.
    """
    rows_report: dict[str, object] = {}
    sampled_arrays: list[np.ndarray] = []

    for name, row in device.rows.items():
        sampled = sample_axis(
            row.object_id,
            z_min_mm,
            z_max_mm,
            n_points,
            radia_module=radia_module or device._radia,
        )
        sampled_array = np.asarray(sampled, dtype=float)
        sampled_arrays.append(sampled_array)
        harmonics = analyze_harmonics(
            sampled,
            period_mm,
            harmonic_orders=harmonic_orders,
            central_periods=central_periods,
        )
        bx1, by1 = harmonics[1]
        rows_report[name] = {
            "object_id": row.object_id,
            "motion_coefficient": row.motion_coefficient,
            "current_shift_mm": row.current_shift_mm,
            "harmonics": harmonics_to_dict(harmonics),
            "fundamental_field_ellipse": polarization_metrics(bx1, by1),
            "field_integrals": asdict(field_integrals(sampled)),
        }

    superposition: dict[str, object] = {"checked": False}
    if combined_rows is not None:
        combined = np.asarray(list(combined_rows), dtype=float)
        stacked = np.stack(sampled_arrays, axis=0)
        summed = np.column_stack((stacked[0, :, 0], np.sum(stacked[:, :, 1:4], axis=0)))
        if combined.shape != summed.shape:
            raise ValueError("combined field rows do not match per-row sample shape")
        if not np.allclose(combined[:, 0], summed[:, 0], rtol=0.0, atol=1.0e-12):
            raise ValueError("combined and per-row samples use different z coordinates")
        delta = combined[:, 1:4] - summed[:, 1:4]
        max_abs = np.max(np.abs(delta), axis=0)
        rms = np.sqrt(np.mean(delta * delta, axis=0))
        superposition = {
            "checked": True,
            "max_abs_error_T": {
                "Bx": float(max_abs[0]),
                "By": float(max_abs[1]),
                "Bz": float(max_abs[2]),
            },
            "rms_error_T": {
                "Bx": float(rms[0]),
                "By": float(rms[1]),
                "Bz": float(rms[2]),
            },
        }

    return {
        "rows": rows_report,
        "superposition": superposition,
    }
