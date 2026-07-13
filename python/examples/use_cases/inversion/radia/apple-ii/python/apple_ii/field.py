from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Iterable, Sequence

import numpy as np

from ._util import finite_number, require_radia
from .geometry import Apple2Device, object_id


@dataclass(frozen=True)
class FieldIntegral:
    first_t_mm: tuple[float, float, float]
    second_t_mm2: tuple[float, float, float]
    z_start_mm: float
    z_end_mm: float


def _as_field_array(rows: Iterable[Sequence[float]]) -> np.ndarray:
    array = np.asarray(list(rows), dtype=float)
    if array.ndim != 2 or array.shape[1] < 4 or array.shape[0] < 2:
        raise ValueError("field rows must be an N x 4 array with N >= 2")
    array = array[:, :4]
    if not np.isfinite(array).all():
        raise ValueError("field rows must contain only finite values")
    if not np.all(np.diff(array[:, 0]) > 0):
        raise ValueError("z coordinates must be strictly increasing")
    return array


def sample_axis(
    device: Apple2Device | int,
    z_min_mm: float,
    z_max_mm: float,
    n_points: int,
    *,
    radia_module: ModuleType | object | None = None,
) -> list[tuple[float, float, float, float]]:
    if not isinstance(n_points, int) or isinstance(n_points, bool) or n_points < 2:
        raise ValueError("n_points must be an integer >= 2")
    z_min_mm = finite_number("z_min_mm", z_min_mm)
    z_max_mm = finite_number("z_max_mm", z_max_mm)
    if z_max_mm <= z_min_mm:
        raise ValueError("z_max_mm must be greater than z_min_mm")
    radia = require_radia(radia_module or getattr(device, "_radia", None))
    z_values = np.linspace(z_min_mm, z_max_mm, n_points)
    result: list[tuple[float, float, float, float]] = []
    for z_mm in z_values:
        field = radia.Fld(object_id(device), "b", [0.0, 0.0, float(z_mm)])
        if len(field) < 3:
            raise ValueError("RADIA Fld returned fewer than three components")
        row = (float(z_mm), float(field[0]), float(field[1]), float(field[2]))
        if not np.isfinite(row).all():
            raise ValueError("RADIA Fld returned a non-finite value")
        result.append(row)
    return result


def field_integrals(rows: Iterable[Sequence[float]]) -> FieldIntegral:
    """Numerically integrate B and its cumulative first integral over sampled z.

    Coordinates are in millimetres, so the returned units are T·mm and T·mm².
    The second integral uses zero cumulative first integral at ``z_start_mm``.
    """
    array = _as_field_array(rows)
    z = array[:, 0]
    fields = array[:, 1:4]
    dz = np.diff(z)
    segments = 0.5 * (fields[:-1] + fields[1:]) * dz[:, None]
    cumulative_first = np.vstack((np.zeros(3), np.cumsum(segments, axis=0)))
    first = cumulative_first[-1]
    second_segments = 0.5 * (cumulative_first[:-1] + cumulative_first[1:]) * dz[:, None]
    second = np.sum(second_segments, axis=0)
    return FieldIntegral(
        tuple(float(value) for value in first),
        tuple(float(value) for value in second),
        float(z[0]),
        float(z[-1]),
    )


def write_csv(rows: Iterable[Sequence[float]], path: str | Path) -> None:
    array = _as_field_array(rows)
    with Path(path).open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(("z_mm", "Bx_T", "By_T", "Bz_T"))
        writer.writerows(array.tolist())
