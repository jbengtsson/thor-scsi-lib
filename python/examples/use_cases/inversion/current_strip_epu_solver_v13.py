#!/usr/bin/env python3
"""Current-strip correction solver with a full 2D fitting objective.

This refactor preserves the V12 command-line interface and numerical model,
but separates configuration, geometry, solving, diagnostics, acceptance, and
output publication into small testable units.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy as np
from scipy.optimize import lsq_linear

MU0 = 4.0e-7 * np.pi
SPEED_OF_LIGHT_M_S = 299_792_458.0
ELECTRON_REST_ENERGY_EV = 510_998.95
RELEASE_ID = "CURRENT-STRIP-SOLVER-24-REFACTORED"
PLOT_ELEVATION_DEG = 20.0
PLOT_AZIMUTH_DEG = -45.0


@dataclass(frozen=True)
class KickMap:
    length_m: float
    x_m: np.ndarray
    y_m: np.ndarray
    theta_x_urad: np.ndarray
    theta_y_urad: np.ndarray
    b2_t2m: np.ndarray | None = None

@dataclass(frozen=True)
class Strip:
    """One longitudinal strip with a rectangular transverse cross-section."""

    x_m: float
    y_m: float
    size_x_m: float
    size_y_m: float
    length_m: float
    family: str
    sheet: str
    coordinate_m: float

@dataclass(frozen=True)
class StripGeometry:
    layout: str
    horizontal_count_per_plane: int
    vertical_count_per_plane: int
    horizontal_half_span_m: float
    vertical_half_span_m: float
    horizontal_pitch_m: float
    vertical_pitch_m: float
    top_bottom_clear_gap_m: float
    left_right_clear_gap_m: float
    side_gap_source: str

@dataclass(frozen=True)
class ZeroCut:
    values: np.ndarray
    method: str
    lower_coordinate_m: float
    upper_coordinate_m: float

@dataclass(frozen=True)
class CurrentSolution:
    currents_a: np.ndarray
    status: int
    message: str
    cost: float
    optimality: float
    iterations: int | None

@dataclass(frozen=True)
class AffineDecomposition:
    """Least-squares affine decomposition along one coordinate axis."""

    fitted: np.ndarray
    residual: np.ndarray
    offset: np.ndarray
    slope_per_m: np.ndarray

@dataclass(frozen=True)
class AffinePlaneDecomposition:
    """Least-squares affine-plane decomposition over transverse points."""

    fitted: np.ndarray
    residual: np.ndarray
    offset: np.ndarray
    slope_x_per_m: np.ndarray
    slope_y_per_m: np.ndarray

def _number_after(lines: list[str], label: str) -> float:
    pattern = re.compile(label, re.IGNORECASE)
    number = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?")

    for i, line in enumerate(lines):
        if not pattern.search(line):
            continue
        for j in range(i + 1, min(i + 4, len(lines))):
            match = number.search(lines[j])
            if match:
                return float(match.group(0))
        raise ValueError(f"Missing numeric value after header {label!r}")

    raise ValueError(f"Header {label!r} not found")

def _read_block(
    lines: list[str], start_index: int, nx: int, ny: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray([float(v) for v in lines[start_index + 1].split()])
    if x.size != nx:
        raise ValueError(f"x header contains {x.size} values; expected {nx}")

    y_values: list[float] = []
    rows: list[list[float]] = []

    for line in lines[start_index + 2 : start_index + 2 + ny]:
        values = [float(v) for v in line.split()]
        if len(values) != nx + 1:
            raise ValueError(
                f"Kick-map row contains {len(values)} values; expected {nx + 1}"
            )
        y_values.append(values[0])
        rows.append(values[1:])

    if len(rows) != ny:
        raise ValueError("Truncated kick-map block")

    return x, np.asarray(y_values), np.asarray(rows)

def parse_kick_units(value: str) -> str:
    """Return the canonical kick-unit name accepted by the solver.

    The ASCII spelling ``T^2m^2`` is canonical.  A few visually equivalent
    spellings are accepted so copied command lines using superscript 2 or a
    Unicode circumflex do not fail unexpectedly.
    """
    compact = value.strip().replace(" ", "")
    compact = compact.replace("²", "^2").replace("ˆ", "^")
    lowered = compact.lower()
    if lowered in {"microrad", "micro-rad", "urad", "µrad", "μrad"}:
        return "microrad"
    if lowered in {"t^2m^2", "t2m2", "t^2m2", "t2m^2"}:
        return "T^2m^2"
    raise argparse.ArgumentTypeError(
        "units must be 'microrad' or 'T^2m^2'"
    )

def electron_beam_rigidity(beam_energy_ev: float) -> tuple[float, float]:
    """Return electron momentum [eV/c] and |B rho| [T m].

    ``beam_energy_ev`` is the electron total energy.  The exact relativistic
    relation p*c = sqrt(E^2 - (m_e*c^2)^2) is used rather than the
    ultra-relativistic approximation p*c = E.
    """
    if not np.isfinite(beam_energy_ev):
        raise ValueError("beam-energy must be finite")
    if beam_energy_ev <= ELECTRON_REST_ENERGY_EV:
        raise ValueError(
            "beam-energy must exceed the electron rest energy "
            f"({ELECTRON_REST_ENERGY_EV:.8g} eV)"
        )
    momentum_ev_c = float(
        np.sqrt(
            (beam_energy_ev - ELECTRON_REST_ENERGY_EV)
            * (beam_energy_ev + ELECTRON_REST_ENERGY_EV)
        )
    )
    rigidity_tm = momentum_ev_c / SPEED_OF_LIGHT_M_S
    if not np.isfinite(rigidity_tm) or rigidity_tm <= 0.0:
        raise ValueError("computed beam rigidity is invalid")
    return momentum_ev_c, rigidity_tm

def kick_input_to_urad_scale(
    units: str,
    beam_rigidity_tm: float,
    kick_scale: float,
) -> float:
    """Return the multiplier converting input kick values to microradians."""
    if units == "microrad":
        unit_scale = 1.0
    elif units == "T^2m^2":
        # Second-order kick factor C [T^2 m^2]: theta [rad] = C/(B rho)^2.
        unit_scale = 1.0e6 / (beam_rigidity_tm * beam_rigidity_tm)
    else:
        raise ValueError(f"unsupported kick units: {units!r}")
    return kick_scale * unit_scale

def load_kick_map(path: Path) -> KickMap:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    nx = int(_number_after(lines, r"Number of Horizontal Points"))
    ny = int(_number_after(lines, r"Number of Vertical Points"))
    length_m = float(_number_after(lines, r"Undulator Length"))

    starts = [
        i for i, line in enumerate(lines)
        if line.strip().upper() == "START"
    ]
    if len(starts) < 2:
        raise ValueError("Horizontal and vertical kick blocks are required")

    x, y, theta_x = _read_block(lines, starts[0], nx, ny)
    x2, y2, theta_y = _read_block(lines, starts[1], nx, ny)

    if not (np.array_equal(x, x2) and np.array_equal(y, y2)):
        raise ValueError("Horizontal and vertical kick grids differ")

    b2 = None
    if len(starts) >= 3:
        x3, y3, b2 = _read_block(lines, starts[2], nx, ny)
        if not (np.array_equal(x, x3) and np.array_equal(y, y3)):
            raise ValueError("Integrated-B² grid differs from the kick grid")

    return KickMap(
        length_m=length_m,
        x_m=x,
        y_m=y,
        theta_x_urad=theta_x,
        theta_y_urad=theta_y,
        b2_t2m=b2,
    )

def zero_coordinate_cut(
    coordinates_m: np.ndarray,
    values: np.ndarray,
    axis: int,
    coordinate_name: str,
) -> ZeroCut:
    """Return ``values`` evaluated at coordinate zero along ``axis``.

    If zero is present on the grid within floating-point tolerance, the
    exact grid plane is returned. Otherwise, zero must be bracketed by two
    coordinates and a linear interpolation is performed. Extrapolation is
    deliberately rejected because it would silently change the physical
    diagnostic requested by the solver.
    """
    coordinates = np.asarray(coordinates_m, dtype=float)
    array = np.asarray(values, dtype=float)

    if coordinates.ndim != 1 or coordinates.size < 1:
        raise ValueError(f"{coordinate_name} coordinates must be one-dimensional")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"{coordinate_name} coordinates must be finite")
    if axis < 0:
        axis += array.ndim
    if axis < 0 or axis >= array.ndim:
        raise ValueError("axis is outside the values array")
    if array.shape[axis] != coordinates.size:
        raise ValueError(
            f"values axis {axis} has length {array.shape[axis]}; "
            f"expected {coordinates.size} from {coordinate_name} coordinates"
        )

    order = np.argsort(coordinates, kind="stable")
    sorted_coordinates = coordinates[order]
    sorted_values = np.take(array, order, axis=axis)

    if np.any(np.diff(sorted_coordinates) <= 0.0):
        raise ValueError(f"{coordinate_name} coordinates must be unique")

    scale = max(1.0, float(np.max(np.abs(sorted_coordinates))))
    tolerance = 32.0 * np.finfo(float).eps * scale
    exact_indices = np.flatnonzero(np.abs(sorted_coordinates) <= tolerance)
    if exact_indices.size:
        exact_index = int(
            exact_indices[
                np.argmin(np.abs(sorted_coordinates[exact_indices]))
            ]
        )
        coordinate = float(sorted_coordinates[exact_index])
        return ZeroCut(
            values=np.take(sorted_values, exact_index, axis=axis),
            method="exact",
            lower_coordinate_m=coordinate,
            upper_coordinate_m=coordinate,
        )

    if sorted_coordinates[0] > 0.0 or sorted_coordinates[-1] < 0.0:
        raise ValueError(
            f"{coordinate_name}=0 is outside the kick-map grid; "
            "zero-cut extrapolation is not allowed"
        )

    upper_index = int(np.searchsorted(sorted_coordinates, 0.0, side="right"))
    lower_index = upper_index - 1
    lower_coordinate = float(sorted_coordinates[lower_index])
    upper_coordinate = float(sorted_coordinates[upper_index])
    fraction = -lower_coordinate / (upper_coordinate - lower_coordinate)

    lower_values = np.take(sorted_values, lower_index, axis=axis)
    upper_values = np.take(sorted_values, upper_index, axis=axis)
    interpolated = (1.0 - fraction) * lower_values + fraction * upper_values

    return ZeroCut(
        values=interpolated,
        method="linear_interpolation",
        lower_coordinate_m=lower_coordinate,
        upper_coordinate_m=upper_coordinate,
    )

def _uniform_strip_centres(
    half_span_m: float,
    count: int,
    strip_size_along_array_m: float,
    label: str,
) -> tuple[np.ndarray, float]:
    """Return uniformly spaced centres and pitch for one strip array."""
    if count < 1:
        raise ValueError(f"{label} strip count must be positive")
    if not np.isfinite(half_span_m) or half_span_m <= 0.0:
        raise ValueError(f"{label} half span must be finite and positive")
    if (
        not np.isfinite(strip_size_along_array_m)
        or strip_size_along_array_m <= 0.0
    ):
        raise ValueError(f"{label} strip size must be finite and positive")

    if count == 1:
        return np.asarray([0.0], dtype=float), 0.0

    centres = np.linspace(-half_span_m, +half_span_m, count, dtype=float)
    pitch_m = float(2.0 * half_span_m / (count - 1))
    if strip_size_along_array_m > pitch_m + 1.0e-15:
        raise ValueError(
            f"{label} strip size exceeds the derived pitch; reduce the "
            "strip count or width, or increase the corresponding range"
        )
    return centres, pitch_m

def make_strip_arrays(
    layout: str,
    horizontal_half_span_m: float,
    horizontal_count_per_plane: int,
    vertical_half_span_m: float,
    vertical_count_per_plane: int,
    width_m: float,
    thickness_m: float,
    top_bottom_clear_gap_m: float,
    left_right_clear_gap_m: float,
    length_m: float,
    side_gap_source: str,
) -> tuple[list[Strip], StripGeometry]:
    """Build horizontal top/bottom and/or vertical left/right strip sheets.

    All conductors run longitudinally along +s.  For the horizontal family,
    strip centres vary in x and the broad cross-section dimension is x.  For
    the vertical family, strip centres vary in y and the same rectangular
    cross-section is rotated by 90 degrees, so its broad dimension is y.
    """
    if layout not in {"horizontal", "vertical", "both"}:
        raise ValueError(f"unknown strip layout: {layout!r}")
    if min(width_m, thickness_m, length_m) <= 0.0:
        raise ValueError("strip width, thickness, and length must be positive")
    if not all(np.isfinite(v) for v in (width_m, thickness_m, length_m)):
        raise ValueError("strip dimensions must be finite")

    strips: list[Strip] = []
    horizontal_pitch_m = float("nan")
    vertical_pitch_m = float("nan")
    horizontal_count = 0
    vertical_count = 0

    if layout in {"horizontal", "both"}:
        if (
            not np.isfinite(top_bottom_clear_gap_m)
            or top_bottom_clear_gap_m <= 0.0
        ):
            raise ValueError("top/bottom clear gap must be finite and positive")
        x_centres, horizontal_pitch_m = _uniform_strip_centres(
            horizontal_half_span_m,
            horizontal_count_per_plane,
            width_m,
            "horizontal-array",
        )
        horizontal_count = horizontal_count_per_plane
        y_centre = top_bottom_clear_gap_m / 2.0 + thickness_m / 2.0
        for sheet, sign in (("top", +1.0), ("bottom", -1.0)):
            for coordinate in x_centres:
                strips.append(
                    Strip(
                        x_m=float(coordinate),
                        y_m=float(sign * y_centre),
                        size_x_m=width_m,
                        size_y_m=thickness_m,
                        length_m=length_m,
                        family="horizontal",
                        sheet=sheet,
                        coordinate_m=float(coordinate),
                    )
                )

    if layout in {"vertical", "both"}:
        if (
            not np.isfinite(left_right_clear_gap_m)
            or left_right_clear_gap_m <= 0.0
        ):
            raise ValueError("left/right clear gap must be finite and positive")
        y_centres, vertical_pitch_m = _uniform_strip_centres(
            vertical_half_span_m,
            vertical_count_per_plane,
            width_m,
            "vertical-array",
        )
        vertical_count = vertical_count_per_plane
        x_centre = left_right_clear_gap_m / 2.0 + thickness_m / 2.0
        for sheet, sign in (("right", +1.0), ("left", -1.0)):
            for coordinate in y_centres:
                strips.append(
                    Strip(
                        x_m=float(sign * x_centre),
                        y_m=float(coordinate),
                        size_x_m=thickness_m,
                        size_y_m=width_m,
                        length_m=length_m,
                        family="vertical",
                        sheet=sheet,
                        coordinate_m=float(coordinate),
                    )
                )

    if not strips:
        raise ValueError("strip layout produced no conductors")

    geometry = StripGeometry(
        layout=layout,
        horizontal_count_per_plane=horizontal_count,
        vertical_count_per_plane=vertical_count,
        horizontal_half_span_m=(
            horizontal_half_span_m if horizontal_count else float("nan")
        ),
        vertical_half_span_m=(
            vertical_half_span_m if vertical_count else float("nan")
        ),
        horizontal_pitch_m=horizontal_pitch_m,
        vertical_pitch_m=vertical_pitch_m,
        top_bottom_clear_gap_m=(
            top_bottom_clear_gap_m if horizontal_count else float("nan")
        ),
        left_right_clear_gap_m=(
            left_right_clear_gap_m if vertical_count else float("nan")
        ),
        side_gap_source=side_gap_source if vertical_count else "not_applicable",
    )
    return strips, geometry

def response_matrix(
    strips: list[Strip],
    points_xy_m: np.ndarray,
    component: str,
    quadrature_order: int = 4,
) -> np.ndarray:
    """
    Return an integrated-field response matrix in T m / A.

    For a +s-directed filament current in right-handed x,y,s coordinates:
        B_x = -mu0 I (y-y0) / (2 pi r²)
        B_y = +mu0 I (x-x0) / (2 pi r²)

    The rectangular cross-section is integrated by normalized
    Gauss-Legendre quadrature.
    """
    points = np.asarray(points_xy_m, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points_xy_m must have shape (n_points, 2)")
    if component not in {"Ix", "Iy"}:
        raise ValueError("component must be 'Ix' or 'Iy'")
    if quadrature_order < 1:
        raise ValueError("quadrature_order must be positive")

    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    weights = weights / 2.0  # normalized average over [-1, 1]

    response = np.zeros((len(points), len(strips)), dtype=float)

    for j, strip in enumerate(strips):
        source_x = strip.x_m + 0.5 * strip.size_x_m * nodes
        source_y = strip.y_m + 0.5 * strip.size_y_m * nodes
        integrated_field = np.zeros(len(points), dtype=float)

        for ix, wx in enumerate(weights):
            dx = points[:, 0] - source_x[ix]
            for iy, wy in enumerate(weights):
                dy = points[:, 1] - source_y[iy]
                radius_sq = dx * dx + dy * dy

                if np.any(radius_sq == 0.0):
                    raise ValueError("Observation point lies inside a filament node")

                if component == "Ix":
                    field_per_amp = -MU0 * dy / (2.0 * np.pi * radius_sq)
                else:
                    field_per_amp = +MU0 * dx / (2.0 * np.pi * radius_sq)

                integrated_field += (
                    wx * wy * field_per_amp * strip.length_m
                )

        response[:, j] = integrated_field

    return response

def solve_currents(
    response: np.ndarray,
    target: np.ndarray,
    ridge: float,
    current_limit_a: float,
) -> CurrentSolution:
    if ridge < 0.0:
        raise ValueError("ridge must be non-negative")
    if current_limit_a <= 0.0:
        raise ValueError("current_limit_a must be positive")

    column_scale = float(np.median(np.linalg.norm(response, axis=0)))
    if not np.isfinite(column_scale) or column_scale <= 0.0:
        raise ValueError("Response matrix has zero or invalid column norm")

    if ridge > 0.0:
        augmented_response = np.vstack(
            [
                response,
                np.sqrt(ridge)
                * column_scale
                * np.eye(response.shape[1]),
            ]
        )
        augmented_target = np.concatenate(
            [target, np.zeros(response.shape[1])]
        )
    else:
        augmented_response = response
        augmented_target = target

    result = lsq_linear(
        augmented_response,
        augmented_target,
        bounds=(-current_limit_a, current_limit_a),
        lsmr_tol="auto",
        max_iter=1000,
    )

    if not result.success:
        raise RuntimeError(f"Bounded least-squares solve failed: {result.message}")

    return CurrentSolution(
        currents_a=result.x,
        status=int(result.status),
        message=str(result.message).replace("\n", " ").strip(),
        cost=float(result.cost),
        optimality=float(result.optimality),
        iterations=(None if result.nit is None else int(result.nit)),
    )

def affine_decomposition(
    coordinate_m: np.ndarray,
    values: np.ndarray,
) -> AffineDecomposition:
    """Split values into least-squares affine and nonlinear parts.

    The first axis of ``values`` must correspond to ``coordinate_m``.  For a
    vector, this fits ``a + b*q``.  For a matrix, the same projection is
    applied independently to every column.  Removing the affine part from
    both the target and each response column makes the nonlinear objective
    invariant to constant steering and linear-gradient components.
    """
    coordinate = np.asarray(coordinate_m, dtype=float)
    array = np.asarray(values, dtype=float)

    if coordinate.ndim != 1 or coordinate.size < 2:
        raise ValueError("affine projection requires at least two coordinates")
    if not np.all(np.isfinite(coordinate)) or not np.all(np.isfinite(array)):
        raise ValueError("affine projection inputs must be finite")
    if array.shape[0] != coordinate.size:
        raise ValueError(
            "the first values axis must match the affine-fit coordinate"
        )
    if np.ptp(coordinate) <= 0.0:
        raise ValueError("affine-fit coordinates must span a nonzero interval")

    coordinate_scale_m = float(np.max(np.abs(coordinate)))
    if coordinate_scale_m == 0.0:
        coordinate_scale_m = float(np.ptp(coordinate))
    normalized_coordinate = coordinate / coordinate_scale_m
    design = np.column_stack(
        [np.ones(coordinate.size, dtype=float), normalized_coordinate]
    )

    original_shape = array.shape
    flattened = array.reshape(coordinate.size, -1)
    coefficients, _, _, _ = np.linalg.lstsq(design, flattened, rcond=None)
    fitted_flat = design @ coefficients
    residual_flat = flattened - fitted_flat

    # An exactly affine target should become an exact numerical zero rather
    # than a machine-epsilon residual that can make a relative error ill posed.
    column_scale = np.max(np.abs(flattened), axis=0)
    tolerance = (
        256.0
        * np.finfo(float).eps
        * np.maximum(column_scale, np.finfo(float).tiny)
    )
    residual_flat[np.abs(residual_flat) <= tolerance] = 0.0

    trailing_shape = original_shape[1:]
    return AffineDecomposition(
        fitted=fitted_flat.reshape(original_shape),
        residual=residual_flat.reshape(original_shape),
        offset=coefficients[0].reshape(trailing_shape),
        slope_per_m=(coefficients[1] / coordinate_scale_m).reshape(
            trailing_shape
        ),
    )

def affine_plane_decomposition(
    points_xy_m: np.ndarray,
    values: np.ndarray,
) -> AffinePlaneDecomposition:
    """Split values into a 2D affine plane and its nonlinear residual.

    The first axis of ``values`` must match ``points_xy_m``.  For a vector,
    this fits ``a + b*x + c*y``.  For a matrix, the same orthogonal projection
    is applied independently to every column.  Coordinate normalization keeps
    the least-squares design well conditioned while the returned slopes are
    expressed per metre.
    """
    points = np.asarray(points_xy_m, dtype=float)
    array = np.asarray(values, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or points.shape[0] < 3:
        raise ValueError("affine-plane projection requires at least three 2D points")
    if array.shape[0] != points.shape[0]:
        raise ValueError("the first values axis must match the affine-plane points")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(array)):
        raise ValueError("affine-plane projection inputs must be finite")

    x = points[:, 0]
    y = points[:, 1]
    x_scale_m = max(float(np.max(np.abs(x))), float(np.ptp(x)))
    y_scale_m = max(float(np.max(np.abs(y))), float(np.ptp(y)))
    if x_scale_m <= 0.0 or y_scale_m <= 0.0:
        raise ValueError("affine-plane points must span both x and y")

    design = np.column_stack(
        [
            np.ones(points.shape[0], dtype=float),
            x / x_scale_m,
            y / y_scale_m,
        ]
    )
    original_shape = array.shape
    flattened = array.reshape(points.shape[0], -1)
    coefficients, _, rank, _ = np.linalg.lstsq(design, flattened, rcond=None)
    if rank < 3:
        raise ValueError("affine-plane design is rank deficient")
    fitted_flat = design @ coefficients
    residual_flat = flattened - fitted_flat

    column_scale = np.max(np.abs(flattened), axis=0)
    tolerance = (
        256.0
        * np.finfo(float).eps
        * np.maximum(column_scale, np.finfo(float).tiny)
    )
    residual_flat[np.abs(residual_flat) <= tolerance] = 0.0

    trailing_shape = original_shape[1:]
    return AffinePlaneDecomposition(
        fitted=fitted_flat.reshape(original_shape),
        residual=residual_flat.reshape(original_shape),
        offset=coefficients[0].reshape(trailing_shape),
        slope_x_per_m=(coefficients[1] / x_scale_m).reshape(trailing_shape),
        slope_y_per_m=(coefficients[2] / y_scale_m).reshape(trailing_shape),
    )

def root_mean_square(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError("RMS input must contain finite values")
    return float(np.sqrt(np.mean(array * array)))

def safe_reduction_factor(raw_value: float, corrected_value: float) -> tuple[float, str]:
    """Return raw/corrected without raising on an exact-zero denominator."""
    if not np.isfinite(raw_value) or raw_value < 0.0:
        raise ValueError("raw reduction-factor value must be finite and non-negative")
    if not np.isfinite(corrected_value) or corrected_value < 0.0:
        raise ValueError(
            "corrected reduction-factor value must be finite and non-negative"
        )
    if corrected_value == 0.0:
        if raw_value == 0.0:
            return 1.0, "both_zero"
        return float("inf"), "corrected_zero"
    return raw_value / corrected_value, "finite"

def relative_rms_residual(residual: np.ndarray, target: np.ndarray) -> float:
    residual_rms = root_mean_square(residual)
    target_rms = root_mean_square(target)
    if target_rms == 0.0:
        return 0.0 if residual_rms == 0.0 else float("inf")
    return residual_rms / target_rms

def publish_staged_outputs(
    staged_paths: list[Path],
    final_paths: list[Path],
) -> None:
    """Publish completed files by atomic same-filesystem replacement."""
    if len(staged_paths) != len(final_paths):
        raise ValueError("staged and final output lists differ in length")
    missing = [path for path in staged_paths if not path.is_file()]
    if missing:
        raise RuntimeError(
            "Output staging is incomplete: " + ", ".join(str(path) for path in missing)
        )
    for staged_path, final_path in zip(staged_paths, final_paths):
        os.replace(staged_path, final_path)

def write_corrected_kick_map(
    path: Path,
    kick_map: KickMap,
    theta_x_corrected_urad: np.ndarray,
    theta_y_corrected_urad: np.ndarray,
    *,
    source_units: str,
    beam_energy_ev: float,
    beam_rigidity_tm: float,
    kick_scale: float,
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# Corrected current-strip kick map\n")
        stream.write("# Right-handed x,y,s convention; beam along +s\n")
        stream.write(f"# Source kick units: {source_units}\n")
        stream.write(f"# Electron total beam energy [eV]: {beam_energy_ev:.12g}\n")
        stream.write(f"# Computed beam rigidity [T m]: {beam_rigidity_tm:.12g}\n")
        stream.write(f"# Additional kick scale: {kick_scale:.12g}\n")
        stream.write("# Corrected kick blocks are written in micro-rad\n")
        stream.write("# Undulator Length [m]\n")
        stream.write(f"{kick_map.length_m:g}\n")
        stream.write("# Number of Horizontal Points\n")
        stream.write(f"{kick_map.x_m.size}\n")
        stream.write("# Number of Vertical Points\n")
        stream.write(f"{kick_map.y_m.size}\n")

        def write_block(label: str, values: np.ndarray) -> None:
            stream.write(f"# {label}\n")
            stream.write("START\n")
            stream.write(
                "\t" + "\t".join(f"{value:.10g}" for value in kick_map.x_m)
                + "\n"
            )
            for y_value, row in zip(kick_map.y_m, values):
                stream.write(
                    f"{y_value:.10g}\t"
                    + "\t".join(f"{value:.10g}" for value in row)
                    + "\n"
                )

        write_block("Horizontal Kick [micro-rad]", theta_x_corrected_urad)
        write_block("Vertical Kick [micro-rad]", theta_y_corrected_urad)

        if kick_map.b2_t2m is not None:
            write_block(
                "Longitudinally Integrated Squared Transverse "
                "Magnetic Field [T2m]",
                kick_map.b2_t2m,
            )

def crop_horizontal_range(
    x_m: np.ndarray,
    values: np.ndarray,
    x_range_m: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Crop arrays to |x| <= x_range_m, with the limit expressed in metres.

    The solver and exported kick map always retain the complete input grid.
    """
    if x_range_m is None:
        return x_m, values
    if x_range_m <= 0.0:
        raise ValueError("x_range_m must be positive")

    limit_m = x_range_m
    mask = np.abs(x_m) <= limit_m + 1.0e-15
    if not np.any(mask):
        raise ValueError(
            f"No kick-map x coordinates lie within ±{x_range_m:g} m"
        )
    return x_m[mask], values[..., mask]

def crop_vertical_range(
    y_m: np.ndarray,
    values: np.ndarray,
    y_range_m: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Crop arrays to |y| <= y_range_m, with the limit expressed in metres.

    The solver and exported kick map always retain the complete input grid.
    The first axis of values is assumed to correspond to y.
    """
    if y_range_m is None:
        return y_m, values
    if y_range_m <= 0.0:
        raise ValueError("y_range_m must be positive")

    limit_m = y_range_m
    mask = np.abs(y_m) <= limit_m + 1.0e-15
    if not np.any(mask):
        raise ValueError(
            f"No kick-map y coordinates lie within ±{y_range_m:g} m"
        )
    return y_m[mask], values[mask, ...]

def canonicalize_plot_grid(
    x_m: np.ndarray,
    y_m: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an ascending x/y grid with values reordered consistently.

    RADIA kick maps commonly list y from positive to negative values.  A
    canonical ascending grid avoids reverse polygon ordering in mplot3d and
    keeps the rendered surface orientation consistent with the axis frame.
    """
    x = np.asarray(x_m, dtype=float)
    y = np.asarray(y_m, dtype=float)
    z = np.asarray(values, dtype=float)

    if z.shape != (y.size, x.size):
        raise ValueError(
            "values must have shape (len(y_m), len(x_m)); "
            f"got {z.shape}, expected {(y.size, x.size)}"
        )

    x_order = np.argsort(x, kind="stable")
    y_order = np.argsort(y, kind="stable")

    x_sorted = x[x_order]
    y_sorted = y[y_order]
    z_sorted = z[np.ix_(y_order, x_order)]

    if np.any(np.diff(x_sorted) <= 0.0):
        raise ValueError("x coordinates must be strictly increasing")
    if np.any(np.diff(y_sorted) <= 0.0):
        raise ValueError("y coordinates must be strictly increasing")

    return x_sorted, y_sorted, z_sorted

def _make_2d_figure(keep_open: bool = False):
    """Create a fully fresh 2D figure.

    In batch mode we avoid pyplot figure state entirely by constructing a
    standalone Figure with its own Agg canvas.  This guarantees that the
    next plot starts from a clean graphics state.
    """
    if keep_open:
        figure = plt.figure(figsize=(8, 5))
        axes = figure.add_subplot(111)
        return figure, axes

    figure = Figure(figsize=(8, 5))
    FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    return figure, axes

def _make_3d_pair_figure(keep_open: bool = False):
    """Create a fresh paired 3D figure with two independent panels."""
    if keep_open:
        figure = plt.figure(figsize=(13.5, 5.6))
    else:
        figure = Figure(figsize=(13.5, 5.6))
        FigureCanvasAgg(figure)

    theta_x_axes = figure.add_axes(
        (0.035, 0.10, 0.39, 0.82), projection="3d"
    )
    theta_x_color_axes = figure.add_axes((0.435, 0.24, 0.015, 0.52))

    theta_y_axes = figure.add_axes(
        (0.535, 0.10, 0.39, 0.82), projection="3d"
    )
    theta_y_color_axes = figure.add_axes((0.935, 0.24, 0.015, 0.52))

    return (
        figure,
        theta_x_axes,
        theta_x_color_axes,
        theta_y_axes,
        theta_y_color_axes,
    )

def plot_cut(
    coordinate_mm: np.ndarray,
    raw_urad: np.ndarray,
    corrected_urad: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
    output_path: Path,
) -> None:
    figure, axes = _make_2d_figure()
    axes.plot(coordinate_mm, raw_urad, label="raw")
    axes.plot(coordinate_mm, corrected_urad, label="corrected")
    axes.set_title(title)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.grid(alpha=0.3)
    axes.legend()
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    try:
        plt.close(figure)
    except Exception:
        pass

def write_strip_currents_csv(
    path: Path,
    strips: list[Strip],
    currents_a: np.ndarray,
) -> None:
    """Write currents together with complete physical strip metadata."""
    currents = np.asarray(currents_a, dtype=float)
    if currents.shape != (len(strips),):
        raise ValueError("Current vector does not match strip geometry")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "index",
                "current_A",
                "family",
                "sheet",
                "x_m",
                "y_m",
                "size_x_m",
                "size_y_m",
                "length_m",
            ]
        )
        for index, (strip, current) in enumerate(zip(strips, currents)):
            writer.writerow(
                [
                    index,
                    f"{current:.10e}",
                    strip.family,
                    strip.sheet,
                    f"{strip.x_m:.10e}",
                    f"{strip.y_m:.10e}",
                    f"{strip.size_x_m:.10e}",
                    f"{strip.size_y_m:.10e}",
                    f"{strip.length_m:.10e}",
                ]
            )

def _representative_plane_data(
    strips: list[Strip],
    currents_a: np.ndarray,
    family: str,
    positive_sheet: str,
    negative_sheet: str,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    currents = np.asarray(currents_a, dtype=float)
    positive_indices = [
        i
        for i, strip in enumerate(strips)
        if strip.family == family and strip.sheet == positive_sheet
    ]
    negative_indices = [
        i
        for i, strip in enumerate(strips)
        if strip.family == family and strip.sheet == negative_sheet
    ]
    if not positive_indices or len(positive_indices) != len(negative_indices):
        raise ValueError(f"Incomplete symmetric {family} strip sheets")

    positive_order = sorted(
        positive_indices, key=lambda i: strips[i].coordinate_m
    )
    negative_order = sorted(
        negative_indices, key=lambda i: strips[i].coordinate_m
    )
    positive_coordinate = np.asarray(
        [strips[i].coordinate_m for i in positive_order], dtype=float
    )
    negative_coordinate = np.asarray(
        [strips[i].coordinate_m for i in negative_order], dtype=float
    )
    if not np.allclose(
        positive_coordinate, negative_coordinate, rtol=0.0, atol=1.0e-15
    ):
        raise ValueError(f"The two {family} sheets use different grids")

    positive_current = currents[positive_order]
    negative_current = currents[negative_order]
    mismatch = positive_current - negative_current
    return (
        positive_coordinate,
        positive_current,
        float(np.max(np.abs(mismatch))),
        float(np.sqrt(np.mean(mismatch * mismatch))),
    )

def plot_strip_currents_representative_planes(
    strips: list[Strip],
    currents_a: np.ndarray,
    output_path: Path,
    horizontal_range_mm: float,
    vertical_range_mm: float,
    keep_open: bool = False,
):
    """Plot one representative plane for each enabled array family."""
    families = {strip.family for strip in strips}
    panel_specs: list[tuple[str, str, str, str, float]] = []
    if "horizontal" in families:
        panel_specs.append(
            ("horizontal", "top", "bottom", "strip-centre x [mm]", horizontal_range_mm)
        )
    if "vertical" in families:
        panel_specs.append(
            ("vertical", "right", "left", "strip-centre y [mm]", vertical_range_mm)
        )
    if not panel_specs:
        raise ValueError("No strip-current family is available for plotting")

    if keep_open:
        figure = plt.figure(figsize=(8.0 if len(panel_specs) == 1 else 13.5, 5.0))
    else:
        figure = Figure(figsize=(8.0 if len(panel_specs) == 1 else 13.5, 5.0))
        FigureCanvasAgg(figure)

    symmetry: dict[str, tuple[float, float]] = {}
    for panel_index, (family, positive, negative, xlabel, limit_mm) in enumerate(
        panel_specs, start=1
    ):
        axes = figure.add_subplot(1, len(panel_specs), panel_index)
        coordinate_m, plotted_current, max_mismatch, rms_mismatch = (
            _representative_plane_data(
                strips, currents_a, family, positive, negative
            )
        )
        axes.plot(
            coordinate_m * 1.0e3,
            plotted_current,
            marker="o",
            markersize=3.0,
            linewidth=1.2,
        )
        axes.axhline(0.0, linewidth=0.8)
        axes.set_title(
            f"{family.capitalize()} array — {positive} plane"
        )
        axes.set_xlabel(xlabel)
        axes.set_ylabel("current [A]")
        axes.set_xlim(-limit_mm, limit_mm)
        axes.grid(alpha=0.3)
        symmetry[family] = (max_mismatch, rms_mismatch)

    figure.suptitle("Fitted strip currents — representative planes")
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)

    if keep_open:
        return figure, symmetry
    try:
        plt.close(figure)
    except Exception:
        pass
    return None, symmetry

def colour_and_z_limits(
    values: np.ndarray,
    shared_limit: float,
    mode: str,
):
    """
    Return a Matplotlib normalizer and z-axis limits.

    auto:
        Colour limits follow the finite displayed data.  If the data span
        zero, TwoSlopeNorm keeps zero at the neutral centre of the diverging
        colour map.  The z axis follows the same data with 5% visual padding.
    shared:
        Use the earlier symmetric component-wise limit, shared between raw
        and corrected figures.
    """
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("Cannot autoscale a plot with no finite values")

    if mode == "shared":
        limit = float(shared_limit)
        if not np.isfinite(limit) or limit <= 0.0:
            limit = max(float(np.max(np.abs(finite))), np.finfo(float).eps)
        normalizer = matplotlib.colors.Normalize(vmin=-limit, vmax=limit)
        return normalizer, -limit, limit

    if mode != "auto":
        raise ValueError(f"Unknown colour-scale mode: {mode!r}")

    data_min = float(np.min(finite))
    data_max = float(np.max(finite))

    if data_min == data_max:
        reference = max(abs(data_min), 1.0e-12)
        colour_pad = 0.05 * reference
        vmin = data_min - colour_pad
        vmax = data_max + colour_pad
    else:
        vmin = data_min
        vmax = data_max

    if vmin < 0.0 < vmax:
        normalizer = matplotlib.colors.TwoSlopeNorm(
            vmin=vmin,
            vcenter=0.0,
            vmax=vmax,
        )
    else:
        normalizer = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    span = vmax - vmin
    z_pad = 0.05 * span if span > 0.0 else max(abs(vmin), 1.0) * 0.05
    return normalizer, vmin - z_pad, vmax + z_pad

def _prepare_3d_component(
    x_m: np.ndarray,
    y_m: np.ndarray,
    values_urad: np.ndarray,
    x_range_m: float | None,
    y_range_m: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Crop, sort, and convert one kick component for 3D plotting."""
    x_plot_m, values_plot_urad = crop_horizontal_range(
        x_m, values_urad, x_range_m
    )
    y_plot_m, values_plot_urad = crop_vertical_range(
        y_m, values_plot_urad, y_range_m
    )
    x_plot_m, y_plot_m, values_plot_urad = canonicalize_plot_grid(
        x_plot_m,
        y_plot_m,
        values_plot_urad,
    )
    x_grid_mm, y_grid_mm = np.meshgrid(
        x_plot_m * 1.0e3,
        y_plot_m * 1.0e3,
        indexing="xy",
    )
    return x_grid_mm, y_grid_mm, values_plot_urad * 1.0e-3

def _draw_3d_component(
    figure,
    axes,
    color_axes,
    x_grid_mm: np.ndarray,
    y_grid_mm: np.ndarray,
    values_mrad: np.ndarray,
    title: str,
    shared_limit_urad: float,
    color_scale: str,
    elevation_deg: float,
    azimuth_deg: float,
    projection: str,
    box_aspect_x: float,
    box_aspect_y: float,
    box_aspect_z: float,
) -> None:
    """Draw one kick component into an already-created 3D panel."""
    normalizer, z_min_mrad, z_max_mrad = colour_and_z_limits(
        values_mrad,
        shared_limit=shared_limit_urad * 1.0e-3,
        mode=color_scale,
    )

    surface = axes.plot_surface(
        x_grid_mm,
        y_grid_mm,
        values_mrad,
        cmap="turbo",
        norm=normalizer,
        linewidth=0.0,
        antialiased=True,
        rcount=values_mrad.shape[0],
        ccount=values_mrad.shape[1],
    )

    axes.set_title(title, pad=2.0, fontsize=13)
    axes.set_xlabel("x [mm]", labelpad=3.0)
    axes.set_ylabel("y [mm]", labelpad=3.0)
    axes.set_zlabel("")
    axes.set_xlim(float(np.min(x_grid_mm)), float(np.max(x_grid_mm)))
    axes.set_ylim(float(np.min(y_grid_mm)), float(np.max(y_grid_mm)))
    axes.set_zlim(z_min_mrad, z_max_mrad)
    axes.set_proj_type(projection)
    axes.view_init(
        elev=elevation_deg,
        azim=azimuth_deg,
        roll=0.0,
        vertical_axis="z",
    )
    axes.set_box_aspect((box_aspect_x, box_aspect_y, box_aspect_z))
    axes.tick_params(axis="both", labelsize=8, pad=0)

    for axis in (axes.xaxis, axes.yaxis, axes.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.pane.set_edgecolor((0.55, 0.55, 0.55, 0.55))

    colorbar = figure.colorbar(surface, cax=color_axes)
    colorbar.ax.tick_params(labelsize=8, pad=2)

def plot_kickmaps_3d_pair(
    x_m: np.ndarray,
    y_m: np.ndarray,
    theta_x_urad: np.ndarray,
    theta_y_urad: np.ndarray,
    output_path: Path,
    theta_x_shared_limit_urad: float,
    theta_y_shared_limit_urad: float,
    figure_title: str,
    x_range_m: float | None = None,
    y_range_m: float | None = None,
    color_scale: str = "auto",
    elevation_deg: float = 20.0,
    azimuth_deg: float = -45.0,
    projection: str = "ortho",
    box_aspect_x: float = 2.5,
    box_aspect_y: float = 1.8,
    box_aspect_z: float = 0.85,
    keep_open: bool = False,
):
    """Save paired horizontal and vertical 3D kick-map panels."""
    theta_x_grid = _prepare_3d_component(
        x_m,
        y_m,
        theta_x_urad,
        x_range_m,
        y_range_m,
    )
    theta_y_grid = _prepare_3d_component(
        x_m,
        y_m,
        theta_y_urad,
        x_range_m,
        y_range_m,
    )

    (
        figure,
        theta_x_axes,
        theta_x_color_axes,
        theta_y_axes,
        theta_y_color_axes,
    ) = _make_3d_pair_figure(keep_open=keep_open)

    _draw_3d_component(
        figure,
        theta_x_axes,
        theta_x_color_axes,
        *theta_x_grid,
        title=r"$\theta_x$ [mrad]",
        shared_limit_urad=theta_x_shared_limit_urad,
        color_scale=color_scale,
        elevation_deg=elevation_deg,
        azimuth_deg=azimuth_deg,
        projection=projection,
        box_aspect_x=box_aspect_x,
        box_aspect_y=box_aspect_y,
        box_aspect_z=box_aspect_z,
    )
    _draw_3d_component(
        figure,
        theta_y_axes,
        theta_y_color_axes,
        *theta_y_grid,
        title=r"$\theta_y$ [mrad]",
        shared_limit_urad=theta_y_shared_limit_urad,
        color_scale=color_scale,
        elevation_deg=elevation_deg,
        azimuth_deg=azimuth_deg,
        projection=projection,
        box_aspect_x=box_aspect_x,
        box_aspect_y=box_aspect_y,
        box_aspect_z=box_aspect_z,
    )

    figure.suptitle(figure_title, y=0.985, fontsize=13)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")

    if keep_open:
        return figure

    try:
        plt.close(figure)
    except Exception:
        pass
    return None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fit longitudinal current strips to a 2D EPU kick map"
    )

    io = parser.add_argument_group("input and output")
    io.add_argument("--kick", type=Path, required=True)
    io.add_argument("--out-prefix", type=Path, default=Path("epu57_fixed"))
    io.add_argument(
        "--units", type=parse_kick_units, choices=("microrad", "T^2m^2"),
        default="microrad",
        help="input kick units: microrad (default) or T^2m^2",
    )
    io.add_argument(
        "--beam-energy", "--beam_energy", "-beam_energy",
        dest="beam_energy_ev", type=float, required=True, metavar="EV",
        help="electron total energy in eV",
    )
    io.add_argument(
        "--kick-scale", type=float, default=1.0, metavar="FACTOR",
        help="additional dimensionless kick multiplier (default: 1)",
    )

    objective = parser.add_argument_group("fit objective")
    objective.add_argument("--horizontal-weight", type=float, default=1.0)
    objective.add_argument("--vertical-weight", type=float, default=0.0)
    objective.add_argument(
        "--target-mode", choices=("nonlinear", "full"), default="nonlinear",
        help="nonlinear removes a+b*x+c*y; full fits the complete 2D maps",
    )
    objective.add_argument("--ridge", type=float, default=1.0e-5)
    objective.add_argument("--limit", type=float, default=6.0)
    objective.add_argument("--quadrature-order", type=int, default=4)

    geometry = parser.add_argument_group("strip geometry")
    geometry.add_argument(
        "--strip-layout", choices=("horizontal", "vertical", "both"),
        default="horizontal",
    )
    geometry.add_argument("--gap", type=float, default=0.016)
    geometry.add_argument("--side-gap", type=float, default=None, metavar="M")
    geometry.add_argument("--length", type=float, default=None)
    geometry.add_argument(
        "--strips-per-plane", "--horizontal-strips-per-plane",
        dest="strips_per_plane", type=int, default=21, metavar="N",
    )
    geometry.add_argument(
        "--horizontal-strip-range-mm", type=float, default=20.0, metavar="MM"
    )
    geometry.add_argument(
        "--vertical-strips-per-plane", type=int, default=None, metavar="N"
    )
    geometry.add_argument(
        "--vertical-strip-range-mm", type=float, default=None, metavar="MM"
    )
    geometry.add_argument("--width", type=float, default=2.0e-3)
    geometry.add_argument("--thickness", type=float, default=0.3e-3)

    aperture = parser.add_argument_group("model aperture")
    aperture.add_argument("--x-range-m", type=float, default=0.020, metavar="M")
    aperture.add_argument("--y-range-m", type=float, default=None, metavar="M")

    acceptance = parser.add_argument_group("acceptance")
    acceptance.add_argument(
        "--max-fit-relative-rms", type=float, default=None, metavar="FRACTION"
    )
    acceptance.add_argument(
        "--min-horizontal-rms-reduction-factor", type=float, default=1.0
    )
    acceptance.add_argument(
        "--min-vertical-rms-reduction-factor", type=float, default=1.0
    )
    acceptance.add_argument(
        "--min-full-map-rms-reduction-factor", type=float, default=1.0
    )
    acceptance.add_argument(
        "--min-full-map-peak-reduction-factor", type=float, default=0.0
    )
    acceptance.add_argument(
        "--current-bound-tolerance", type=float, default=1.0e-6
    )

    plot = parser.add_argument_group("plotting")
    plot.add_argument("--color-scale", choices=("auto", "shared"), default="auto")
    plot.add_argument("--elevation-deg", type=float, default=PLOT_ELEVATION_DEG)
    plot.add_argument("--azimuth-deg", type=float, default=PLOT_AZIMUTH_DEG)
    plot.add_argument("--projection", choices=("ortho", "persp"), default="ortho")
    plot.add_argument("--box-aspect-x", type=float, default=2.5)
    plot.add_argument("--box-aspect-y", type=float, default=1.8)
    plot.add_argument("--box-aspect-z", type=float, default=0.85)
    plot.add_argument("--show", action="store_true")
    return parser


def _require_finite(
    name: str,
    value: float | None,
    *,
    positive: bool = False,
    nonnegative: bool = False,
    allow_none: bool = False,
) -> None:
    if value is None:
        if allow_none:
            return
        raise ValueError(f"{name} is required")
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if positive and value <= 0.0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and value < 0.0:
        raise ValueError(f"{name} must be non-negative")


def validate_arguments(args: argparse.Namespace) -> argparse.Namespace:
    for name in ("kick_scale", "x_range_m", "gap", "width", "thickness", "limit"):
        _require_finite(name.replace("_", "-"), getattr(args, name), positive=True)
    _require_finite("y-range-m", args.y_range_m, positive=True, allow_none=True)
    _require_finite("side-gap", args.side_gap, positive=True, allow_none=True)
    _require_finite("length", args.length, positive=True, allow_none=True)
    _require_finite("ridge", args.ridge, nonnegative=True)
    _require_finite("horizontal-weight", args.horizontal_weight, nonnegative=True)
    _require_finite("vertical-weight", args.vertical_weight, nonnegative=True)
    if args.horizontal_weight == args.vertical_weight == 0.0:
        raise ValueError("at least one kick weight must be positive")
    if args.strips_per_plane < 1 or args.quadrature_order < 1:
        raise ValueError("strip counts and quadrature order must be positive")

    args.effective_vertical_strips_per_plane = (
        args.strips_per_plane
        if args.vertical_strips_per_plane is None
        else args.vertical_strips_per_plane
    )
    if args.effective_vertical_strips_per_plane < 1:
        raise ValueError("vertical-strips-per-plane must be positive")

    _require_finite(
        "horizontal-strip-range-mm", args.horizontal_strip_range_mm, positive=True
    )
    args.effective_horizontal_strip_range_mm = args.horizontal_strip_range_mm
    args.effective_vertical_strip_range_mm = (
        args.horizontal_strip_range_mm
        if args.vertical_strip_range_mm is None
        else args.vertical_strip_range_mm
    )
    _require_finite(
        "vertical-strip-range-mm", args.effective_vertical_strip_range_mm,
        positive=True,
    )

    _require_finite(
        "max-fit-relative-rms", args.max_fit_relative_rms,
        nonnegative=True, allow_none=True,
    )
    if args.max_fit_relative_rms is None:
        args.effective_max_fit_relative_rms = (
            0.15 if args.target_mode == "nonlinear" else 0.05
        )
        args.max_fit_relative_rms_source = "target_mode_default"
    else:
        args.effective_max_fit_relative_rms = args.max_fit_relative_rms
        args.max_fit_relative_rms_source = "user"

    for name in (
        "min_horizontal_rms_reduction_factor",
        "min_vertical_rms_reduction_factor",
        "min_full_map_rms_reduction_factor",
    ):
        _require_finite(name.replace("_", "-"), getattr(args, name), positive=True)
    _require_finite(
        "min-full-map-peak-reduction-factor",
        args.min_full_map_peak_reduction_factor,
        nonnegative=True,
    )
    _require_finite(
        "current-bound-tolerance", args.current_bound_tolerance,
        nonnegative=True,
    )
    if args.current_bound_tolerance >= 1.0:
        raise ValueError("current-bound-tolerance must lie in [0, 1)")

    for name in ("elevation_deg", "azimuth_deg"):
        _require_finite(name.replace("_", "-"), getattr(args, name))
    for name in ("box_aspect_x", "box_aspect_y", "box_aspect_z"):
        _require_finite(name.replace("_", "-"), getattr(args, name), positive=True)
    return args


def _crop_2d(
    x_m: np.ndarray,
    y_m: np.ndarray,
    values: np.ndarray,
    x_limit_m: float | None,
    y_limit_m: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_crop, values_crop = crop_horizontal_range(x_m, values, x_limit_m)
    y_crop, values_crop = crop_vertical_range(y_m, values_crop, y_limit_m)
    return x_crop, y_crop, values_crop


def _format_value(value: Any) -> str:
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.12g}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return str(value)


def _kv_text(items: OrderedDict[str, Any] | dict[str, Any]) -> str:
    return "".join(f"{key}={_format_value(value)}\n" for key, value in items.items())


class CurrentStripWorkflow:
    """Orchestrate one correction run while keeping physics helpers pure."""

    OUTPUT_SUFFIXES = (
        "_strip_currents.csv",
        "_strip_currents_one_plane.png",
        "_corrected_kickmap.dat",
        "_metrics.txt",
        "_plot_manifest.txt",
        "_theta_x_midplane_cut.png",
        "_theta_y_centerline_cut.png",
        "_raw_kickmaps_3d.png",
        "_corrected_kickmaps_3d.png",
    )

    def __init__(self, args: argparse.Namespace):
        self.args = validate_arguments(args)
        self.metrics: OrderedDict[str, Any] = OrderedDict()
        self.acceptance_failures: list[str] = []

    def run(self) -> None:
        self._prepare()
        self._solve()
        self._analyse()
        self._build_metrics()
        metrics_text = _kv_text(self.metrics)
        if self.acceptance_failures:
            print(metrics_text, end="", file=sys.stderr)
            raise RuntimeError(
                "Correction acceptance failed: " + "; ".join(self.acceptance_failures)
            )
        manifest_text = self._write_outputs(metrics_text)
        print(_kv_text(self.metrics), end="")
        print(manifest_text, end="")

    def _prepare(self) -> None:
        a = self.args
        self.beam_momentum_ev_c, self.beam_rigidity_tm = electron_beam_rigidity(
            a.beam_energy_ev
        )
        self.input_to_urad_scale = kick_input_to_urad_scale(
            a.units, self.beam_rigidity_tm, a.kick_scale
        )
        raw = load_kick_map(a.kick)
        self.kick_map = KickMap(
            raw.length_m,
            raw.x_m,
            raw.y_m,
            self.input_to_urad_scale * raw.theta_x_urad,
            self.input_to_urad_scale * raw.theta_y_urad,
            raw.b2_t2m,
        )
        self.maximum_map_x_m = float(np.max(np.abs(self.kick_map.x_m)))
        self.maximum_map_y_m = float(np.max(np.abs(self.kick_map.y_m)))
        self.model_y_half_span_m = (
            self.maximum_map_y_m if a.y_range_m is None else a.y_range_m
        )
        self.model_y_range_source = (
            "complete_input_y_grid" if a.y_range_m is None else "user"
        )

        if a.strip_layout in {"horizontal", "both"} and a.gap / 2.0 <= self.model_y_half_span_m:
            raise ValueError(
                "The top/bottom gap must enclose the model y aperture"
            )
        if a.side_gap is None:
            self.effective_side_gap_m = 2.0 * (a.x_range_m + a.width)
            self.side_gap_source = "auto_model_x_range_plus_one_width_margin"
        else:
            self.effective_side_gap_m = a.side_gap
            self.side_gap_source = "user"
        if a.strip_layout in {"vertical", "both"} and self.effective_side_gap_m / 2.0 <= a.x_range_m:
            raise ValueError(
                "The left/right gap must enclose the model x aperture"
            )

        strip_length_m = self.kick_map.length_m if a.length is None else a.length
        self.strips, self.strip_geometry = make_strip_arrays(
            layout=a.strip_layout,
            horizontal_half_span_m=a.effective_horizontal_strip_range_mm * 1.0e-3,
            horizontal_count_per_plane=a.strips_per_plane,
            vertical_half_span_m=a.effective_vertical_strip_range_mm * 1.0e-3,
            vertical_count_per_plane=a.effective_vertical_strips_per_plane,
            width_m=a.width,
            thickness_m=a.thickness,
            top_bottom_clear_gap_m=a.gap,
            left_right_clear_gap_m=self.effective_side_gap_m,
            length_m=strip_length_m,
            side_gap_source=self.side_gap_source,
        )
        self.horizontal_indices = np.asarray(
            [i for i, strip in enumerate(self.strips) if strip.family == "horizontal"],
            dtype=int,
        )
        self.vertical_indices = np.asarray(
            [i for i, strip in enumerate(self.strips) if strip.family == "vertical"],
            dtype=int,
        )

        self.model_x_m, self.model_y_m, self.raw_theta_x_model = _crop_2d(
            self.kick_map.x_m, self.kick_map.y_m, self.kick_map.theta_x_urad,
            a.x_range_m, a.y_range_m,
        )
        _, _, self.raw_theta_y_model = _crop_2d(
            self.kick_map.x_m, self.kick_map.y_m, self.kick_map.theta_y_urad,
            a.x_range_m, a.y_range_m,
        )
        if self.model_x_m.size < 2 or self.model_y_m.size < 2:
            raise ValueError("The model aperture must contain at least two x and y points")

        model_x_grid, model_y_grid = np.meshgrid(self.model_x_m, self.model_y_m)
        self.model_points = np.column_stack(
            [model_x_grid.ravel(), model_y_grid.ravel()]
        )
        self.response_iy_model = response_matrix(
            self.strips, self.model_points, "Iy", a.quadrature_order
        )
        self.response_ix_model = response_matrix(
            self.strips, self.model_points, "Ix", a.quadrature_order
        )
        self.raw_theta_x_flat = self.raw_theta_x_model.ravel()
        self.raw_theta_y_flat = self.raw_theta_y_model.ravel()
        self.target_iy_tm = self.beam_rigidity_tm * self.raw_theta_x_flat * 1.0e-6
        self.target_ix_tm = -self.beam_rigidity_tm * self.raw_theta_y_flat * 1.0e-6
        self.raw_theta_x_plane = affine_plane_decomposition(
            self.model_points, self.raw_theta_x_flat
        )
        self.raw_theta_y_plane = affine_plane_decomposition(
            self.model_points, self.raw_theta_y_flat
        )
        target_iy_plane = affine_plane_decomposition(self.model_points, self.target_iy_tm)
        target_ix_plane = affine_plane_decomposition(self.model_points, self.target_ix_tm)
        if a.target_mode == "nonlinear":
            self.fit_target_iy_tm = target_iy_plane.residual
            self.fit_target_ix_tm = target_ix_plane.residual
            self.horizontal_fit_basis = "full_2d_constant_and_linear_plane_removed"
            self.vertical_fit_basis = "full_2d_constant_and_linear_plane_removed"
        else:
            self.fit_target_iy_tm = self.target_iy_tm
            self.fit_target_ix_tm = self.target_ix_tm
            self.horizontal_fit_basis = "full_2d_complete_map"
            self.vertical_fit_basis = "full_2d_complete_map"
        self.fit_response_iy = self.response_iy_model
        self.fit_response_ix = self.response_ix_model

        weight_scale = max(a.horizontal_weight, a.vertical_weight)
        self.normalized_horizontal_weight = a.horizontal_weight / weight_scale
        self.normalized_vertical_weight = a.vertical_weight / weight_scale
        self.total_normalized_fit_weight = (
            self.normalized_horizontal_weight + self.normalized_vertical_weight
        )

    def _solve(self) -> None:
        a = self.args
        response_blocks: list[np.ndarray] = []
        target_blocks: list[np.ndarray] = []
        for weight, response, target in (
            (self.normalized_horizontal_weight, self.fit_response_iy, self.fit_target_iy_tm),
            (self.normalized_vertical_weight, self.fit_response_ix, self.fit_target_ix_tm),
        ):
            if weight > 0.0:
                scale = np.sqrt(weight / target.size)
                response_blocks.append(scale * response)
                target_blocks.append(scale * target)
        self.current_solution = solve_currents(
            np.vstack(response_blocks), np.concatenate(target_blocks), a.ridge, a.limit
        )
        self.currents_a = self.current_solution.currents_a

        x_grid, y_grid = np.meshgrid(self.kick_map.x_m, self.kick_map.y_m)
        all_points = np.column_stack([x_grid.ravel(), y_grid.ravel()])
        response_iy_full = response_matrix(
            self.strips, all_points, "Iy", a.quadrature_order
        )
        response_ix_full = response_matrix(
            self.strips, all_points, "Ix", a.quadrature_order
        )
        strip_iy_tm = (response_iy_full @ self.currents_a).reshape(x_grid.shape)
        strip_ix_tm = (response_ix_full @ self.currents_a).reshape(x_grid.shape)
        self.theta_x_corrected_urad = (
            self.kick_map.theta_x_urad - 1.0e6 * strip_iy_tm / self.beam_rigidity_tm
        )
        self.theta_y_corrected_urad = (
            self.kick_map.theta_y_urad + 1.0e6 * strip_ix_tm / self.beam_rigidity_tm
        )

    def _model_cut(self, values: np.ndarray, component: str) -> ZeroCut:
        a = self.args
        if component == "x":
            complete = zero_coordinate_cut(self.kick_map.y_m, values, 0, "y")
            _, cropped = crop_horizontal_range(
                self.kick_map.x_m, complete.values, a.x_range_m
            )
        else:
            complete = zero_coordinate_cut(self.kick_map.x_m, values, 1, "x")
            _, column = crop_vertical_range(
                self.kick_map.y_m, complete.values[:, np.newaxis], a.y_range_m
            )
            cropped = column[:, 0]
        return ZeroCut(
            cropped, complete.method,
            complete.lower_coordinate_m, complete.upper_coordinate_m,
        )

    def _analyse(self) -> None:
        a = self.args
        _, _, self.corrected_theta_x_model = _crop_2d(
            self.kick_map.x_m, self.kick_map.y_m, self.theta_x_corrected_urad,
            a.x_range_m, a.y_range_m,
        )
        _, _, self.corrected_theta_y_model = _crop_2d(
            self.kick_map.x_m, self.kick_map.y_m, self.theta_y_corrected_urad,
            a.x_range_m, a.y_range_m,
        )
        self.corrected_theta_x_flat = self.corrected_theta_x_model.ravel()
        self.corrected_theta_y_flat = self.corrected_theta_y_model.ravel()
        self.corrected_theta_x_plane = affine_plane_decomposition(
            self.model_points, self.corrected_theta_x_flat
        )
        self.corrected_theta_y_plane = affine_plane_decomposition(
            self.model_points, self.corrected_theta_y_flat
        )

        self.raw_theta_x_y0 = self._model_cut(self.kick_map.theta_x_urad, "x")
        self.raw_theta_y_x0 = self._model_cut(self.kick_map.theta_y_urad, "y")
        self.corrected_theta_x_y0 = self._model_cut(self.theta_x_corrected_urad, "x")
        self.corrected_theta_y_x0 = self._model_cut(self.theta_y_corrected_urad, "y")
        self.raw_theta_x_affine = affine_decomposition(
            self.model_x_m, self.raw_theta_x_y0.values
        )
        self.raw_theta_y_affine = affine_decomposition(
            self.model_y_m, self.raw_theta_y_x0.values
        )
        self.corrected_theta_x_affine = affine_decomposition(
            self.model_x_m, self.corrected_theta_x_y0.values
        )
        self.corrected_theta_y_affine = affine_decomposition(
            self.model_y_m, self.corrected_theta_y_x0.values
        )

        self.fit_prediction_iy_tm = self.fit_response_iy @ self.currents_a
        self.fit_prediction_ix_tm = self.fit_response_ix @ self.currents_a
        self.fit_residual_iy_tm = self.fit_prediction_iy_tm - self.fit_target_iy_tm
        self.fit_residual_ix_tm = self.fit_prediction_ix_tm - self.fit_target_ix_tm
        full_residual_iy = self.response_iy_model @ self.currents_a - self.target_iy_tm
        full_residual_ix = self.response_ix_model @ self.currents_a - self.target_ix_tm

        self.horizontal_fit_target_rms_tm = root_mean_square(self.fit_target_iy_tm)
        self.horizontal_fit_residual_rms_tm = root_mean_square(self.fit_residual_iy_tm)
        self.horizontal_fit_relative_rms = relative_rms_residual(
            self.fit_residual_iy_tm, self.fit_target_iy_tm
        )
        self.vertical_fit_target_rms_tm = root_mean_square(self.fit_target_ix_tm)
        self.vertical_fit_residual_rms_tm = root_mean_square(self.fit_residual_ix_tm)
        self.vertical_fit_relative_rms = relative_rms_residual(
            self.fit_residual_ix_tm, self.fit_target_ix_tm
        )
        self.horizontal_full_fit_relative_rms = relative_rms_residual(
            full_residual_iy, self.target_iy_tm
        )
        self.vertical_full_fit_relative_rms = relative_rms_residual(
            full_residual_ix, self.target_ix_tm
        )
        self.fit_target_rms_tm = self._weighted_rms(
            self.horizontal_fit_target_rms_tm, self.vertical_fit_target_rms_tm
        )
        self.fit_residual_rms_tm = self._weighted_rms(
            self.horizontal_fit_residual_rms_tm, self.vertical_fit_residual_rms_tm
        )
        self.fit_relative_rms = (
            0.0 if self.fit_target_rms_tm == self.fit_residual_rms_tm == 0.0
            else float("inf") if self.fit_target_rms_tm == 0.0
            else self.fit_residual_rms_tm / self.fit_target_rms_tm
        )

        self._compute_map_metrics()
        self._compute_centerline_metrics()
        self._evaluate_acceptance()

    def _weighted_rms(self, horizontal: float, vertical: float) -> float:
        return float(np.sqrt(
            (
                self.normalized_horizontal_weight * horizontal**2
                + self.normalized_vertical_weight * vertical**2
            ) / self.total_normalized_fit_weight
        ))

    def _reduction(self, raw: float, corrected: float, prefix: str) -> None:
        factor, status = safe_reduction_factor(raw, corrected)
        setattr(self, f"{prefix}_factor", factor)
        setattr(self, f"{prefix}_status", status)

    def _compute_map_metrics(self) -> None:
        self.raw_theta_x_model_rms = root_mean_square(self.raw_theta_x_flat)
        self.corrected_theta_x_model_rms = root_mean_square(self.corrected_theta_x_flat)
        self.raw_theta_y_model_rms = root_mean_square(self.raw_theta_y_flat)
        self.corrected_theta_y_model_rms = root_mean_square(self.corrected_theta_y_flat)
        self._reduction(
            self.raw_theta_x_model_rms, self.corrected_theta_x_model_rms,
            "full_map_theta_x_rms_reduction",
        )
        self._reduction(
            self.raw_theta_y_model_rms, self.corrected_theta_y_model_rms,
            "full_map_theta_y_rms_reduction",
        )
        self.raw_theta_x_model_peak = float(np.max(np.abs(self.raw_theta_x_flat)))
        self.corrected_theta_x_model_peak = float(
            np.max(np.abs(self.corrected_theta_x_flat))
        )
        self.raw_theta_y_model_peak = float(np.max(np.abs(self.raw_theta_y_flat)))
        self.corrected_theta_y_model_peak = float(
            np.max(np.abs(self.corrected_theta_y_flat))
        )
        self._reduction(
            self.raw_theta_x_model_peak, self.corrected_theta_x_model_peak,
            "full_map_theta_x_peak_reduction",
        )
        self._reduction(
            self.raw_theta_y_model_peak, self.corrected_theta_y_model_peak,
            "full_map_theta_y_peak_reduction",
        )
        self.raw_theta_x_nonlinear_rms = root_mean_square(
            self.raw_theta_x_plane.residual
        )
        self.corrected_theta_x_nonlinear_rms = root_mean_square(
            self.corrected_theta_x_plane.residual
        )
        self.raw_theta_y_nonlinear_rms = root_mean_square(
            self.raw_theta_y_plane.residual
        )
        self.corrected_theta_y_nonlinear_rms = root_mean_square(
            self.corrected_theta_y_plane.residual
        )
        self._reduction(
            self.raw_theta_x_nonlinear_rms, self.corrected_theta_x_nonlinear_rms,
            "horizontal_nonlinear_rms_reduction",
        )
        self._reduction(
            self.raw_theta_y_nonlinear_rms, self.corrected_theta_y_nonlinear_rms,
            "vertical_nonlinear_rms_reduction",
        )

    def _compute_centerline_metrics(self) -> None:
        self.raw_theta_x_rms = root_mean_square(self.raw_theta_x_y0.values)
        self.corrected_theta_x_rms = root_mean_square(
            self.corrected_theta_x_y0.values
        )
        self._reduction(
            self.raw_theta_x_rms, self.corrected_theta_x_rms,
            "horizontal_rms_reduction",
        )
        self.raw_theta_x_y0_nonlinear_rms = root_mean_square(
            self.raw_theta_x_affine.residual
        )
        self.corrected_theta_x_y0_nonlinear_rms = root_mean_square(
            self.corrected_theta_x_affine.residual
        )
        self._reduction(
            self.raw_theta_x_y0_nonlinear_rms,
            self.corrected_theta_x_y0_nonlinear_rms,
            "horizontal_y0_nonlinear_rms_reduction",
        )
        self.raw_theta_y_peak = float(np.max(np.abs(self.raw_theta_y_x0.values)))
        self.corrected_theta_y_peak = float(
            np.max(np.abs(self.corrected_theta_y_x0.values))
        )
        self._reduction(
            self.raw_theta_y_peak, self.corrected_theta_y_peak,
            "vertical_peak_reduction",
        )
        self.raw_theta_y_nonlinear_peak = float(
            np.max(np.abs(self.raw_theta_y_affine.residual))
        )
        self.corrected_theta_y_nonlinear_peak = float(
            np.max(np.abs(self.corrected_theta_y_affine.residual))
        )
        self._reduction(
            self.raw_theta_y_nonlinear_peak,
            self.corrected_theta_y_nonlinear_peak,
            "vertical_nonlinear_peak_reduction",
        )

    def _evaluate_acceptance(self) -> None:
        a = self.args
        if a.target_mode == "nonlinear":
            self.horizontal_acceptance_reduction_factor = (
                self.horizontal_nonlinear_rms_reduction_factor
            )
            self.vertical_acceptance_reduction_factor = (
                self.vertical_nonlinear_rms_reduction_factor
            )
            self.horizontal_acceptance_reduction_status = (
                self.horizontal_nonlinear_rms_reduction_status
            )
            self.vertical_acceptance_reduction_status = (
                self.vertical_nonlinear_rms_reduction_status
            )
            self.horizontal_acceptance_rms_basis = "full_2d_affine_plane_removed"
            self.vertical_acceptance_rms_basis = "full_2d_affine_plane_removed"
        else:
            self.horizontal_acceptance_reduction_factor = (
                self.full_map_theta_x_rms_reduction_factor
            )
            self.vertical_acceptance_reduction_factor = (
                self.full_map_theta_y_rms_reduction_factor
            )
            self.horizontal_acceptance_reduction_status = (
                self.full_map_theta_x_rms_reduction_status
            )
            self.vertical_acceptance_reduction_status = (
                self.full_map_theta_y_rms_reduction_status
            )
            self.horizontal_acceptance_rms_basis = "full_2d_complete_map"
            self.vertical_acceptance_rms_basis = "full_2d_complete_map"

        bound_threshold = a.limit * (1.0 - a.current_bound_tolerance)
        active = np.abs(self.currents_a) >= bound_threshold
        self.current_bound_active_count = int(np.count_nonzero(active))
        self.current_bound_active_fraction = float(
            self.current_bound_active_count / self.currents_a.size
        )

        if self.fit_relative_rms > a.effective_max_fit_relative_rms:
            self.acceptance_failures.append(
                f"fit relative RMS {self.fit_relative_rms:.6g} exceeds "
                f"{a.effective_max_fit_relative_rms:.6g}"
            )
        component_checks = (
            ("horizontal", a.horizontal_weight,
             self.horizontal_acceptance_reduction_factor,
             a.min_horizontal_rms_reduction_factor),
            ("vertical", a.vertical_weight,
             self.vertical_acceptance_reduction_factor,
             a.min_vertical_rms_reduction_factor),
        )
        for name, weight, value, minimum in component_checks:
            if weight > 0.0 and value < minimum:
                self.acceptance_failures.append(
                    f"{name} active-basis 2D RMS reduction factor "
                    f"{value:.6g} is below {minimum:.6g}"
                )
        for name, value in (
            ("horizontal", self.full_map_theta_x_rms_reduction_factor),
            ("vertical", self.full_map_theta_y_rms_reduction_factor),
        ):
            if value < a.min_full_map_rms_reduction_factor:
                self.acceptance_failures.append(
                    f"{name} full-map RMS reduction factor {value:.6g} is below "
                    f"{a.min_full_map_rms_reduction_factor:.6g}"
                )
        if a.min_full_map_peak_reduction_factor > 0.0:
            for name, value in (
                ("horizontal", self.full_map_theta_x_peak_reduction_factor),
                ("vertical", self.full_map_theta_y_peak_reduction_factor),
            ):
                if value < a.min_full_map_peak_reduction_factor:
                    self.acceptance_failures.append(
                        f"{name} full-map peak reduction factor {value:.6g} is below "
                        f"{a.min_full_map_peak_reduction_factor:.6g}"
                    )
        self.acceptance_status = (
            "FAILED" if self.acceptance_failures
            else "ACCEPTED_LIMIT_SATURATED" if self.current_bound_active_count
            else "ACCEPTED"
        )

    def _build_metrics(self) -> None:
        a, g, s = self.args, self.strip_geometry, self.current_solution
        m = self.metrics
        add = m.__setitem__
        base_items = [
            ("release_id", RELEASE_ID), ("input_file", a.kick),
            ("grid_shape", self.kick_map.theta_x_urad.shape),
            ("model_domain_source", "x_range_m_and_y_range_m"),
            ("model_x_requested_half_span_m", a.x_range_m),
            ("model_x_point_count", self.model_x_m.size),
            ("model_x_min_m", np.min(self.model_x_m)),
            ("model_x_max_m", np.max(self.model_x_m)),
            ("model_y_requested_half_span_m", self.model_y_half_span_m),
            ("model_y_range_source", self.model_y_range_source),
            ("model_y_point_count", self.model_y_m.size),
            ("model_y_min_m", np.min(self.model_y_m)),
            ("model_y_max_m", np.max(self.model_y_m)),
            ("export_grid_scope", "complete_input_grid"),
            ("length_m", self.kick_map.length_m),
            ("input_kick_units", a.units), ("beam_energy_eV", a.beam_energy_ev),
            ("electron_rest_energy_eV", ELECTRON_REST_ENERGY_EV),
            ("beam_momentum_eV_c", self.beam_momentum_ev_c),
            ("beam_rigidity_Tm", self.beam_rigidity_tm),
            ("input_to_urad_scale", self.input_to_urad_scale),
            ("kick_scale", a.kick_scale),
            ("horizontal_weight", a.horizontal_weight),
            ("vertical_weight", a.vertical_weight), ("target_mode", a.target_mode),
            ("horizontal_fit_basis", self.horizontal_fit_basis),
            ("vertical_fit_basis", self.vertical_fit_basis),
            ("fit_domain", "full_2d_model_aperture"),
            ("fit_point_count", self.model_points.shape[0]),
            ("affine_projection_basis_2d", "constant_plus_x_plus_y"),
            ("normalized_horizontal_weight", self.normalized_horizontal_weight),
            ("normalized_vertical_weight", self.normalized_vertical_weight),
            ("fit_weighting", "max_normalized_component_weight_times_mean_squared_residual"),
            ("color_scale", a.color_scale), ("plot_elevation_deg", a.elevation_deg),
            ("plot_azimuth_deg", a.azimuth_deg), ("projection", a.projection),
            ("box_aspect_x", a.box_aspect_x), ("box_aspect_y", a.box_aspect_y),
            ("box_aspect_z", a.box_aspect_z),
        ]
        for key, value in base_items:
            add(key, value)

        for prefix, cut in (("theta_x_y0", self.raw_theta_x_y0),
                            ("theta_y_x0", self.raw_theta_y_x0)):
            add(f"{prefix}_cut_method", cut.method)
            add(f"{prefix}_bracket_lower_m", cut.lower_coordinate_m)
            add(f"{prefix}_bracket_upper_m", cut.upper_coordinate_m)
        for key, value in (
            ("raw_theta_x_y0_min_urad", np.min(self.raw_theta_x_y0.values)),
            ("raw_theta_x_y0_max_urad", np.max(self.raw_theta_x_y0.values)),
            ("raw_theta_y_x0_min_urad", np.min(self.raw_theta_y_x0.values)),
            ("raw_theta_y_x0_max_urad", np.max(self.raw_theta_y_x0.values)),
        ):
            add(key, value)

        affine_items = [
            ("raw_theta_x_affine", self.raw_theta_x_affine),
            ("corrected_theta_x_affine", self.corrected_theta_x_affine),
            ("raw_theta_y_affine", self.raw_theta_y_affine),
            ("corrected_theta_y_affine", self.corrected_theta_y_affine),
        ]
        for prefix, decomposition in affine_items:
            add(f"{prefix}_offset_urad", float(decomposition.offset))
            add(
                f"{prefix}_gradient_urad_per_mm",
                float(decomposition.slope_per_m) * 1.0e-3,
            )

        diagnostic_items = [
            ("raw_theta_x_y0_rms_urad", self.raw_theta_x_rms),
            ("corrected_theta_x_y0_rms_urad", self.corrected_theta_x_rms),
            ("horizontal_rms_reduction_factor", self.horizontal_rms_reduction_factor),
            ("horizontal_rms_reduction_status", self.horizontal_rms_reduction_status),
            ("raw_theta_x_y0_nonlinear_rms_urad", self.raw_theta_x_y0_nonlinear_rms),
            ("corrected_theta_x_y0_nonlinear_rms_urad", self.corrected_theta_x_y0_nonlinear_rms),
            ("horizontal_y0_nonlinear_rms_reduction_factor", self.horizontal_y0_nonlinear_rms_reduction_factor),
            ("horizontal_y0_nonlinear_rms_reduction_status", self.horizontal_y0_nonlinear_rms_reduction_status),
            ("horizontal_acceptance_rms_basis", self.horizontal_acceptance_rms_basis),
            ("horizontal_acceptance_reduction_factor", self.horizontal_acceptance_reduction_factor),
            ("horizontal_acceptance_reduction_status", self.horizontal_acceptance_reduction_status),
            ("vertical_acceptance_rms_basis", self.vertical_acceptance_rms_basis),
            ("vertical_acceptance_reduction_factor", self.vertical_acceptance_reduction_factor),
            ("vertical_acceptance_reduction_status", self.vertical_acceptance_reduction_status),
            ("raw_theta_y_x0_peak_urad", self.raw_theta_y_peak),
            ("corrected_theta_y_x0_peak_urad", self.corrected_theta_y_peak),
            ("vertical_peak_reduction_factor", self.vertical_peak_reduction_factor),
            ("vertical_peak_reduction_status", self.vertical_peak_reduction_status),
            ("raw_theta_y_x0_nonlinear_peak_urad", self.raw_theta_y_nonlinear_peak),
            ("corrected_theta_y_x0_nonlinear_peak_urad", self.corrected_theta_y_nonlinear_peak),
            ("vertical_nonlinear_peak_reduction_factor", self.vertical_nonlinear_peak_reduction_factor),
            ("vertical_nonlinear_peak_reduction_status", self.vertical_nonlinear_peak_reduction_status),
        ]
        for key, value in diagnostic_items:
            add(key, value)

        map_items = [
            ("raw_theta_x_model_2d_rms_urad", self.raw_theta_x_model_rms),
            ("corrected_theta_x_model_2d_rms_urad", self.corrected_theta_x_model_rms),
            ("full_map_theta_x_rms_reduction_factor", self.full_map_theta_x_rms_reduction_factor),
            ("full_map_theta_x_rms_reduction_status", self.full_map_theta_x_rms_reduction_status),
            ("raw_theta_y_model_2d_rms_urad", self.raw_theta_y_model_rms),
            ("corrected_theta_y_model_2d_rms_urad", self.corrected_theta_y_model_rms),
            ("full_map_theta_y_rms_reduction_factor", self.full_map_theta_y_rms_reduction_factor),
            ("full_map_theta_y_rms_reduction_status", self.full_map_theta_y_rms_reduction_status),
            ("raw_theta_x_model_2d_peak_urad", self.raw_theta_x_model_peak),
            ("corrected_theta_x_model_2d_peak_urad", self.corrected_theta_x_model_peak),
            ("full_map_theta_x_peak_reduction_factor", self.full_map_theta_x_peak_reduction_factor),
            ("full_map_theta_x_peak_reduction_status", self.full_map_theta_x_peak_reduction_status),
            ("raw_theta_y_model_2d_peak_urad", self.raw_theta_y_model_peak),
            ("corrected_theta_y_model_2d_peak_urad", self.corrected_theta_y_model_peak),
            ("full_map_theta_y_peak_reduction_factor", self.full_map_theta_y_peak_reduction_factor),
            ("full_map_theta_y_peak_reduction_status", self.full_map_theta_y_peak_reduction_status),
            ("raw_theta_x_model_2d_nonlinear_rms_urad", self.raw_theta_x_nonlinear_rms),
            ("corrected_theta_x_model_2d_nonlinear_rms_urad", self.corrected_theta_x_nonlinear_rms),
            ("raw_theta_y_model_2d_nonlinear_rms_urad", self.raw_theta_y_nonlinear_rms),
            ("corrected_theta_y_model_2d_nonlinear_rms_urad", self.corrected_theta_y_nonlinear_rms),
        ]
        for key, value in map_items:
            add(key, value)

        for component, raw_plane, corrected_plane in (
            ("theta_x", self.raw_theta_x_plane, self.corrected_theta_x_plane),
            ("theta_y", self.raw_theta_y_plane, self.corrected_theta_y_plane),
        ):
            for label, plane in (("raw", raw_plane), ("corrected", corrected_plane)):
                prefix = f"{label}_{component}_affine_plane"
                add(f"{prefix}_offset_urad", float(plane.offset))
                add(f"{prefix}_gradient_x_urad_per_mm", float(plane.slope_x_per_m) * 1.0e-3)
                add(f"{prefix}_gradient_y_urad_per_mm", float(plane.slope_y_per_m) * 1.0e-3)

        solve_items = [
            ("solver_status", s.status), ("solver_message", s.message),
            ("solver_cost", s.cost), ("solver_optimality", s.optimality),
            ("solver_iterations", s.iterations),
            ("horizontal_fit_target_rms_Tm", self.horizontal_fit_target_rms_tm),
            ("horizontal_fit_residual_rms_Tm", self.horizontal_fit_residual_rms_tm),
            ("horizontal_fit_relative_rms", self.horizontal_fit_relative_rms),
            ("vertical_fit_target_rms_Tm", self.vertical_fit_target_rms_tm),
            ("vertical_fit_residual_rms_Tm", self.vertical_fit_residual_rms_tm),
            ("vertical_fit_relative_rms", self.vertical_fit_relative_rms),
            ("horizontal_full_fit_relative_rms", self.horizontal_full_fit_relative_rms),
            ("vertical_full_fit_relative_rms", self.vertical_full_fit_relative_rms),
            ("fit_target_rms_Tm", self.fit_target_rms_tm),
            ("fit_residual_rms_Tm", self.fit_residual_rms_tm),
            ("fit_relative_rms", self.fit_relative_rms),
            ("max_fit_relative_rms", a.effective_max_fit_relative_rms),
            ("max_fit_relative_rms_source", a.max_fit_relative_rms_source),
            ("min_horizontal_rms_reduction_factor", a.min_horizontal_rms_reduction_factor),
            ("min_vertical_rms_reduction_factor", a.min_vertical_rms_reduction_factor),
            ("min_full_map_rms_reduction_factor", a.min_full_map_rms_reduction_factor),
            ("min_full_map_peak_reduction_factor", a.min_full_map_peak_reduction_factor),
            ("current_bound_tolerance", a.current_bound_tolerance),
            ("current_bound_active_count", self.current_bound_active_count),
            ("current_bound_active_fraction", self.current_bound_active_fraction),
            ("acceptance_status", self.acceptance_status),
        ]
        for key, value in solve_items:
            add(key, value)

        geometry_items = [
            ("strip_layout", g.layout),
            ("horizontal_array_enabled", int(self.horizontal_indices.size > 0)),
            ("vertical_array_enabled", int(self.vertical_indices.size > 0)),
            ("horizontal_strips_per_plane", g.horizontal_count_per_plane),
            ("horizontal_strip_range_option_mm", a.effective_horizontal_strip_range_mm),
            ("vertical_strips_per_plane", g.vertical_count_per_plane),
            ("horizontal_strip_half_span_mm", g.horizontal_half_span_m * 1.0e3),
            ("vertical_strip_half_span_mm", g.vertical_half_span_m * 1.0e3),
            ("horizontal_strip_pitch_mm", g.horizontal_pitch_m * 1.0e3),
            ("vertical_strip_pitch_mm", g.vertical_pitch_m * 1.0e3),
            ("top_bottom_clear_gap_mm", g.top_bottom_clear_gap_m * 1.0e3),
            ("left_right_clear_gap_mm", g.left_right_clear_gap_m * 1.0e3),
            ("side_gap_source", g.side_gap_source),
            ("strip_width_mm", a.width * 1.0e3),
            ("strip_thickness_mm", a.thickness * 1.0e3),
            ("total_strip_count", len(self.strips)),
            ("max_abs_current_A", np.max(np.abs(self.currents_a))),
            ("max_abs_horizontal_array_current_A", self._family_max(self.horizontal_indices)),
            ("max_abs_vertical_array_current_A", self._family_max(self.vertical_indices)),
        ]
        for key, value in geometry_items:
            add(key, value)

    def _family_max(self, indices: np.ndarray) -> float:
        return float(np.max(np.abs(self.currents_a[indices]))) if indices.size else float("nan")

    def _write_outputs(self, metrics_text: str) -> str:
        a = self.args
        prefix = a.out_prefix
        prefix.parent.mkdir(parents=True, exist_ok=True)
        final_paths = [prefix.with_name(prefix.name + suffix) for suffix in self.OUTPUT_SUFFIXES]
        with tempfile.TemporaryDirectory(
            dir=prefix.parent, prefix=f".{prefix.name}_staging_"
        ) as staging_directory:
            staging_prefix = Path(staging_directory) / prefix.name
            staged_paths = [
                staging_prefix.with_name(staging_prefix.name + suffix)
                for suffix in self.OUTPUT_SUFFIXES
            ]
            (
                currents_path, currents_plot_path, corrected_map_path,
                metrics_path, manifest_path, theta_x_cut_path,
                theta_y_cut_path, raw_3d_path, corrected_3d_path,
            ) = staged_paths
            write_strip_currents_csv(currents_path, self.strips, self.currents_a)
            currents_figure, symmetry = plot_strip_currents_representative_planes(
                self.strips, self.currents_a, currents_plot_path,
                a.effective_horizontal_strip_range_mm,
                a.effective_vertical_strip_range_mm,
                keep_open=a.show,
            )
            write_corrected_kick_map(
                corrected_map_path, self.kick_map,
                self.theta_x_corrected_urad, self.theta_y_corrected_urad,
                source_units=a.units, beam_energy_ev=a.beam_energy_ev,
                beam_rigidity_tm=self.beam_rigidity_tm, kick_scale=a.kick_scale,
            )

            if a.target_mode == "nonlinear":
                x_raw, x_corrected = (
                    self.raw_theta_x_affine.residual,
                    self.corrected_theta_x_affine.residual,
                )
                y_raw, y_corrected = (
                    self.raw_theta_y_affine.residual,
                    self.corrected_theta_y_affine.residual,
                )
                x_title = (
                    "Horizontal nonlinear centreline diagnostic on y = 0 "
                    "(constant and gradient removed)"
                )
                y_title = (
                    "Vertical nonlinear centreline diagnostic on x = 0 "
                    "(constant and gradient removed)"
                )
            else:
                x_raw, x_corrected = self.raw_theta_x_y0.values, self.corrected_theta_x_y0.values
                y_raw, y_corrected = self.raw_theta_y_x0.values, self.corrected_theta_y_x0.values
                x_title = "Horizontal kick on the horizontal mid-plane (y = 0)"
                y_title = "Vertical kick on the vertical centreline (x = 0)"
            plot_cut(
                self.model_x_m * 1.0e3, x_raw, x_corrected, x_title,
                "x [mm]", "θx [µrad]", theta_x_cut_path,
            )
            plot_cut(
                self.model_y_m * 1.0e3, y_raw, y_corrected, y_title,
                "y [mm]", "θy [µrad]", theta_y_cut_path,
            )

            theta_x_limit = max(
                float(np.max(np.abs(self.raw_theta_x_model))),
                float(np.max(np.abs(self.corrected_theta_x_model))),
            )
            theta_y_limit = max(
                float(np.max(np.abs(self.raw_theta_y_model))),
                float(np.max(np.abs(self.corrected_theta_y_model))),
            )
            plot_args = dict(
                x_range_m=a.x_range_m, y_range_m=a.y_range_m,
                color_scale=a.color_scale, elevation_deg=a.elevation_deg,
                azimuth_deg=a.azimuth_deg, projection=a.projection,
                box_aspect_x=a.box_aspect_x, box_aspect_y=a.box_aspect_y,
                box_aspect_z=a.box_aspect_z, keep_open=a.show,
            )
            raw_figure = plot_kickmaps_3d_pair(
                self.kick_map.x_m, self.kick_map.y_m,
                self.kick_map.theta_x_urad, self.kick_map.theta_y_urad,
                raw_3d_path, theta_x_limit, theta_y_limit,
                "Kick maps before current-strip correction", **plot_args,
            )
            corrected_figure = plot_kickmaps_3d_pair(
                self.kick_map.x_m, self.kick_map.y_m,
                self.theta_x_corrected_urad, self.theta_y_corrected_urad,
                corrected_3d_path, theta_x_limit, theta_y_limit,
                "Kick maps after current-strip correction", **plot_args,
            )

            for family in ("horizontal", "vertical"):
                max_mismatch, rms_mismatch = symmetry.get(
                    family, (float("nan"), float("nan"))
                )
                self.metrics[f"{family}_sheet_symmetry_max_abs_mismatch_A"] = max_mismatch
                self.metrics[f"{family}_sheet_symmetry_rms_mismatch_A"] = rms_mismatch
            metrics_text = _kv_text(self.metrics)
            metrics_path.write_text(metrics_text, encoding="utf-8", newline="\n")

            manifest = self._build_manifest()
            manifest_text = _kv_text(manifest)
            manifest_path.write_text(manifest_text, encoding="utf-8", newline="\n")
            if a.show:
                plt.show()
                for figure in (currents_figure, raw_figure, corrected_figure):
                    if figure is not None:
                        plt.close(figure)
            publish_staged_outputs(staged_paths, final_paths)
        return manifest_text

    def _build_manifest(self) -> OrderedDict[str, Any]:
        a, g = self.args, self.strip_geometry
        return OrderedDict([
            ("release_id", RELEASE_ID), ("acceptance_status", self.acceptance_status),
            ("model_x_half_span_m", a.x_range_m),
            ("model_y_half_span_m", self.model_y_half_span_m),
            ("model_x_point_count", self.model_x_m.size),
            ("model_y_point_count", self.model_y_m.size),
            ("export_grid_scope", "complete_input_grid"),
            ("input_kick_units", a.units), ("beam_energy_eV", a.beam_energy_ev),
            ("beam_rigidity_Tm", self.beam_rigidity_tm),
            ("input_to_urad_scale", self.input_to_urad_scale),
            ("kick_scale", a.kick_scale),
            ("horizontal_weight", a.horizontal_weight),
            ("vertical_weight", a.vertical_weight),
            ("target_mode", a.target_mode),
            ("fit_domain", "full_2d_model_aperture"),
            ("fit_point_count", self.model_points.shape[0]),
            ("horizontal_fit_basis", self.horizontal_fit_basis),
            ("vertical_fit_basis", self.vertical_fit_basis),
            ("normalized_horizontal_weight", self.normalized_horizontal_weight),
            ("normalized_vertical_weight", self.normalized_vertical_weight),
            ("strip_layout", g.layout),
            ("strip_currents_representative_planes_file", f"{a.out_prefix.name}_strip_currents_one_plane.png"),
            ("horizontal_strip_half_span_mm", g.horizontal_half_span_m * 1.0e3),
            ("vertical_strip_half_span_mm", g.vertical_half_span_m * 1.0e3),
            ("horizontal_strip_pitch_mm", g.horizontal_pitch_m * 1.0e3),
            ("vertical_strip_pitch_mm", g.vertical_pitch_m * 1.0e3),
            ("horizontal_strips_per_plane", g.horizontal_count_per_plane),
            ("vertical_strips_per_plane", g.vertical_count_per_plane),
            ("top_bottom_clear_gap_mm", g.top_bottom_clear_gap_m * 1.0e3),
            ("left_right_clear_gap_mm", g.left_right_clear_gap_m * 1.0e3),
            ("raw_kickmaps_file", f"{a.out_prefix.name}_raw_kickmaps_3d.png"),
            ("corrected_kickmaps_file", f"{a.out_prefix.name}_corrected_kickmaps_3d.png"),
            ("plot_elevation_deg", a.elevation_deg),
            ("plot_azimuth_deg", a.azimuth_deg),
        ])


def main() -> None:
    CurrentStripWorkflow(build_parser().parse_args()).run()


if __name__ == "__main__":
    main()
