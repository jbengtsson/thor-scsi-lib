#!/usr/bin/env python3
"""Current-strip correction solver with a full 2D fitting objective.

This refactor preserves the V12 numerical model while separating configuration,
geometry, solving, diagnostics, acceptance, and output publication into small
testable units. It supports both standard ``--help`` and a dedicated
``help [OPTION]`` command. All command-line range and geometry lengths are
expressed in metres. ``--x-range`` and ``--y-range`` each control both
the corresponding fitted model half-aperture and strip-centre half-span.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import tempfile
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
RELEASE_ID = "CURRENT-STRIP-SOLVER-28-CONSOLIDATED-RANGES"
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


@dataclass
class ComponentState:
    """Static component metadata plus computed per-component run state."""

    key: str
    label: str
    raw_map: np.ndarray
    response_component: str
    target_sign: float
    correction_sign: float
    weight: float

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

def _append_strip_family(
    strips: list[Strip],
    centres_m: np.ndarray,
    *,
    family: str,
    sheets: tuple[tuple[str, float], tuple[str, float]],
    fixed_centre_m: float,
    width_m: float,
    thickness_m: float,
    length_m: float,
) -> None:
    horizontal = family == "horizontal"
    for sheet, sign in sheets:
        for coordinate in centres_m:
            coordinate_m = float(coordinate)
            strips.append(Strip(
                x_m=coordinate_m if horizontal else float(sign * fixed_centre_m),
                y_m=float(sign * fixed_centre_m) if horizontal else coordinate_m,
                size_x_m=width_m if horizontal else thickness_m,
                size_y_m=thickness_m if horizontal else width_m,
                length_m=length_m,
                family=family,
                sheet=sheet,
                coordinate_m=coordinate_m,
            ))


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
    """Build horizontal top/bottom and/or vertical left/right strip sheets."""
    if layout not in {"horizontal", "vertical", "both"}:
        raise ValueError(f"unknown strip layout: {layout!r}")
    dimensions = (width_m, thickness_m, length_m)
    if not all(np.isfinite(value) and value > 0.0 for value in dimensions):
        raise ValueError("strip width, thickness, and length must be finite and positive")

    strips: list[Strip] = []
    horizontal_count = vertical_count = 0
    horizontal_pitch_m = vertical_pitch_m = float("nan")

    if layout in {"horizontal", "both"}:
        if not np.isfinite(top_bottom_clear_gap_m) or top_bottom_clear_gap_m <= 0.0:
            raise ValueError("top/bottom clear gap must be finite and positive")
        centres, horizontal_pitch_m = _uniform_strip_centres(
            horizontal_half_span_m, horizontal_count_per_plane, width_m,
            "horizontal-array",
        )
        horizontal_count = horizontal_count_per_plane
        _append_strip_family(
            strips,
            centres,
            family="horizontal",
            sheets=(("top", +1.0), ("bottom", -1.0)),
            fixed_centre_m=top_bottom_clear_gap_m / 2.0 + thickness_m / 2.0,
            width_m=width_m,
            thickness_m=thickness_m,
            length_m=length_m,
        )

    if layout in {"vertical", "both"}:
        if not np.isfinite(left_right_clear_gap_m) or left_right_clear_gap_m <= 0.0:
            raise ValueError("left/right clear gap must be finite and positive")
        centres, vertical_pitch_m = _uniform_strip_centres(
            vertical_half_span_m, vertical_count_per_plane, width_m,
            "vertical-array",
        )
        vertical_count = vertical_count_per_plane
        _append_strip_family(
            strips,
            centres,
            family="vertical",
            sheets=(("right", +1.0), ("left", -1.0)),
            fixed_centre_m=left_right_clear_gap_m / 2.0 + thickness_m / 2.0,
            width_m=width_m,
            thickness_m=thickness_m,
            length_m=length_m,
        )

    if not strips:
        raise ValueError("strip layout produced no conductors")
    nan = float("nan")
    return strips, StripGeometry(
        layout=layout,
        horizontal_count_per_plane=horizontal_count,
        vertical_count_per_plane=vertical_count,
        horizontal_half_span_m=horizontal_half_span_m if horizontal_count else nan,
        vertical_half_span_m=vertical_half_span_m if vertical_count else nan,
        horizontal_pitch_m=horizontal_pitch_m,
        vertical_pitch_m=vertical_pitch_m,
        top_bottom_clear_gap_m=top_bottom_clear_gap_m if horizontal_count else nan,
        left_right_clear_gap_m=left_right_clear_gap_m if vertical_count else nan,
        side_gap_source=side_gap_source if vertical_count else "not_applicable",
    )

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

def _project_least_squares(
    design: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Project the leading values axis onto a least-squares design matrix."""
    array = np.asarray(values, dtype=float)
    if array.shape[0] != design.shape[0]:
        raise ValueError("the first values axis must match the projection design")
    if not np.all(np.isfinite(design)) or not np.all(np.isfinite(array)):
        raise ValueError("projection inputs must be finite")

    original_shape = array.shape
    flattened = array.reshape(design.shape[0], -1)
    coefficients, _, rank, _ = np.linalg.lstsq(design, flattened, rcond=None)
    fitted_flat = design @ coefficients
    residual_flat = flattened - fitted_flat
    column_scale = np.max(np.abs(flattened), axis=0)
    tolerance = (
        256.0
        * np.finfo(float).eps
        * np.maximum(column_scale, np.finfo(float).tiny)
    )
    residual_flat[np.abs(residual_flat) <= tolerance] = 0.0
    return (
        fitted_flat.reshape(original_shape),
        residual_flat.reshape(original_shape),
        coefficients,
        int(rank),
    )


def affine_decomposition(
    coordinate_m: np.ndarray,
    values: np.ndarray,
) -> AffineDecomposition:
    """Split values into a least-squares line and nonlinear residual."""
    coordinate = np.asarray(coordinate_m, dtype=float)
    if coordinate.ndim != 1 or coordinate.size < 2:
        raise ValueError("affine projection requires at least two coordinates")
    if not np.all(np.isfinite(coordinate)) or np.ptp(coordinate) <= 0.0:
        raise ValueError("affine-fit coordinates must be finite and span an interval")

    scale_m = float(np.max(np.abs(coordinate))) or float(np.ptp(coordinate))
    design = np.column_stack(
        [np.ones(coordinate.size, dtype=float), coordinate / scale_m]
    )
    fitted, residual, coefficients, _ = _project_least_squares(design, values)
    trailing_shape = np.asarray(values).shape[1:]
    return AffineDecomposition(
        fitted=fitted,
        residual=residual,
        offset=coefficients[0].reshape(trailing_shape),
        slope_per_m=(coefficients[1] / scale_m).reshape(trailing_shape),
    )


def affine_plane_decomposition(
    points_xy_m: np.ndarray,
    values: np.ndarray,
) -> AffinePlaneDecomposition:
    """Split values into ``a + b*x + c*y`` and a nonlinear residual."""
    points = np.asarray(points_xy_m, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or points.shape[0] < 3:
        raise ValueError("affine-plane projection requires at least three 2D points")
    x, y = points.T
    x_scale_m = max(float(np.max(np.abs(x))), float(np.ptp(x)))
    y_scale_m = max(float(np.max(np.abs(y))), float(np.ptp(y)))
    if x_scale_m <= 0.0 or y_scale_m <= 0.0:
        raise ValueError("affine-plane points must span both x and y")

    design = np.column_stack(
        [np.ones(points.shape[0], dtype=float), x / x_scale_m, y / y_scale_m]
    )
    fitted, residual, coefficients, rank = _project_least_squares(design, values)
    if rank < 3:
        raise ValueError("affine-plane design is rank deficient")
    trailing_shape = np.asarray(values).shape[1:]
    return AffinePlaneDecomposition(
        fitted=fitted,
        residual=residual,
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
        staged_path.replace(final_path)

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

def _crop_coordinate_range(
    coordinates_m: np.ndarray,
    values: np.ndarray,
    limit_m: float | None,
    *,
    axis: int,
    coordinate_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Crop one values axis to a symmetric coordinate interval in metres."""
    if limit_m is None:
        return coordinates_m, values
    if limit_m <= 0.0:
        raise ValueError(f"{coordinate_name}_range must be positive")
    mask = np.abs(coordinates_m) <= limit_m + 1.0e-15
    if not np.any(mask):
        raise ValueError(
            f"No kick-map {coordinate_name} coordinates lie within ±{limit_m:g} m"
        )
    selection = [slice(None)] * np.ndim(values)
    selection[axis] = mask
    return coordinates_m[mask], values[tuple(selection)]


def crop_horizontal_range(
    x_m: np.ndarray, values: np.ndarray, x_range_m: float | None
) -> tuple[np.ndarray, np.ndarray]:
    return _crop_coordinate_range(
        x_m, values, x_range_m, axis=-1, coordinate_name="x"
    )


def crop_vertical_range(
    y_m: np.ndarray, values: np.ndarray, y_range_m: float | None
) -> tuple[np.ndarray, np.ndarray]:
    return _crop_coordinate_range(
        y_m, values, y_range_m, axis=0, coordinate_name="y"
    )

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
    plt.close(figure)

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
    horizontal_range_m: float,
    vertical_range_m: float,
    keep_open: bool = False,
):
    """Plot one representative plane for each enabled array family."""
    families = {strip.family for strip in strips}
    panel_specs: list[tuple[str, str, str, str, float]] = []
    if "horizontal" in families:
        panel_specs.append(
            ("horizontal", "top", "bottom", "strip-centre x [mm]", horizontal_range_m)
        )
    if "vertical" in families:
        panel_specs.append(
            ("vertical", "right", "left", "strip-centre y [mm]", vertical_range_m)
        )
    if not panel_specs:
        raise ValueError("No strip-current family is available for plotting")

    if keep_open:
        figure = plt.figure(figsize=(8.0 if len(panel_specs) == 1 else 13.5, 5.0))
    else:
        figure = Figure(figsize=(8.0 if len(panel_specs) == 1 else 13.5, 5.0))
        FigureCanvasAgg(figure)

    symmetry: dict[str, tuple[float, float]] = {}
    for panel_index, (family, positive, negative, xlabel, limit_m) in enumerate(
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
        axes.set_xlim(-limit_m * 1.0e3, limit_m * 1.0e3)
        axes.grid(alpha=0.3)
        symmetry[family] = (max_mismatch, rms_mismatch)

    figure.suptitle("Fitted strip currents — representative planes")
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)

    if keep_open:
        return figure, symmetry
    plt.close(figure)
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

    plt.close(figure)
    return None


def _add_arguments(
    group: argparse._ArgumentGroup,
    specs: list[tuple[tuple[str, ...], dict[str, Any]]],
) -> None:
    for flags, options in specs:
        group.add_argument(*flags, **options)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fit longitudinal current strips to a 2D EPU kick map",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Help commands:\n"
            "  %(prog)s help                 show all command-line options\n"
            "  %(prog)s help target-mode     show help for one option\n"
            "  %(prog)s --help               standard argparse help"
        ),
    )
    groups = {
        name: parser.add_argument_group(name)
        for name in (
            "input and output",
            "fit objective",
            "strip geometry",
            "model aperture and strip span",
            "acceptance",
            "plotting",
        )
    }
    _add_arguments(groups["input and output"], [
        (("--kick",), dict(type=Path, required=True)),
        (("--out-prefix",), dict(type=Path, default=Path("epu57_fixed"))),
        (("--units",), dict(
            type=parse_kick_units, choices=("microrad", "T^2m^2"),
            default="microrad",
            help="input kick units: microrad (default) or T^2m^2",
        )),
        (("--beam-energy", "--beam_energy", "-beam_energy"), dict(
            dest="beam_energy_ev", type=float, required=True, metavar="EV",
            help="electron total energy in eV",
        )),
        (("--kick-scale",), dict(
            type=float, default=1.0, metavar="FACTOR",
            help="additional dimensionless kick multiplier (default: 1)",
        )),
    ])
    _add_arguments(groups["fit objective"], [
        (("--horizontal-weight",), dict(type=float, default=1.0)),
        (("--vertical-weight",), dict(type=float, default=0.0)),
        (("--target-mode",), dict(
            choices=("nonlinear", "full"), default="nonlinear",
            help="nonlinear removes a+b*x+c*y; full fits the complete 2D maps",
        )),
        (("--ridge",), dict(type=float, default=1.0e-5)),
        (("--limit",), dict(type=float, default=6.0)),
        (("--quadrature-order",), dict(type=int, default=4)),
    ])
    _add_arguments(groups["strip geometry"], [
        (("--strip-layout",), dict(choices=("horizontal", "vertical", "both"), default="horizontal")),
        (("--gap",), dict(type=float, default=0.016)),
        (("--side-gap",), dict(type=float, default=None, metavar="M")),
        (("--length",), dict(type=float, default=None)),
        (("--strips-per-plane", "--horizontal-strips-per-plane"), dict(
            dest="strips_per_plane", type=int, default=21, metavar="N",
        )),
        (("--vertical-strips-per-plane",), dict(type=int, default=None, metavar="N")),
        (("--width",), dict(type=float, default=2.0e-3)),
        (("--thickness",), dict(type=float, default=0.3e-3)),
    ])
    _add_arguments(groups["model aperture and strip span"], [
        (("--x-range",), dict(
            type=float, default=0.020, metavar="M",
            help=("horizontal half-range in metres for both the fitted model "
                  "aperture and the horizontal strip-centre span"),
        )),
        (("--y-range",), dict(
            type=float, default=None, metavar="M",
            help=("vertical half-range in metres for both the fitted model "
                  "aperture and the vertical strip-centre span; default uses "
                  "the full input y extent"),
        )),
    ])
    _add_arguments(groups["acceptance"], [
        (("--max-fit-relative-rms",), dict(type=float, default=None, metavar="FRACTION")),
        (("--min-horizontal-rms-reduction-factor",), dict(type=float, default=1.0)),
        (("--min-vertical-rms-reduction-factor",), dict(type=float, default=1.0)),
        (("--min-full-map-rms-reduction-factor",), dict(type=float, default=1.0)),
        (("--min-full-map-peak-reduction-factor",), dict(type=float, default=0.0)),
        (("--current-bound-tolerance",), dict(type=float, default=1.0e-6)),
    ])
    _add_arguments(groups["plotting"], [
        (("--color-scale",), dict(choices=("auto", "shared"), default="auto")),
        (("--elevation-deg",), dict(type=float, default=PLOT_ELEVATION_DEG)),
        (("--azimuth-deg",), dict(type=float, default=PLOT_AZIMUTH_DEG)),
        (("--projection",), dict(choices=("ortho", "persp"), default="ortho")),
        (("--box-aspect-x",), dict(type=float, default=2.5)),
        (("--box-aspect-y",), dict(type=float, default=1.8)),
        (("--box-aspect-z",), dict(type=float, default=0.85)),
        (("--show",), dict(action="store_true")),
    ])
    return parser

def _option_action(
    parser: argparse.ArgumentParser, topic: str
) -> argparse.Action | None:
    """Resolve a help topic to one parser option action."""
    normalized = topic.strip().lstrip("-").replace("_", "-")
    for action in parser._actions:
        for option in action.option_strings:
            if option.lstrip("-").replace("_", "-") == normalized:
                return action
    return None


def print_option_help(parser: argparse.ArgumentParser, topic: str) -> None:
    """Print focused help for one command-line option."""
    action = _option_action(parser, topic)
    if action is None:
        parser.error(
            f"unknown help topic {topic!r}; use '{parser.prog} help' "
            "to list all options"
        )

    option_names = ", ".join(action.option_strings)
    metavar = action.metavar
    if metavar is None:
        if action.choices is not None:
            metavar = "{" + ",".join(str(value) for value in action.choices) + "}"
        elif action.nargs == 0:
            metavar = ""
        else:
            metavar = action.dest.upper()
    heading = f"{option_names} {metavar}".rstrip()

    print(heading)
    print("-" * len(heading))
    print(action.help or "No additional help is available for this option.")
    if action.choices is not None:
        print("choices: " + ", ".join(str(value) for value in action.choices))
    if action.default is not argparse.SUPPRESS:
        print(f"default: {action.default}")
    if getattr(action, "required", False):
        print("required: yes")


def parse_command_line(
    argv: list[str] | None = None,
) -> argparse.Namespace:
    """Parse normal options or handle the dedicated ``help`` command."""
    parser = build_parser()
    arguments = list(sys.argv[1:] if argv is None else argv)

    if arguments and arguments[0] == "help":
        if len(arguments) == 1:
            parser.print_help()
        elif len(arguments) == 2:
            print_option_help(parser, arguments[1])
        else:
            parser.error("help accepts at most one option name")
        raise SystemExit(0)

    return parser.parse_args(arguments)

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


def _validate_fields(
    args: argparse.Namespace,
    names: tuple[str, ...],
    **constraints: bool,
) -> None:
    for name in names:
        _require_finite(
            name.replace("_", "-"), getattr(args, name), **constraints
        )


def validate_arguments(args: argparse.Namespace) -> argparse.Namespace:
    _validate_fields(
        args, ("kick_scale", "x_range", "gap", "width", "thickness", "limit"),
        positive=True,
    )
    _validate_fields(args, ("y_range", "side_gap", "length"), positive=True, allow_none=True)
    _validate_fields(args, ("ridge",), nonnegative=True)
    _validate_fields(args, ("horizontal_weight", "vertical_weight"), nonnegative=True)
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

    _validate_fields(args, ("max_fit_relative_rms",), nonnegative=True, allow_none=True)
    args.effective_max_fit_relative_rms = (
        0.15 if args.target_mode == "nonlinear" else 0.05
    ) if args.max_fit_relative_rms is None else args.max_fit_relative_rms
    args.max_fit_relative_rms_source = (
        "target_mode_default" if args.max_fit_relative_rms is None else "user"
    )

    _validate_fields(args, (
        "min_horizontal_rms_reduction_factor",
        "min_vertical_rms_reduction_factor",
        "min_full_map_rms_reduction_factor",
    ), positive=True)
    _validate_fields(args, (
        "min_full_map_peak_reduction_factor", "current_bound_tolerance"
    ), nonnegative=True)
    if args.current_bound_tolerance >= 1.0:
        raise ValueError("current-bound-tolerance must lie in [0, 1)")

    _validate_fields(args, ("elevation_deg", "azimuth_deg"))
    _validate_fields(args, ("box_aspect_x", "box_aspect_y", "box_aspect_z"), positive=True)
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


def _kv_text(items: dict[str, Any]) -> str:
    return "".join(f"{key}={_format_value(value)}\n" for key, value in items.items())


class CurrentStripWorkflow:
    """Run one correction while keeping the two kick components symmetric."""

    OUTPUT_SUFFIXES = {
        "currents": "_strip_currents.csv",
        "currents_plot": "_strip_currents_one_plane.png",
        "corrected_map": "_corrected_kickmap.dat",
        "metrics": "_metrics.txt",
        "manifest": "_plot_manifest.txt",
        "theta_x_cut": "_theta_x_midplane_cut.png",
        "theta_y_cut": "_theta_y_centerline_cut.png",
        "raw_3d": "_raw_kickmaps_3d.png",
        "corrected_3d": "_corrected_kickmaps_3d.png",
    }

    def __init__(self, args: argparse.Namespace):
        self.args = validate_arguments(args)
        self.metrics: dict[str, Any] = {}
        self.acceptance_failures: list[str] = []

    @property
    def components(self) -> tuple[ComponentState, ComponentState]:
        return self.theta_x, self.theta_y

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
        manifest_text = self._write_outputs()
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
        self.theta_x = ComponentState(
            key="x",
            label="horizontal",
            raw_map=self.kick_map.theta_x_urad,
            response_component="Iy",
            target_sign=+1.0,
            correction_sign=-1.0,
            weight=a.horizontal_weight,
        )
        self.theta_y = ComponentState(
            key="y",
            label="vertical",
            raw_map=self.kick_map.theta_y_urad,
            response_component="Ix",
            target_sign=-1.0,
            correction_sign=+1.0,
            weight=a.vertical_weight,
        )

        maximum_map_y_m = float(np.max(np.abs(self.kick_map.y_m)))
        self.model_y_half_span_m = (
            maximum_map_y_m if a.y_range is None else a.y_range
        )
        self.model_y_range_source = (
            "complete_input_y_grid" if a.y_range is None else "user"
        )
        if (
            a.strip_layout in {"horizontal", "both"}
            and a.gap / 2.0 <= self.model_y_half_span_m
        ):
            raise ValueError("The top/bottom gap must enclose the model y aperture")

        if a.side_gap is None:
            self.effective_side_gap_m = 2.0 * (a.x_range + a.width)
            self.side_gap_source = "auto_model_x_range_plus_one_width_margin"
        else:
            self.effective_side_gap_m = a.side_gap
            self.side_gap_source = "user"
        if (
            a.strip_layout in {"vertical", "both"}
            and self.effective_side_gap_m / 2.0 <= a.x_range
        ):
            raise ValueError("The left/right gap must enclose the model x aperture")

        strip_length_m = self.kick_map.length_m if a.length is None else a.length
        self.strips, self.strip_geometry = make_strip_arrays(
            layout=a.strip_layout,
            horizontal_half_span_m=a.x_range,
            horizontal_count_per_plane=a.strips_per_plane,
            vertical_half_span_m=self.model_y_half_span_m,
            vertical_count_per_plane=a.effective_vertical_strips_per_plane,
            width_m=a.width,
            thickness_m=a.thickness,
            top_bottom_clear_gap_m=a.gap,
            left_right_clear_gap_m=self.effective_side_gap_m,
            length_m=strip_length_m,
            side_gap_source=self.side_gap_source,
        )
        self.horizontal_indices = self._strip_indices("horizontal")
        self.vertical_indices = self._strip_indices("vertical")

        self.model_x_m, self.model_y_m, self.theta_x.model_map = _crop_2d(
            self.kick_map.x_m,
            self.kick_map.y_m,
            self.theta_x.raw_map,
            a.x_range,
            a.y_range,
        )
        _, _, self.theta_y.model_map = _crop_2d(
            self.kick_map.x_m,
            self.kick_map.y_m,
            self.theta_y.raw_map,
            a.x_range,
            a.y_range,
        )
        if self.model_x_m.size < 2 or self.model_y_m.size < 2:
            raise ValueError("The model aperture must contain at least two x and y points")

        model_x_grid, model_y_grid = np.meshgrid(self.model_x_m, self.model_y_m)
        self.model_points = np.column_stack(
            [model_x_grid.ravel(), model_y_grid.ravel()]
        )
        weight_scale = max(component.weight for component in self.components)
        self.total_normalized_fit_weight = 0.0
        for component in self.components:
            component.response_model = response_matrix(
                self.strips,
                self.model_points,
                component.response_component,
                a.quadrature_order,
            )
            component.raw_flat = component.model_map.ravel()
            component.target_tm = (
                component.target_sign
                * self.beam_rigidity_tm
                * component.raw_flat
                * 1.0e-6
            )
            component.raw_plane = affine_plane_decomposition(
                self.model_points, component.raw_flat
            )
            target_plane = affine_plane_decomposition(
                self.model_points, component.target_tm
            )
            if a.target_mode == "nonlinear":
                component.fit_target_tm = target_plane.residual
                component.fit_basis = "full_2d_constant_and_linear_plane_removed"
            else:
                component.fit_target_tm = component.target_tm
                component.fit_basis = "full_2d_complete_map"
            component.normalized_weight = component.weight / weight_scale
            self.total_normalized_fit_weight += component.normalized_weight

    def _strip_indices(self, family: str) -> np.ndarray:
        return np.asarray(
            [i for i, strip in enumerate(self.strips) if strip.family == family],
            dtype=int,
        )

    def _solve(self) -> None:
        a = self.args
        response_blocks: list[np.ndarray] = []
        target_blocks: list[np.ndarray] = []
        for component in self.components:
            if component.normalized_weight > 0.0:
                scale = np.sqrt(
                    component.normalized_weight / component.fit_target_tm.size
                )
                response_blocks.append(scale * component.response_model)
                target_blocks.append(scale * component.fit_target_tm)
        self.current_solution = solve_currents(
            np.vstack(response_blocks),
            np.concatenate(target_blocks),
            a.ridge,
            a.limit,
        )
        self.currents_a = self.current_solution.currents_a

        x_grid, y_grid = np.meshgrid(self.kick_map.x_m, self.kick_map.y_m)
        all_points = np.column_stack([x_grid.ravel(), y_grid.ravel()])
        for component in self.components:
            response_full = response_matrix(
                self.strips,
                all_points,
                component.response_component,
                a.quadrature_order,
            )
            strip_field_tm = (response_full @ self.currents_a).reshape(x_grid.shape)
            component.corrected_map = component.raw_map + (
                component.correction_sign
                * 1.0e6
                * strip_field_tm
                / self.beam_rigidity_tm
            )

    def _model_cut(self, values: np.ndarray, component: ComponentState) -> ZeroCut:
        a = self.args
        if component.key == "x":
            complete = zero_coordinate_cut(self.kick_map.y_m, values, 0, "y")
            _, cropped = crop_horizontal_range(
                self.kick_map.x_m, complete.values, a.x_range
            )
        else:
            complete = zero_coordinate_cut(self.kick_map.x_m, values, 1, "x")
            _, column = crop_vertical_range(
                self.kick_map.y_m, complete.values[:, np.newaxis], a.y_range
            )
            cropped = column[:, 0]
        return ZeroCut(
            cropped,
            complete.method,
            complete.lower_coordinate_m,
            complete.upper_coordinate_m,
        )

    def _analyse(self) -> None:
        a = self.args
        for component in self.components:
            _, _, component.corrected_model = _crop_2d(
                self.kick_map.x_m,
                self.kick_map.y_m,
                component.corrected_map,
                a.x_range,
                a.y_range,
            )
            component.corrected_flat = component.corrected_model.ravel()
            component.corrected_plane = affine_plane_decomposition(
                self.model_points, component.corrected_flat
            )
            component.raw_cut = self._model_cut(component.raw_map, component)
            component.corrected_cut = self._model_cut(
                component.corrected_map, component
            )
            coordinate = self.model_x_m if component.key == "x" else self.model_y_m
            component.raw_affine = affine_decomposition(
                coordinate, component.raw_cut.values
            )
            component.corrected_affine = affine_decomposition(
                coordinate, component.corrected_cut.values
            )

            component.fit_prediction_tm = (
                component.response_model @ self.currents_a
            )
            component.fit_residual_tm = (
                component.fit_prediction_tm - component.fit_target_tm
            )
            full_residual = (
                component.response_model @ self.currents_a - component.target_tm
            )
            component.fit_target_rms_tm = root_mean_square(
                component.fit_target_tm
            )
            component.fit_residual_rms_tm = root_mean_square(
                component.fit_residual_tm
            )
            component.fit_relative_rms = relative_rms_residual(
                component.fit_residual_tm, component.fit_target_tm
            )
            component.full_fit_relative_rms = relative_rms_residual(
                full_residual, component.target_tm
            )

        self.fit_target_rms_tm = self._weighted_rms(
            self.theta_x.fit_target_rms_tm,
            self.theta_y.fit_target_rms_tm,
        )
        self.fit_residual_rms_tm = self._weighted_rms(
            self.theta_x.fit_residual_rms_tm,
            self.theta_y.fit_residual_rms_tm,
        )
        self.fit_relative_rms = (
            0.0
            if self.fit_target_rms_tm == self.fit_residual_rms_tm == 0.0
            else float("inf")
            if self.fit_target_rms_tm == 0.0
            else self.fit_residual_rms_tm / self.fit_target_rms_tm
        )
        self._compute_map_metrics()
        self._compute_centerline_metrics()
        self._evaluate_acceptance()

    def _weighted_rms(self, horizontal: float, vertical: float) -> float:
        return float(
            np.sqrt(
                (
                    self.theta_x.normalized_weight * horizontal**2
                    + self.theta_y.normalized_weight * vertical**2
                )
                / self.total_normalized_fit_weight
            )
        )

    @staticmethod
    def _set_reduction(
        owner: Any, prefix: str, raw: float, corrected: float
    ) -> None:
        factor, status = safe_reduction_factor(raw, corrected)
        setattr(owner, f"{prefix}_factor", factor)
        setattr(owner, f"{prefix}_status", status)

    def _compute_map_metrics(self) -> None:
        for component in self.components:
            component.raw_model_rms = root_mean_square(component.raw_flat)
            component.corrected_model_rms = root_mean_square(
                component.corrected_flat
            )
            self._set_reduction(
                component,
                "map_rms_reduction",
                component.raw_model_rms,
                component.corrected_model_rms,
            )
            component.raw_model_peak = float(np.max(np.abs(component.raw_flat)))
            component.corrected_model_peak = float(
                np.max(np.abs(component.corrected_flat))
            )
            self._set_reduction(
                component,
                "map_peak_reduction",
                component.raw_model_peak,
                component.corrected_model_peak,
            )
            component.raw_nonlinear_rms = root_mean_square(
                component.raw_plane.residual
            )
            component.corrected_nonlinear_rms = root_mean_square(
                component.corrected_plane.residual
            )
            self._set_reduction(
                component,
                "nonlinear_rms_reduction",
                component.raw_nonlinear_rms,
                component.corrected_nonlinear_rms,
            )

    def _compute_centerline_metrics(self) -> None:
        x = self.theta_x
        x.raw_centerline_rms = root_mean_square(x.raw_cut.values)
        x.corrected_centerline_rms = root_mean_square(x.corrected_cut.values)
        self._set_reduction(
            x,
            "centerline_rms_reduction",
            x.raw_centerline_rms,
            x.corrected_centerline_rms,
        )
        x.raw_centerline_nonlinear_rms = root_mean_square(x.raw_affine.residual)
        x.corrected_centerline_nonlinear_rms = root_mean_square(
            x.corrected_affine.residual
        )
        self._set_reduction(
            x,
            "centerline_nonlinear_rms_reduction",
            x.raw_centerline_nonlinear_rms,
            x.corrected_centerline_nonlinear_rms,
        )

        y = self.theta_y
        y.raw_centerline_peak = float(np.max(np.abs(y.raw_cut.values)))
        y.corrected_centerline_peak = float(np.max(np.abs(y.corrected_cut.values)))
        self._set_reduction(
            y,
            "centerline_peak_reduction",
            y.raw_centerline_peak,
            y.corrected_centerline_peak,
        )
        y.raw_centerline_nonlinear_peak = float(
            np.max(np.abs(y.raw_affine.residual))
        )
        y.corrected_centerline_nonlinear_peak = float(
            np.max(np.abs(y.corrected_affine.residual))
        )
        self._set_reduction(
            y,
            "centerline_nonlinear_peak_reduction",
            y.raw_centerline_nonlinear_peak,
            y.corrected_centerline_nonlinear_peak,
        )

    def _evaluate_acceptance(self) -> None:
        a = self.args
        for component in self.components:
            if a.target_mode == "nonlinear":
                component.acceptance_reduction_factor = (
                    component.nonlinear_rms_reduction_factor
                )
                component.acceptance_reduction_status = (
                    component.nonlinear_rms_reduction_status
                )
                component.acceptance_rms_basis = "full_2d_affine_plane_removed"
            else:
                component.acceptance_reduction_factor = (
                    component.map_rms_reduction_factor
                )
                component.acceptance_reduction_status = (
                    component.map_rms_reduction_status
                )
                component.acceptance_rms_basis = "full_2d_complete_map"

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
        component_minima = (
            (self.theta_x, a.min_horizontal_rms_reduction_factor),
            (self.theta_y, a.min_vertical_rms_reduction_factor),
        )
        for component, minimum in component_minima:
            if (
                component.weight > 0.0
                and component.acceptance_reduction_factor < minimum
            ):
                self.acceptance_failures.append(
                    f"{component.label} active-basis 2D RMS reduction factor "
                    f"{component.acceptance_reduction_factor:.6g} is below "
                    f"{minimum:.6g}"
                )
        for component in self.components:
            if component.map_rms_reduction_factor < a.min_full_map_rms_reduction_factor:
                self.acceptance_failures.append(
                    f"{component.label} full-map RMS reduction factor "
                    f"{component.map_rms_reduction_factor:.6g} is below "
                    f"{a.min_full_map_rms_reduction_factor:.6g}"
                )
        if a.min_full_map_peak_reduction_factor > 0.0:
            for component in self.components:
                if (
                    component.map_peak_reduction_factor
                    < a.min_full_map_peak_reduction_factor
                ):
                    self.acceptance_failures.append(
                        f"{component.label} full-map peak reduction factor "
                        f"{component.map_peak_reduction_factor:.6g} is below "
                        f"{a.min_full_map_peak_reduction_factor:.6g}"
                    )
        self.acceptance_status = (
            "FAILED"
            if self.acceptance_failures
            else "ACCEPTED_LIMIT_SATURATED"
            if self.current_bound_active_count
            else "ACCEPTED"
        )

    def _build_metrics(self) -> None:
        a, g, s = self.args, self.strip_geometry, self.current_solution
        x, y = self.components
        m = self.metrics
        m.update([
            ("release_id", RELEASE_ID),
            ("input_file", a.kick),
            ("grid_shape", self.kick_map.theta_x_urad.shape),
            ("model_domain_source", "x_range_and_y_range"),
            ("model_x_requested_half_span_m", a.x_range),
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
            ("input_kick_units", a.units),
            ("beam_energy_eV", a.beam_energy_ev),
            ("electron_rest_energy_eV", ELECTRON_REST_ENERGY_EV),
            ("beam_momentum_eV_c", self.beam_momentum_ev_c),
            ("beam_rigidity_Tm", self.beam_rigidity_tm),
            ("input_to_urad_scale", self.input_to_urad_scale),
            ("kick_scale", a.kick_scale),
            ("horizontal_weight", a.horizontal_weight),
            ("vertical_weight", a.vertical_weight),
            ("target_mode", a.target_mode),
            ("horizontal_fit_basis", x.fit_basis),
            ("vertical_fit_basis", y.fit_basis),
            ("fit_domain", "full_2d_model_aperture"),
            ("fit_point_count", self.model_points.shape[0]),
            ("affine_projection_basis_2d", "constant_plus_x_plus_y"),
            ("normalized_horizontal_weight", x.normalized_weight),
            ("normalized_vertical_weight", y.normalized_weight),
            ("fit_weighting", "max_normalized_component_weight_times_mean_squared_residual"),
            ("color_scale", a.color_scale),
            ("plot_elevation_deg", a.elevation_deg),
            ("plot_azimuth_deg", a.azimuth_deg),
            ("projection", a.projection),
            ("box_aspect_x", a.box_aspect_x),
            ("box_aspect_y", a.box_aspect_y),
            ("box_aspect_z", a.box_aspect_z),
        ])

        for prefix, cut in (("theta_x_y0", x.raw_cut), ("theta_y_x0", y.raw_cut)):
            m.update([
                (f"{prefix}_cut_method", cut.method),
                (f"{prefix}_bracket_lower_m", cut.lower_coordinate_m),
                (f"{prefix}_bracket_upper_m", cut.upper_coordinate_m),
            ])
        m.update([
            ("raw_theta_x_y0_min_urad", np.min(x.raw_cut.values)),
            ("raw_theta_x_y0_max_urad", np.max(x.raw_cut.values)),
            ("raw_theta_y_x0_min_urad", np.min(y.raw_cut.values)),
            ("raw_theta_y_x0_max_urad", np.max(y.raw_cut.values)),
        ])

        for prefix, decomposition in (
            ("raw_theta_x_affine", x.raw_affine),
            ("corrected_theta_x_affine", x.corrected_affine),
            ("raw_theta_y_affine", y.raw_affine),
            ("corrected_theta_y_affine", y.corrected_affine),
        ):
            m.update([
                (f"{prefix}_offset_urad", float(decomposition.offset)),
                (f"{prefix}_gradient_urad_per_mm", float(decomposition.slope_per_m) * 1.0e-3),
            ])

        m.update([
            ("raw_theta_x_y0_rms_urad", x.raw_centerline_rms),
            ("corrected_theta_x_y0_rms_urad", x.corrected_centerline_rms),
            ("horizontal_rms_reduction_factor", x.centerline_rms_reduction_factor),
            ("horizontal_rms_reduction_status", x.centerline_rms_reduction_status),
            ("raw_theta_x_y0_nonlinear_rms_urad", x.raw_centerline_nonlinear_rms),
            ("corrected_theta_x_y0_nonlinear_rms_urad", x.corrected_centerline_nonlinear_rms),
            ("horizontal_y0_nonlinear_rms_reduction_factor", x.centerline_nonlinear_rms_reduction_factor),
            ("horizontal_y0_nonlinear_rms_reduction_status", x.centerline_nonlinear_rms_reduction_status),
            ("horizontal_acceptance_rms_basis", x.acceptance_rms_basis),
            ("horizontal_acceptance_reduction_factor", x.acceptance_reduction_factor),
            ("horizontal_acceptance_reduction_status", x.acceptance_reduction_status),
            ("vertical_acceptance_rms_basis", y.acceptance_rms_basis),
            ("vertical_acceptance_reduction_factor", y.acceptance_reduction_factor),
            ("vertical_acceptance_reduction_status", y.acceptance_reduction_status),
            ("raw_theta_y_x0_peak_urad", y.raw_centerline_peak),
            ("corrected_theta_y_x0_peak_urad", y.corrected_centerline_peak),
            ("vertical_peak_reduction_factor", y.centerline_peak_reduction_factor),
            ("vertical_peak_reduction_status", y.centerline_peak_reduction_status),
            ("raw_theta_y_x0_nonlinear_peak_urad", y.raw_centerline_nonlinear_peak),
            ("corrected_theta_y_x0_nonlinear_peak_urad", y.corrected_centerline_nonlinear_peak),
            ("vertical_nonlinear_peak_reduction_factor", y.centerline_nonlinear_peak_reduction_factor),
            ("vertical_nonlinear_peak_reduction_status", y.centerline_nonlinear_peak_reduction_status),
        ])

        for measure, raw_attr, corrected_attr, factor_attr in (
            ("rms", "raw_model_rms", "corrected_model_rms", "map_rms_reduction"),
            ("peak", "raw_model_peak", "corrected_model_peak", "map_peak_reduction"),
        ):
            for axis, component in (("x", x), ("y", y)):
                m.update([
                    (f"raw_theta_{axis}_model_2d_{measure}_urad", getattr(component, raw_attr)),
                    (f"corrected_theta_{axis}_model_2d_{measure}_urad", getattr(component, corrected_attr)),
                    (f"full_map_theta_{axis}_{measure}_reduction_factor", getattr(component, f"{factor_attr}_factor")),
                    (f"full_map_theta_{axis}_{measure}_reduction_status", getattr(component, f"{factor_attr}_status")),
                ])
        for axis, component in (("x", x), ("y", y)):
            m.update([
                (f"raw_theta_{axis}_model_2d_nonlinear_rms_urad", component.raw_nonlinear_rms),
                (f"corrected_theta_{axis}_model_2d_nonlinear_rms_urad", component.corrected_nonlinear_rms),
            ])

        for axis, component in (("x", x), ("y", y)):
            for label, plane in (("raw", component.raw_plane), ("corrected", component.corrected_plane)):
                prefix = f"{label}_theta_{axis}_affine_plane"
                m.update([
                    (f"{prefix}_offset_urad", float(plane.offset)),
                    (f"{prefix}_gradient_x_urad_per_mm", float(plane.slope_x_per_m) * 1.0e-3),
                    (f"{prefix}_gradient_y_urad_per_mm", float(plane.slope_y_per_m) * 1.0e-3),
                ])

        m.update([
            ("solver_status", s.status),
            ("solver_message", s.message),
            ("solver_cost", s.cost),
            ("solver_optimality", s.optimality),
            ("solver_iterations", s.iterations),
            ("horizontal_fit_target_rms_Tm", x.fit_target_rms_tm),
            ("horizontal_fit_residual_rms_Tm", x.fit_residual_rms_tm),
            ("horizontal_fit_relative_rms", x.fit_relative_rms),
            ("vertical_fit_target_rms_Tm", y.fit_target_rms_tm),
            ("vertical_fit_residual_rms_Tm", y.fit_residual_rms_tm),
            ("vertical_fit_relative_rms", y.fit_relative_rms),
            ("horizontal_full_fit_relative_rms", x.full_fit_relative_rms),
            ("vertical_full_fit_relative_rms", y.full_fit_relative_rms),
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
            ("strip_layout", g.layout),
            ("horizontal_array_enabled", int(self.horizontal_indices.size > 0)),
            ("vertical_array_enabled", int(self.vertical_indices.size > 0)),
            ("horizontal_strips_per_plane", g.horizontal_count_per_plane),
            ("strip_span_source", "model_aperture_ranges"),
            ("vertical_strips_per_plane", g.vertical_count_per_plane),
            ("horizontal_strip_half_span_m", g.horizontal_half_span_m),
            ("vertical_strip_half_span_m", g.vertical_half_span_m),
            ("horizontal_strip_pitch_m", g.horizontal_pitch_m),
            ("vertical_strip_pitch_m", g.vertical_pitch_m),
            ("top_bottom_clear_gap_m", g.top_bottom_clear_gap_m),
            ("left_right_clear_gap_m", g.left_right_clear_gap_m),
            ("side_gap_source", g.side_gap_source),
            ("strip_width_m", a.width),
            ("strip_thickness_m", a.thickness),
            ("total_strip_count", len(self.strips)),
            ("max_abs_current_A", np.max(np.abs(self.currents_a))),
            ("max_abs_horizontal_array_current_A", self._family_max(self.horizontal_indices)),
            ("max_abs_vertical_array_current_A", self._family_max(self.vertical_indices)),
        ])

    def _family_max(self, indices: np.ndarray) -> float:
        return (
            float(np.max(np.abs(self.currents_a[indices])))
            if indices.size
            else float("nan")
        )

    def _output_paths(self, prefix: Path) -> dict[str, Path]:
        return {
            name: prefix.with_name(prefix.name + suffix)
            for name, suffix in self.OUTPUT_SUFFIXES.items()
        }

    def _write_outputs(self) -> str:
        a = self.args
        x, y = self.components
        prefix = a.out_prefix
        prefix.parent.mkdir(parents=True, exist_ok=True)
        final_paths = self._output_paths(prefix)
        with tempfile.TemporaryDirectory(
            dir=prefix.parent, prefix=f".{prefix.name}_staging_"
        ) as staging_directory:
            paths = self._output_paths(Path(staging_directory) / prefix.name)
            write_strip_currents_csv(paths["currents"], self.strips, self.currents_a)
            currents_figure, symmetry = plot_strip_currents_representative_planes(
                self.strips,
                self.currents_a,
                paths["currents_plot"],
                a.x_range,
                self.model_y_half_span_m,
                keep_open=a.show,
            )
            write_corrected_kick_map(
                paths["corrected_map"],
                self.kick_map,
                x.corrected_map,
                y.corrected_map,
                source_units=a.units,
                beam_energy_ev=a.beam_energy_ev,
                beam_rigidity_tm=self.beam_rigidity_tm,
                kick_scale=a.kick_scale,
            )

            if a.target_mode == "nonlinear":
                cut_values = (
                    x.raw_affine.residual,
                    x.corrected_affine.residual,
                    y.raw_affine.residual,
                    y.corrected_affine.residual,
                )
                titles = (
                    "Horizontal nonlinear centreline diagnostic on y = 0 "
                    "(constant and gradient removed)",
                    "Vertical nonlinear centreline diagnostic on x = 0 "
                    "(constant and gradient removed)",
                )
            else:
                cut_values = (
                    x.raw_cut.values,
                    x.corrected_cut.values,
                    y.raw_cut.values,
                    y.corrected_cut.values,
                )
                titles = (
                    "Horizontal kick on the horizontal mid-plane (y = 0)",
                    "Vertical kick on the vertical centreline (x = 0)",
                )
            for coordinate, raw_values, corrected_values, title, xlabel, ylabel, path in (
                (self.model_x_m * 1.0e3, *cut_values[:2], titles[0], "x [mm]", "θx [µrad]", paths["theta_x_cut"]),
                (self.model_y_m * 1.0e3, *cut_values[2:], titles[1], "y [mm]", "θy [µrad]", paths["theta_y_cut"]),
            ):
                plot_cut(
                    coordinate, raw_values, corrected_values,
                    title, xlabel, ylabel, path,
                )

            shared_limits = (
                max(float(np.max(np.abs(x.model_map))), float(np.max(np.abs(x.corrected_model)))),
                max(float(np.max(np.abs(y.model_map))), float(np.max(np.abs(y.corrected_model)))),
            )
            plot_args = dict(
                x_range_m=a.x_range,
                y_range_m=a.y_range,
                color_scale=a.color_scale,
                elevation_deg=a.elevation_deg,
                azimuth_deg=a.azimuth_deg,
                projection=a.projection,
                box_aspect_x=a.box_aspect_x,
                box_aspect_y=a.box_aspect_y,
                box_aspect_z=a.box_aspect_z,
                keep_open=a.show,
            )
            map_figures = []
            for name, theta_x, theta_y, title in (
                ("raw_3d", x.raw_map, y.raw_map, "Kick maps before current-strip correction"),
                ("corrected_3d", x.corrected_map, y.corrected_map, "Kick maps after current-strip correction"),
            ):
                map_figures.append(plot_kickmaps_3d_pair(
                    self.kick_map.x_m,
                    self.kick_map.y_m,
                    theta_x,
                    theta_y,
                    paths[name],
                    *shared_limits,
                    title,
                    **plot_args,
                ))

            for family in ("horizontal", "vertical"):
                max_mismatch, rms_mismatch = symmetry.get(
                    family, (float("nan"), float("nan"))
                )
                self.metrics[f"{family}_sheet_symmetry_max_abs_mismatch_A"] = max_mismatch
                self.metrics[f"{family}_sheet_symmetry_rms_mismatch_A"] = rms_mismatch
            paths["metrics"].write_text(
                _kv_text(self.metrics), encoding="utf-8", newline="\n"
            )
            manifest_text = _kv_text(self._build_manifest())
            paths["manifest"].write_text(
                manifest_text, encoding="utf-8", newline="\n"
            )
            if a.show:
                plt.show()
                for figure in (currents_figure, *map_figures):
                    if figure is not None:
                        plt.close(figure)
            publish_staged_outputs(
                list(paths.values()), list(final_paths.values())
            )
        return manifest_text

    def _build_manifest(self) -> dict[str, Any]:
        a, g = self.args, self.strip_geometry
        x, y = self.components
        return dict([
                ("release_id", RELEASE_ID),
                ("acceptance_status", self.acceptance_status),
                ("model_x_half_span_m", a.x_range),
                ("model_y_half_span_m", self.model_y_half_span_m),
                ("model_x_point_count", self.model_x_m.size),
                ("model_y_point_count", self.model_y_m.size),
                ("export_grid_scope", "complete_input_grid"),
                ("input_kick_units", a.units),
                ("beam_energy_eV", a.beam_energy_ev),
                ("beam_rigidity_Tm", self.beam_rigidity_tm),
                ("input_to_urad_scale", self.input_to_urad_scale),
                ("kick_scale", a.kick_scale),
                ("horizontal_weight", a.horizontal_weight),
                ("vertical_weight", a.vertical_weight),
                ("target_mode", a.target_mode),
                ("fit_domain", "full_2d_model_aperture"),
                ("fit_point_count", self.model_points.shape[0]),
                ("horizontal_fit_basis", x.fit_basis),
                ("vertical_fit_basis", y.fit_basis),
                ("normalized_horizontal_weight", x.normalized_weight),
                ("normalized_vertical_weight", y.normalized_weight),
                ("strip_layout", g.layout),
                (
                    "strip_currents_representative_planes_file",
                    f"{a.out_prefix.name}_strip_currents_one_plane.png",
                ),
                ("horizontal_strip_half_span_m", g.horizontal_half_span_m),
                ("vertical_strip_half_span_m", g.vertical_half_span_m),
                ("horizontal_strip_pitch_m", g.horizontal_pitch_m),
                ("vertical_strip_pitch_m", g.vertical_pitch_m),
                ("horizontal_strips_per_plane", g.horizontal_count_per_plane),
                ("vertical_strips_per_plane", g.vertical_count_per_plane),
                ("top_bottom_clear_gap_m", g.top_bottom_clear_gap_m),
                ("left_right_clear_gap_m", g.left_right_clear_gap_m),
                ("raw_kickmaps_file", f"{a.out_prefix.name}_raw_kickmaps_3d.png"),
                (
                    "corrected_kickmaps_file",
                    f"{a.out_prefix.name}_corrected_kickmaps_3d.png",
                ),
                ("plot_elevation_deg", a.elevation_deg),
                ("plot_azimuth_deg", a.azimuth_deg),
        ])

def main() -> None:
    CurrentStripWorkflow(parse_command_line()).run()


if __name__ == "__main__":
    main()
