#!/usr/bin/env python3
"""
current_strip_epu_57_solver_16.py

Self-contained current-strip correction prototype for the supplied EPU57
RADIA kick map.

Key corrections relative to current_strip_epu_57_solver_2.py
-------------------------------------------------------------
1. The vertical kick is no longer sampled on y = 0, where it vanishes by
   symmetry.  Vertical diagnostics use the full 2D map and the x = 0 cut.
2. One physical strip-current vector is fitted to a configurable weighted
   combination of the horizontal y = 0 and vertical x = 0 targets, and is
   then applied to BOTH kick components.
3. Right-handed Lorentz-force signs are applied consistently:
       theta_x(strip) = - integral(B_y ds) / (B rho)
       theta_y(strip) = + integral(B_x ds) / (B rho)
4. The complete corrected 2D kick map is exported.
5. Horizontal and vertical kick maps are combined side-by-side in one raw
   3D figure and one corrected 3D figure; --show displays both figures.
6. Filled-contour 2D kick-map PNGs are not generated. The two 1D diagnostic
   cuts are retained.
7. --x-range-mm sets both the displayed horizontal range and the outermost
   strip-centre positions. --strips-per-plane sets the number of strips on
   each sheet (default 21). The strip pitch is derived uniformly from -X to
   +X. With X=20 mm and 21 strips, the pitch is 2.0 mm. --y-range-mm limits
   only the displayed vertical range. The exported kick map remains full-range.
8. --kick-scale multiplies both input kick components before fitting,
   correction, plotting, metrics, and export; its default is 1.0.
9. The 3D colour and z ranges autoscale independently from the displayed
   data by default. --color-scale shared restores the earlier common,
   symmetric raw/corrected component scales.
10. Both panels in both paired figures use one shared oblique camera. The
    defaults are elevation 20 degrees and azimuth -45 degrees, which keep
    both x and y directions visible. --elevation-deg and --azimuth-deg
    override these values and are applied identically to every panel.
11. Before every 3D plot, x and y coordinates are sorted into ascending
    order and the kick matrix is reordered consistently.  This is essential
    because the RADIA file stores y from positive to negative values.
12. Each paired raw/corrected plot is rendered from a fresh, isolated
    Figure object in batch mode. The horizontal and vertical panels are
    independent Axes3D objects inside that figure and share only the explicit
    camera and geometry settings.
13. --projection and --box-aspect-x/y/z control the 3D geometry.
    The defaults use orthographic projection and box aspect
    (2.5, 1.8, 0.85). The larger y aspect reflects Gnuplot's independent
    axis normalization and prevents the theta_y surface from appearing
    artificially edge-on.
14. The fitted current distribution is plotted for the upper strip plane.
    Since the two planes are symmetric, the lower-plane distribution is not
    plotted separately; the row-to-row mismatch is recorded in the metrics.
15. The strip-current plot uses the same ±x-range-mm strip-centre span.
    With the default 2 mm strip width, 21 strips at 2 mm pitch touch without
    overlap and their physical outer edges lie at ±21 mm.
16. The script has no dependency on the missing
   current_strip_poisson_solver.py module.
17. Horizontal and vertical zero-coordinate diagnostics are now evaluated
    at exact x = 0 or y = 0 coordinates.  When the input grid does not
    contain zero, the kick map is linearly interpolated between the two grid
    coordinates that bracket zero; extrapolation is rejected.
18. Solver convergence is separated from physical acceptance.  The script
    reports fit residuals, current-bound activity, and horizontal RMS
    reduction, and it rejects solutions that do not satisfy configurable
    acceptance thresholds.
19. Reduction factors are safe for exact-zero residuals.  All output files
    are staged in a temporary directory and published only after the full run
    and acceptance checks complete successfully.
20. --horizontal-weight and --vertical-weight set independent non-negative
    relative component weights in the joint least-squares objective.  The
    weights are normalized by their maximum and each component is normalized
    by its number of sampled cut points.  Defaults remain (1, 0).
21. --target-mode nonlinear (the default) projects the constant and linear
    components out of both each target cut and every strip-response column
    before fitting.  The optimizer therefore corrects only the nonlinear
    residuals of theta_x(x, y=0) and theta_y(x=0, y).  The exported map is
    still the complete physical map; its remaining affine coefficients are
    recorded for subsequent lattice correction.  --target-mode full restores
    the R15 full-kick objective.

Model limitations
-----------------
The strip response is the free-space 2D Green-function solution for long
longitudinal rectangular current strips, evaluated by Gauss-Legendre
quadrature over each strip cross-section.  This is appropriate for a local
proof-of-principle model, but it omits return conductors, end turns, magnetic
materials, chamber effects, and finite-length end fields.

The input kick-map coordinates are assumed to be the standard right-handed
accelerator coordinates x, y, s with the beam along +s.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib

# Use a non-interactive backend for batch runs. When --show is requested,
# retain the user's normal interactive backend.
if "--show" not in sys.argv:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy as np
from scipy.optimize import lsq_linear

MU0 = 4.0e-7 * np.pi

RELEASE_ID = "EPU57-SOLVER-16-NONLINEAR-TARGET-PROJECTION"
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
    x_m: float
    y_m: float
    width_m: float
    thickness_m: float
    length_m: float


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


def make_two_strip_rows(
    centre_half_span_m: float,
    n_per_row: int,
    width_m: float,
    thickness_m: float,
    gap_m: float,
    length_m: float,
) -> tuple[list[Strip], float]:
    """Build two symmetric strip rows with centres spanning ``±X``.

    The ``n_per_row`` centres are uniformly distributed from
    ``-centre_half_span_m`` to ``+centre_half_span_m``. The derived pitch is
    ``2*X/(n_per_row-1)`` for more than one strip. Strip overlap is rejected;
    touching strips are allowed.
    """
    if n_per_row < 1:
        raise ValueError("strips-per-plane must be positive")
    if min(
        centre_half_span_m,
        width_m,
        thickness_m,
        gap_m,
        length_m,
    ) <= 0.0:
        raise ValueError("All geometric dimensions must be positive")

    if n_per_row == 1:
        x_centres = np.asarray([0.0], dtype=float)
        derived_pitch_m = 0.0
    else:
        x_centres = np.linspace(
            -centre_half_span_m,
            +centre_half_span_m,
            n_per_row,
            dtype=float,
        )
        derived_pitch_m = float(
            2.0 * centre_half_span_m / (n_per_row - 1)
        )
        if width_m > derived_pitch_m + 1.0e-15:
            raise ValueError(
                "strip width exceeds the derived pitch; reduce --width, "
                "reduce --strips-per-plane, or increase --x-range-mm"
            )

    y_centre = gap_m / 2.0 + thickness_m / 2.0
    strips = [
        Strip(
            x_m=float(x_centre),
            y_m=float(row_sign * y_centre),
            width_m=width_m,
            thickness_m=thickness_m,
            length_m=length_m,
        )
        for row_sign in (+1.0, -1.0)
        for x_centre in x_centres
    ]
    return strips, derived_pitch_m

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
        source_x = strip.x_m + 0.5 * strip.width_m * nodes
        source_y = strip.y_m + 0.5 * strip.thickness_m * nodes
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
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# Corrected EPU57 current-strip kick map\n")
        stream.write("# Right-handed x,y,s convention; beam along +s\n")
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
    x_range_mm: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Crop plotting arrays to |x| <= x_range_mm.

    The solver and exported kick map always retain the complete input grid.
    """
    if x_range_mm is None:
        return x_m, values
    if x_range_mm <= 0.0:
        raise ValueError("x_range_mm must be positive")

    limit_m = x_range_mm * 1.0e-3
    mask = np.abs(x_m) <= limit_m + 1.0e-15
    if not np.any(mask):
        raise ValueError(
            f"No kick-map x coordinates lie within ±{x_range_mm:g} mm"
        )
    return x_m[mask], values[..., mask]


def crop_vertical_range(
    y_m: np.ndarray,
    values: np.ndarray,
    y_range_mm: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Crop plotting arrays to |y| <= y_range_mm.

    The solver and exported kick map always retain the complete input grid.
    The first axis of values is assumed to correspond to y.
    """
    if y_range_mm is None:
        return y_m, values
    if y_range_mm <= 0.0:
        raise ValueError("y_range_mm must be positive")

    limit_m = y_range_mm * 1.0e-3
    mask = np.abs(y_m) <= limit_m + 1.0e-15
    if not np.any(mask):
        raise ValueError(
            f"No kick-map y coordinates lie within ±{y_range_mm:g} mm"
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


def plot_map(
    x_m: np.ndarray,
    y_m: np.ndarray,
    values_urad: np.ndarray,
    title: str,
    output_path: Path,
    symmetric_limit_urad: float,
    x_range_mm: float | None = None,
) -> None:
    x_plot_m, values_plot_urad = crop_horizontal_range(
        x_m, values_urad, x_range_mm
    )
    x_grid, y_grid = np.meshgrid(x_plot_m * 1.0e3, y_m * 1.0e3)
    levels = np.linspace(-symmetric_limit_urad, symmetric_limit_urad, 41)

    figure, axes = _make_2d_figure()
    contour = axes.contourf(
        x_grid,
        y_grid,
        values_plot_urad,
        levels=levels,
        extend="both",
    )
    figure.colorbar(contour, ax=axes, label="kick [µrad]")
    axes.set_title(title)
    axes.set_xlabel("x [mm]")
    axes.set_ylabel("y [mm]")
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    try:
        plt.close(figure)
    except Exception:
        pass


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





def plot_strip_currents_one_plane(
    strips: list[Strip],
    currents_a: np.ndarray,
    n_per_row: int,
    output_path: Path,
    x_range_mm: float,
    keep_open: bool = False,
):
    """Plot the fitted current distribution for the upper strip plane.

    ``make_two_strip_rows`` stores the complete upper row first and the
    complete lower row second.  The upper row is plotted against its physical
    strip-centre x coordinate.  The lower row is used only for the symmetry
    check returned to the caller.
    """
    currents = np.asarray(currents_a, dtype=float)
    if n_per_row < 1:
        raise ValueError("n_per_row must be positive")
    if len(strips) != 2 * n_per_row:
        raise ValueError("Expected exactly two strip rows")
    if currents.shape != (2 * n_per_row,):
        raise ValueError("Current vector does not match the two strip rows")

    upper_strips = strips[:n_per_row]
    lower_strips = strips[n_per_row:]
    upper_currents_a = currents[:n_per_row]
    lower_currents_a = currents[n_per_row:]

    upper_x_m = np.asarray([strip.x_m for strip in upper_strips])
    lower_x_m = np.asarray([strip.x_m for strip in lower_strips])
    if not np.allclose(upper_x_m, lower_x_m, rtol=0.0, atol=1.0e-15):
        raise ValueError("Upper and lower strip rows use different x grids")

    order = np.argsort(upper_x_m)
    x_mm = upper_x_m[order] * 1.0e3
    plotted_currents_a = upper_currents_a[order]

    figure, axes = _make_2d_figure(keep_open=keep_open)
    axes.plot(
        x_mm,
        plotted_currents_a,
        marker="o",
        markersize=3.0,
        linewidth=1.2,
    )
    axes.axhline(0.0, linewidth=0.8)
    axes.set_title("Fitted strip currents — one plane")
    axes.set_xlabel("strip-centre x [mm]")
    axes.set_ylabel("current [A]")
    axes.set_xlim(-x_range_mm, x_range_mm)
    axes.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)

    symmetry_difference_a = upper_currents_a - lower_currents_a
    max_abs_mismatch_a = float(np.max(np.abs(symmetry_difference_a)))
    rms_mismatch_a = float(np.sqrt(np.mean(symmetry_difference_a ** 2)))

    if keep_open:
        return figure, max_abs_mismatch_a, rms_mismatch_a

    try:
        plt.close(figure)
    except Exception:
        pass
    return None, max_abs_mismatch_a, rms_mismatch_a


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
    x_range_mm: float | None,
    y_range_mm: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Crop, sort, and convert one kick component for 3D plotting."""
    x_plot_m, values_plot_urad = crop_horizontal_range(
        x_m, values_urad, x_range_mm
    )
    y_plot_m, values_plot_urad = crop_vertical_range(
        y_m, values_plot_urad, y_range_mm
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
    x_range_mm: float | None = None,
    y_range_mm: float | None = None,
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
        x_range_mm,
        y_range_mm,
    )
    theta_y_grid = _prepare_3d_component(
        x_m,
        y_m,
        theta_y_urad,
        x_range_mm,
        y_range_mm,
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

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Current-strip correction prototype for an EPU kick map"
    )
    parser.add_argument("--kick", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, default=Path("epu57_fixed"))
    parser.add_argument("--bRho", type=float, default=10.0)
    parser.add_argument(
        "--kick-scale",
        type=float,
        default=1.0,
        metavar="FACTOR",
        help=(
            "multiply both input kick components by FACTOR before fitting, "
            "correction, plotting, metrics, and export (default: 1.0)"
        ),
    )
    parser.add_argument(
        "--horizontal-weight",
        type=float,
        default=1.0,
        metavar="WEIGHT",
        help=(
            "relative non-negative weight for the horizontal y=0 fit "
            "residual; weights are normalized by their maximum "
            "(default: 1.0)"
        ),
    )
    parser.add_argument(
        "--vertical-weight",
        type=float,
        default=0.0,
        metavar="WEIGHT",
        help=(
            "relative non-negative weight for the vertical x=0 fit "
            "residual; weights are normalized by their maximum; 0 selects "
            "a horizontal-only objective (default: 0.0)"
        ),
    )
    parser.add_argument(
        "--target-mode",
        choices=("nonlinear", "full"),
        default="nonlinear",
        help=(
            "fit target: 'nonlinear' projects out constant and linear "
            "components from each target cut and response column before "
            "optimization (default); 'full' restores the R15 full-kick fit"
        ),
    )
    parser.add_argument("--gap", type=float, default=0.016)
    parser.add_argument(
        "--length",
        type=float,
        default=None,
        help="Strip length in m; default is the kick-map undulator length",
    )
    parser.add_argument(
        "--strips-per-plane",
        type=int,
        default=21,
        metavar="N",
        help=(
            "number of strips on each current sheet; centres are distributed "
            "uniformly from -x-range-mm to +x-range-mm (default: 21)"
        ),
    )
    parser.add_argument("--width", type=float, default=2.0e-3)
    parser.add_argument("--thickness", type=float, default=0.3e-3)
    parser.add_argument("--ridge", type=float, default=1.0e-5)
    parser.add_argument("--limit", type=float, default=6.0)
    parser.add_argument("--quadrature-order", type=int, default=4)
    parser.add_argument(
        "--max-fit-relative-rms",
        type=float,
        default=None,
        metavar="FRACTION",
        help=(
            "maximum accepted weighted RMS fit residual divided by the "
            "weighted RMS target; defaults to 0.15 for nonlinear targets "
            "and 0.05 for full targets"
        ),
    )
    parser.add_argument(
        "--min-horizontal-rms-reduction-factor",
        type=float,
        default=1.0,
        metavar="FACTOR",
        help=(
            "minimum accepted raw/corrected horizontal y=0 RMS ratio on "
            "the active target basis (nonlinear residual or full cut); "
            "default: 1.0"
        ),
    )
    parser.add_argument(
        "--current-bound-tolerance",
        type=float,
        default=1.0e-6,
        metavar="FRACTION",
        help=(
            "relative tolerance for reporting currents at the imposed bound "
            "(default: 1e-6)"
        ),
    )
    parser.add_argument(
        "--x-range-mm",
        type=float,
        default=20.0,
        metavar="MM",
        help=(
            "set the outermost strip-centre positions and limit displayed "
            "kick maps/cuts to x = ±MM (default: 20 mm); the exported kick "
            "map remains full-range"
        ),
    )
    parser.add_argument(
        "--y-range-mm",
        type=float,
        default=None,
        metavar="MM",
        help=(
            "limit displayed 3D kick maps and the vertical centreline cut "
            "to y = ±MM; the fit and exported kick map remain full-range"
        ),
    )
    parser.add_argument(
        "--color-scale",
        choices=("auto", "shared"),
        default="auto",
        help=(
            "3D colour/z scaling: 'auto' uses each displayed panel's data "
            "(default); 'shared' uses common symmetric raw/corrected scales"
        ),
    )
    parser.add_argument(
        "--elevation-deg",
        type=float,
        default=PLOT_ELEVATION_DEG,
        metavar="DEG",
        help=(
            "common elevation angle for all 3D panels "
            f"(default: {PLOT_ELEVATION_DEG:g} degrees)"
        ),
    )
    parser.add_argument(
        "--azimuth-deg",
        type=float,
        default=PLOT_AZIMUTH_DEG,
        metavar="DEG",
        help=(
            "common azimuth angle for all 3D panels "
            f"(default: {PLOT_AZIMUTH_DEG:g} degrees)"
        ),
    )
    parser.add_argument(
        "--projection",
        choices=("ortho", "persp"),
        default="ortho",
        help=(
            "common 3D projection for both panels "
            "(default: ortho, closest to Gnuplot)"
        ),
    )
    parser.add_argument(
        "--box-aspect-x",
        type=float,
        default=2.5,
        metavar="VALUE",
        help="common 3D box aspect in x (default: 2.5)",
    )
    parser.add_argument(
        "--box-aspect-y",
        type=float,
        default=1.8,
        metavar="VALUE",
        help="common 3D box aspect in y (default: 1.8)",
    )
    parser.add_argument(
        "--box-aspect-z",
        type=float,
        default=0.85,
        metavar="VALUE",
        help="common 3D box aspect in z (default: 0.85)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help=(
            "display the current plot and paired raw/corrected 3D figures "
            "after saving them"
        ),
    )
    arguments = parser.parse_args()

    if not np.isfinite(arguments.bRho) or arguments.bRho <= 0.0:
        raise ValueError("bRho must be finite and positive")
    if not np.isfinite(arguments.kick_scale) or arguments.kick_scale <= 0.0:
        raise ValueError("kick-scale must be finite and positive")
    for name, value in (
        ("horizontal-weight", arguments.horizontal_weight),
        ("vertical-weight", arguments.vertical_weight),
    ):
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
    if arguments.horizontal_weight == 0.0 and arguments.vertical_weight == 0.0:
        raise ValueError(
            "at least one of horizontal-weight or vertical-weight must be positive"
        )
    if not np.isfinite(arguments.x_range_mm) or arguments.x_range_mm <= 0.0:
        raise ValueError("x-range-mm must be finite and positive")
    if arguments.strips_per_plane < 1:
        raise ValueError("strips-per-plane must be positive")
    if arguments.y_range_mm is not None and (
        not np.isfinite(arguments.y_range_mm) or arguments.y_range_mm <= 0.0
    ):
        raise ValueError("y-range-mm must be finite and positive")
    if arguments.max_fit_relative_rms is not None and (
        not np.isfinite(arguments.max_fit_relative_rms)
        or arguments.max_fit_relative_rms < 0.0
    ):
        raise ValueError("max-fit-relative-rms must be finite and non-negative")
    if arguments.max_fit_relative_rms is None:
        effective_max_fit_relative_rms = (
            0.15 if arguments.target_mode == "nonlinear" else 0.05
        )
        max_fit_relative_rms_source = "target_mode_default"
    else:
        effective_max_fit_relative_rms = arguments.max_fit_relative_rms
        max_fit_relative_rms_source = "user"
    if not np.isfinite(arguments.min_horizontal_rms_reduction_factor) or (
        arguments.min_horizontal_rms_reduction_factor <= 0.0
    ):
        raise ValueError(
            "min-horizontal-rms-reduction-factor must be finite and positive"
        )
    if not np.isfinite(arguments.current_bound_tolerance) or not (
        0.0 <= arguments.current_bound_tolerance < 1.0
    ):
        raise ValueError("current-bound-tolerance must lie in [0, 1)")
    if not np.isfinite(arguments.elevation_deg):
        raise ValueError("elevation-deg must be finite")
    if not np.isfinite(arguments.azimuth_deg):
        raise ValueError("azimuth-deg must be finite")
    for name, value in (
        ("box-aspect-x", arguments.box_aspect_x),
        ("box-aspect-y", arguments.box_aspect_y),
        ("box-aspect-z", arguments.box_aspect_z),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")

    plot_elevation_deg = arguments.elevation_deg
    plot_azimuth_deg = arguments.azimuth_deg

    kick_map_unscaled = load_kick_map(arguments.kick)
    kick_map = KickMap(
        length_m=kick_map_unscaled.length_m,
        x_m=kick_map_unscaled.x_m,
        y_m=kick_map_unscaled.y_m,
        theta_x_urad=(
            arguments.kick_scale * kick_map_unscaled.theta_x_urad
        ),
        theta_y_urad=(
            arguments.kick_scale * kick_map_unscaled.theta_y_urad
        ),
        # Integrated B^2 is not a kick angle and is intentionally unscaled.
        b2_t2m=kick_map_unscaled.b2_t2m,
    )
    strip_length_m = (
        kick_map.length_m if arguments.length is None else arguments.length
    )

    strips, derived_pitch_m = make_two_strip_rows(
        centre_half_span_m=arguments.x_range_mm * 1.0e-3,
        n_per_row=arguments.strips_per_plane,
        width_m=arguments.width,
        thickness_m=arguments.thickness,
        gap_m=arguments.gap,
        length_m=strip_length_m,
    )
    n_strips_per_row = arguments.strips_per_plane
    upper_strip_centres_m = np.asarray(
        [strip.x_m for strip in strips[:n_strips_per_row]],
        dtype=float,
    )
    strip_outer_half_width_m = float(
        np.max(np.abs(upper_strip_centres_m)) + arguments.width / 2.0
    )

    raw_theta_x_y0 = zero_coordinate_cut(
        kick_map.y_m,
        kick_map.theta_x_urad,
        axis=0,
        coordinate_name="y",
    )
    raw_theta_y_x0 = zero_coordinate_cut(
        kick_map.x_m,
        kick_map.theta_y_urad,
        axis=1,
        coordinate_name="x",
    )

    # Fit one physical current vector to independently weighted horizontal
    # and vertical cancellation targets.  Dividing each block by sqrt(N)
    # makes the weights component-level coefficients rather than implicit
    # multipliers of the number of sampled points:
    #     objective = w_x*mean(r_x^2) + w_y*mean(r_y^2).
    horizontal_fit_points = np.column_stack(
        [kick_map.x_m, np.zeros_like(kick_map.x_m)]
    )
    vertical_fit_points = np.column_stack(
        [np.zeros_like(kick_map.y_m), kick_map.y_m]
    )
    response_iy_midplane = response_matrix(
        strips,
        horizontal_fit_points,
        component="Iy",
        quadrature_order=arguments.quadrature_order,
    )
    response_ix_centerline = response_matrix(
        strips,
        vertical_fit_points,
        component="Ix",
        quadrature_order=arguments.quadrature_order,
    )

    raw_theta_x_midplane_rad = raw_theta_x_y0.values * 1.0e-6
    raw_theta_y_centerline_rad = raw_theta_y_x0.values * 1.0e-6

    # theta_x(strip) = -Iy/(B rho), so cancellation requires
    # Iy_target = +(B rho) * theta_x(raw).
    target_iy_tm = arguments.bRho * raw_theta_x_midplane_rad
    # theta_y(strip) = +Ix/(B rho), so cancellation requires
    # Ix_target = -(B rho) * theta_y(raw).
    target_ix_tm = -arguments.bRho * raw_theta_y_centerline_rad

    raw_theta_x_affine = affine_decomposition(
        kick_map.x_m, raw_theta_x_y0.values
    )
    raw_theta_y_affine = affine_decomposition(
        kick_map.y_m, raw_theta_y_x0.values
    )
    target_iy_affine = affine_decomposition(kick_map.x_m, target_iy_tm)
    target_ix_affine = affine_decomposition(kick_map.y_m, target_ix_tm)
    response_iy_affine = affine_decomposition(
        kick_map.x_m, response_iy_midplane
    )
    response_ix_affine = affine_decomposition(
        kick_map.y_m, response_ix_centerline
    )

    if arguments.target_mode == "nonlinear":
        fit_target_iy_tm = target_iy_affine.residual
        fit_target_ix_tm = target_ix_affine.residual
        fit_response_iy = response_iy_affine.residual
        fit_response_ix = response_ix_affine.residual
    else:
        fit_target_iy_tm = target_iy_tm
        fit_target_ix_tm = target_ix_tm
        fit_response_iy = response_iy_midplane
        fit_response_ix = response_ix_centerline

    fit_weight_scale = max(
        arguments.horizontal_weight,
        arguments.vertical_weight,
    )
    normalized_horizontal_weight = (
        arguments.horizontal_weight / fit_weight_scale
    )
    normalized_vertical_weight = arguments.vertical_weight / fit_weight_scale
    total_normalized_fit_weight = (
        normalized_horizontal_weight + normalized_vertical_weight
    )

    weighted_response_blocks: list[np.ndarray] = []
    weighted_target_blocks: list[np.ndarray] = []
    if normalized_horizontal_weight > 0.0:
        horizontal_scale = np.sqrt(
            normalized_horizontal_weight / fit_target_iy_tm.size
        )
        weighted_response_blocks.append(horizontal_scale * fit_response_iy)
        weighted_target_blocks.append(horizontal_scale * fit_target_iy_tm)
    if normalized_vertical_weight > 0.0:
        vertical_scale = np.sqrt(
            normalized_vertical_weight / fit_target_ix_tm.size
        )
        weighted_response_blocks.append(vertical_scale * fit_response_ix)
        weighted_target_blocks.append(vertical_scale * fit_target_ix_tm)

    weighted_fit_response = np.vstack(weighted_response_blocks)
    weighted_fit_target = np.concatenate(weighted_target_blocks)

    current_solution = solve_currents(
        response=weighted_fit_response,
        target=weighted_fit_target,
        ridge=arguments.ridge,
        current_limit_a=arguments.limit,
    )
    currents_a = current_solution.currents_a

    x_grid, y_grid = np.meshgrid(kick_map.x_m, kick_map.y_m)
    all_points = np.column_stack([x_grid.ravel(), y_grid.ravel()])

    response_iy_full = response_matrix(
        strips,
        all_points,
        component="Iy",
        quadrature_order=arguments.quadrature_order,
    )
    response_ix_full = response_matrix(
        strips,
        all_points,
        component="Ix",
        quadrature_order=arguments.quadrature_order,
    )

    strip_iy_tm = (response_iy_full @ currents_a).reshape(x_grid.shape)
    strip_ix_tm = (response_ix_full @ currents_a).reshape(x_grid.shape)

    # Right-handed Lorentz-force signs:
    # theta_x(strip) = -Iy/(B rho)
    # theta_y(strip) = +Ix/(B rho)
    theta_x_corrected_urad = (
        kick_map.theta_x_urad
        - 1.0e6 * strip_iy_tm / arguments.bRho
    )
    theta_y_corrected_urad = (
        kick_map.theta_y_urad
        + 1.0e6 * strip_ix_tm / arguments.bRho
    )

    corrected_theta_x_y0 = zero_coordinate_cut(
        kick_map.y_m,
        theta_x_corrected_urad,
        axis=0,
        coordinate_name="y",
    )
    corrected_theta_y_x0 = zero_coordinate_cut(
        kick_map.x_m,
        theta_y_corrected_urad,
        axis=1,
        coordinate_name="x",
    )
    corrected_theta_x_affine = affine_decomposition(
        kick_map.x_m, corrected_theta_x_y0.values
    )
    corrected_theta_y_affine = affine_decomposition(
        kick_map.y_m, corrected_theta_y_x0.values
    )

    fit_prediction_iy_tm = fit_response_iy @ currents_a
    fit_residual_iy_tm = fit_prediction_iy_tm - fit_target_iy_tm
    fit_prediction_ix_tm = fit_response_ix @ currents_a
    fit_residual_ix_tm = fit_prediction_ix_tm - fit_target_ix_tm
    full_prediction_iy_tm = response_iy_midplane @ currents_a
    full_residual_iy_tm = full_prediction_iy_tm - target_iy_tm
    full_prediction_ix_tm = response_ix_centerline @ currents_a
    full_residual_ix_tm = full_prediction_ix_tm - target_ix_tm

    horizontal_fit_target_rms_tm = root_mean_square(fit_target_iy_tm)
    horizontal_fit_residual_rms_tm = root_mean_square(fit_residual_iy_tm)
    horizontal_fit_relative_rms = relative_rms_residual(
        fit_residual_iy_tm,
        fit_target_iy_tm,
    )
    vertical_fit_target_rms_tm = root_mean_square(fit_target_ix_tm)
    vertical_fit_residual_rms_tm = root_mean_square(fit_residual_ix_tm)
    vertical_fit_relative_rms = relative_rms_residual(
        fit_residual_ix_tm,
        fit_target_ix_tm,
    )
    horizontal_full_fit_relative_rms = relative_rms_residual(
        full_residual_iy_tm, target_iy_tm
    )
    vertical_full_fit_relative_rms = relative_rms_residual(
        full_residual_ix_tm, target_ix_tm
    )

    fit_target_rms_tm = float(
        np.sqrt(
            (
                normalized_horizontal_weight
                * horizontal_fit_target_rms_tm ** 2
                + normalized_vertical_weight
                * vertical_fit_target_rms_tm ** 2
            )
            / total_normalized_fit_weight
        )
    )
    fit_residual_rms_tm = float(
        np.sqrt(
            (
                normalized_horizontal_weight
                * horizontal_fit_residual_rms_tm ** 2
                + normalized_vertical_weight
                * vertical_fit_residual_rms_tm ** 2
            )
            / total_normalized_fit_weight
        )
    )
    fit_relative_rms = (
        0.0
        if fit_target_rms_tm == 0.0 and fit_residual_rms_tm == 0.0
        else (
            float("inf")
            if fit_target_rms_tm == 0.0
            else fit_residual_rms_tm / fit_target_rms_tm
        )
    )

    raw_theta_x_rms = root_mean_square(raw_theta_x_y0.values)
    corrected_theta_x_rms = root_mean_square(corrected_theta_x_y0.values)
    horizontal_rms_reduction_factor, horizontal_ratio_status = (
        safe_reduction_factor(raw_theta_x_rms, corrected_theta_x_rms)
    )
    raw_theta_x_nonlinear_rms = root_mean_square(
        raw_theta_x_affine.residual
    )
    corrected_theta_x_nonlinear_rms = root_mean_square(
        corrected_theta_x_affine.residual
    )
    (
        horizontal_nonlinear_rms_reduction_factor,
        horizontal_nonlinear_ratio_status,
    ) = safe_reduction_factor(
        raw_theta_x_nonlinear_rms, corrected_theta_x_nonlinear_rms
    )

    raw_theta_y_peak = float(np.max(np.abs(raw_theta_y_x0.values)))
    corrected_theta_y_peak = float(
        np.max(np.abs(corrected_theta_y_x0.values))
    )
    vertical_peak_reduction_factor, vertical_ratio_status = (
        safe_reduction_factor(raw_theta_y_peak, corrected_theta_y_peak)
    )
    raw_theta_y_nonlinear_peak = float(
        np.max(np.abs(raw_theta_y_affine.residual))
    )
    corrected_theta_y_nonlinear_peak = float(
        np.max(np.abs(corrected_theta_y_affine.residual))
    )
    (
        vertical_nonlinear_peak_reduction_factor,
        vertical_nonlinear_ratio_status,
    ) = safe_reduction_factor(
        raw_theta_y_nonlinear_peak, corrected_theta_y_nonlinear_peak
    )

    if arguments.target_mode == "nonlinear":
        horizontal_acceptance_reduction_factor = (
            horizontal_nonlinear_rms_reduction_factor
        )
        horizontal_acceptance_ratio_status = (
            horizontal_nonlinear_ratio_status
        )
        horizontal_acceptance_rms_basis = "nonlinear_residual"
    else:
        horizontal_acceptance_reduction_factor = (
            horizontal_rms_reduction_factor
        )
        horizontal_acceptance_ratio_status = horizontal_ratio_status
        horizontal_acceptance_rms_basis = "full_cut"

    bound_margin_a = arguments.limit * arguments.current_bound_tolerance
    bound_threshold_a = arguments.limit - bound_margin_a
    current_bound_active = np.abs(currents_a) >= bound_threshold_a
    current_bound_active_count = int(np.count_nonzero(current_bound_active))
    current_bound_active_fraction = float(
        current_bound_active_count / currents_a.size
    )

    acceptance_failures: list[str] = []
    if fit_relative_rms > effective_max_fit_relative_rms:
        acceptance_failures.append(
            "fit relative RMS "
            f"{fit_relative_rms:.6g} exceeds "
            f"{effective_max_fit_relative_rms:.6g}"
        )
    if arguments.horizontal_weight > 0.0 and (
        horizontal_acceptance_reduction_factor
        < arguments.min_horizontal_rms_reduction_factor
    ):
        acceptance_failures.append(
            "horizontal RMS reduction factor "
            f"{horizontal_acceptance_reduction_factor:.6g} is below "
            f"{arguments.min_horizontal_rms_reduction_factor:.6g}"
        )

    if acceptance_failures:
        acceptance_status = "FAILED"
    elif current_bound_active_count:
        acceptance_status = "ACCEPTED_LIMIT_SATURATED"
    else:
        acceptance_status = "ACCEPTED"

    prefix = arguments.out_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)

    output_suffixes = [
        "_strip_currents.csv",
        "_strip_currents_one_plane.png",
        "_corrected_kickmap.dat",
        "_metrics.txt",
        "_plot_manifest.txt",
        "_theta_x_midplane_cut.png",
        "_theta_y_centerline_cut.png",
        "_raw_kickmaps_3d.png",
        "_corrected_kickmaps_3d.png",
    ]
    final_paths = [
        prefix.with_name(prefix.name + suffix) for suffix in output_suffixes
    ]

    metrics = (
        f"release_id={RELEASE_ID}\n"
        f"input_file={arguments.kick}\n"
        f"grid_shape={kick_map.theta_x_urad.shape}\n"
        f"length_m={kick_map.length_m:.12g}\n"
        f"kick_scale={arguments.kick_scale:.12g}\n"
        f"horizontal_weight={arguments.horizontal_weight:.12g}\n"
        f"vertical_weight={arguments.vertical_weight:.12g}\n"
        f"target_mode={arguments.target_mode}\n"
        f"affine_projection_basis=constant_plus_linear_coordinate\n"
        f"normalized_horizontal_weight="
        f"{normalized_horizontal_weight:.12g}\n"
        f"normalized_vertical_weight="
        f"{normalized_vertical_weight:.12g}\n"
        f"fit_weighting=max_normalized_component_weight_times_mean_squared_residual\n"
        f"color_scale={arguments.color_scale}\n"
        f"plot_elevation_deg={plot_elevation_deg:.12g}\n"
        f"plot_azimuth_deg={plot_azimuth_deg:.12g}\n"
        f"projection={arguments.projection}\n"
        f"box_aspect_x={arguments.box_aspect_x:.12g}\n"
        f"box_aspect_y={arguments.box_aspect_y:.12g}\n"
        f"box_aspect_z={arguments.box_aspect_z:.12g}\n"
        f"theta_x_y0_cut_method={raw_theta_x_y0.method}\n"
        f"theta_x_y0_bracket_lower_m="
        f"{raw_theta_x_y0.lower_coordinate_m:.12g}\n"
        f"theta_x_y0_bracket_upper_m="
        f"{raw_theta_x_y0.upper_coordinate_m:.12g}\n"
        f"theta_y_x0_cut_method={raw_theta_y_x0.method}\n"
        f"theta_y_x0_bracket_lower_m="
        f"{raw_theta_y_x0.lower_coordinate_m:.12g}\n"
        f"theta_y_x0_bracket_upper_m="
        f"{raw_theta_y_x0.upper_coordinate_m:.12g}\n"
        f"raw_theta_x_y0_min_urad="
        f"{np.min(raw_theta_x_y0.values):.12g}\n"
        f"raw_theta_x_y0_max_urad="
        f"{np.max(raw_theta_x_y0.values):.12g}\n"
        f"raw_theta_y_x0_min_urad="
        f"{np.min(raw_theta_y_x0.values):.12g}\n"
        f"raw_theta_y_x0_max_urad="
        f"{np.max(raw_theta_y_x0.values):.12g}\n"
        f"raw_theta_x_affine_offset_urad="
        f"{float(raw_theta_x_affine.offset):.12g}\n"
        f"raw_theta_x_affine_gradient_urad_per_mm="
        f"{float(raw_theta_x_affine.slope_per_m) * 1.0e-3:.12g}\n"
        f"corrected_theta_x_affine_offset_urad="
        f"{float(corrected_theta_x_affine.offset):.12g}\n"
        f"corrected_theta_x_affine_gradient_urad_per_mm="
        f"{float(corrected_theta_x_affine.slope_per_m) * 1.0e-3:.12g}\n"
        f"raw_theta_y_affine_offset_urad="
        f"{float(raw_theta_y_affine.offset):.12g}\n"
        f"raw_theta_y_affine_gradient_urad_per_mm="
        f"{float(raw_theta_y_affine.slope_per_m) * 1.0e-3:.12g}\n"
        f"corrected_theta_y_affine_offset_urad="
        f"{float(corrected_theta_y_affine.offset):.12g}\n"
        f"corrected_theta_y_affine_gradient_urad_per_mm="
        f"{float(corrected_theta_y_affine.slope_per_m) * 1.0e-3:.12g}\n"
        f"raw_theta_x_y0_rms_urad={raw_theta_x_rms:.12g}\n"
        f"corrected_theta_x_y0_rms_urad="
        f"{corrected_theta_x_rms:.12g}\n"
        f"horizontal_rms_reduction_factor="
        f"{horizontal_rms_reduction_factor:.12g}\n"
        f"horizontal_rms_reduction_status={horizontal_ratio_status}\n"
        f"raw_theta_x_y0_nonlinear_rms_urad="
        f"{raw_theta_x_nonlinear_rms:.12g}\n"
        f"corrected_theta_x_y0_nonlinear_rms_urad="
        f"{corrected_theta_x_nonlinear_rms:.12g}\n"
        f"horizontal_nonlinear_rms_reduction_factor="
        f"{horizontal_nonlinear_rms_reduction_factor:.12g}\n"
        f"horizontal_nonlinear_rms_reduction_status="
        f"{horizontal_nonlinear_ratio_status}\n"
        f"horizontal_acceptance_rms_basis="
        f"{horizontal_acceptance_rms_basis}\n"
        f"horizontal_acceptance_reduction_factor="
        f"{horizontal_acceptance_reduction_factor:.12g}\n"
        f"horizontal_acceptance_reduction_status="
        f"{horizontal_acceptance_ratio_status}\n"
        f"raw_theta_y_x0_peak_urad={raw_theta_y_peak:.12g}\n"
        f"corrected_theta_y_x0_peak_urad="
        f"{corrected_theta_y_peak:.12g}\n"
        f"vertical_peak_reduction_factor="
        f"{vertical_peak_reduction_factor:.12g}\n"
        f"vertical_peak_reduction_status={vertical_ratio_status}\n"
        f"raw_theta_y_x0_nonlinear_peak_urad="
        f"{raw_theta_y_nonlinear_peak:.12g}\n"
        f"corrected_theta_y_x0_nonlinear_peak_urad="
        f"{corrected_theta_y_nonlinear_peak:.12g}\n"
        f"vertical_nonlinear_peak_reduction_factor="
        f"{vertical_nonlinear_peak_reduction_factor:.12g}\n"
        f"vertical_nonlinear_peak_reduction_status="
        f"{vertical_nonlinear_ratio_status}\n"
        f"solver_status={current_solution.status}\n"
        f"solver_message={current_solution.message}\n"
        f"solver_cost={current_solution.cost:.12g}\n"
        f"solver_optimality={current_solution.optimality:.12g}\n"
        f"solver_iterations={current_solution.iterations}\n"
        f"horizontal_fit_target_rms_Tm="
        f"{horizontal_fit_target_rms_tm:.12g}\n"
        f"horizontal_fit_residual_rms_Tm="
        f"{horizontal_fit_residual_rms_tm:.12g}\n"
        f"horizontal_fit_relative_rms="
        f"{horizontal_fit_relative_rms:.12g}\n"
        f"vertical_fit_target_rms_Tm="
        f"{vertical_fit_target_rms_tm:.12g}\n"
        f"vertical_fit_residual_rms_Tm="
        f"{vertical_fit_residual_rms_tm:.12g}\n"
        f"vertical_fit_relative_rms="
        f"{vertical_fit_relative_rms:.12g}\n"
        f"horizontal_full_fit_relative_rms="
        f"{horizontal_full_fit_relative_rms:.12g}\n"
        f"vertical_full_fit_relative_rms="
        f"{vertical_full_fit_relative_rms:.12g}\n"
        f"fit_target_rms_Tm={fit_target_rms_tm:.12g}\n"
        f"fit_residual_rms_Tm={fit_residual_rms_tm:.12g}\n"
        f"fit_relative_rms={fit_relative_rms:.12g}\n"
        f"max_fit_relative_rms="
        f"{effective_max_fit_relative_rms:.12g}\n"
        f"max_fit_relative_rms_source={max_fit_relative_rms_source}\n"
        f"min_horizontal_rms_reduction_factor="
        f"{arguments.min_horizontal_rms_reduction_factor:.12g}\n"
        f"current_bound_tolerance="
        f"{arguments.current_bound_tolerance:.12g}\n"
        f"current_bound_active_count={current_bound_active_count}\n"
        f"current_bound_active_fraction="
        f"{current_bound_active_fraction:.12g}\n"
        f"acceptance_status={acceptance_status}\n"
        f"x_range_mm={arguments.x_range_mm:.12g}\n"
        f"n_strips_per_row={n_strips_per_row}\n"
        f"strip_pitch_mm={derived_pitch_m * 1.0e3:.12g}\n"
        f"strip_width_mm={arguments.width * 1.0e3:.12g}\n"
        f"strip_outer_half_width_mm="
        f"{strip_outer_half_width_m * 1.0e3:.12g}\n"
        f"max_abs_current_A={np.max(np.abs(currents_a)):.12g}\n"
    )

    if acceptance_failures:
        print(metrics, end="", file=sys.stderr)
        raise RuntimeError(
            "Correction acceptance failed: " + "; ".join(acceptance_failures)
        )

    with tempfile.TemporaryDirectory(
        dir=prefix.parent,
        prefix=f".{prefix.name}_staging_",
    ) as staging_directory:
        staging_prefix = Path(staging_directory) / prefix.name
        staged_paths = [
            staging_prefix.with_name(staging_prefix.name + suffix)
            for suffix in output_suffixes
        ]
        (
            currents_path,
            currents_plot_path,
            corrected_map_path,
            metrics_path,
            plot_manifest_path,
            theta_x_cut_path,
            theta_y_cut_path,
            raw_3d_path,
            corrected_3d_path,
        ) = staged_paths

        np.savetxt(
            currents_path,
            np.column_stack([np.arange(currents_a.size), currents_a]),
            delimiter=",",
            header="index,current_A",
            comments="",
            fmt=["%d", "%.10e"],
        )

        (
            currents_figure,
            current_symmetry_max_abs_mismatch_a,
            current_symmetry_rms_mismatch_a,
        ) = plot_strip_currents_one_plane(
            strips,
            currents_a,
            n_strips_per_row,
            currents_plot_path,
            arguments.x_range_mm,
            keep_open=arguments.show,
        )

        write_corrected_kick_map(
            corrected_map_path,
            kick_map,
            theta_x_corrected_urad,
            theta_y_corrected_urad,
        )

        theta_x_limit = float(
            max(
                np.max(np.abs(kick_map.theta_x_urad)),
                np.max(np.abs(theta_x_corrected_urad)),
            )
        )
        theta_y_limit = float(
            max(
                np.max(np.abs(kick_map.theta_y_urad)),
                np.max(np.abs(theta_y_corrected_urad)),
            )
        )

        if arguments.target_mode == "nonlinear":
            theta_x_plot_raw = raw_theta_x_affine.residual
            theta_x_plot_corrected = corrected_theta_x_affine.residual
            theta_x_plot_title = (
                "Nonlinear horizontal kick on y = 0 "
                "(affine component removed)"
            )
        else:
            theta_x_plot_raw = raw_theta_x_y0.values
            theta_x_plot_corrected = corrected_theta_x_y0.values
            theta_x_plot_title = (
                "Horizontal kick on the horizontal mid-plane (y = 0)"
            )
        x_cut_m, raw_theta_x_cut_urad = crop_horizontal_range(
            kick_map.x_m, theta_x_plot_raw, arguments.x_range_mm
        )
        _, corrected_theta_x_cut_urad = crop_horizontal_range(
            kick_map.x_m, theta_x_plot_corrected, arguments.x_range_mm
        )
        plot_cut(
            x_cut_m * 1.0e3,
            raw_theta_x_cut_urad,
            corrected_theta_x_cut_urad,
            theta_x_plot_title,
            "x [mm]",
            "θx [µrad]",
            theta_x_cut_path,
        )

        if arguments.target_mode == "nonlinear":
            theta_y_plot_raw = raw_theta_y_affine.residual
            theta_y_plot_corrected = corrected_theta_y_affine.residual
            theta_y_plot_title = (
                "Nonlinear vertical kick on x = 0 "
                "(affine component removed)"
            )
        else:
            theta_y_plot_raw = raw_theta_y_x0.values
            theta_y_plot_corrected = corrected_theta_y_x0.values
            theta_y_plot_title = (
                "Vertical kick on the vertical centreline (x = 0)"
            )
        y_cut_m, raw_theta_y_cut_urad = crop_vertical_range(
            kick_map.y_m,
            theta_y_plot_raw[:, np.newaxis],
            arguments.y_range_mm,
        )
        _, corrected_theta_y_cut_urad = crop_vertical_range(
            kick_map.y_m,
            theta_y_plot_corrected[:, np.newaxis],
            arguments.y_range_mm,
        )
        plot_cut(
            y_cut_m * 1.0e3,
            raw_theta_y_cut_urad[:, 0],
            corrected_theta_y_cut_urad[:, 0],
            theta_y_plot_title,
            "y [mm]",
            "θy [µrad]",
            theta_y_cut_path,
        )

        raw_kickmaps_3d_figure = plot_kickmaps_3d_pair(
            kick_map.x_m,
            kick_map.y_m,
            kick_map.theta_x_urad,
            kick_map.theta_y_urad,
            raw_3d_path,
            theta_x_limit,
            theta_y_limit,
            "EPU57 kick maps before current-strip correction",
            x_range_mm=arguments.x_range_mm,
            y_range_mm=arguments.y_range_mm,
            color_scale=arguments.color_scale,
            elevation_deg=plot_elevation_deg,
            azimuth_deg=plot_azimuth_deg,
            projection=arguments.projection,
            box_aspect_x=arguments.box_aspect_x,
            box_aspect_y=arguments.box_aspect_y,
            box_aspect_z=arguments.box_aspect_z,
            keep_open=arguments.show,
        )
        corrected_kickmaps_3d_figure = plot_kickmaps_3d_pair(
            kick_map.x_m,
            kick_map.y_m,
            theta_x_corrected_urad,
            theta_y_corrected_urad,
            corrected_3d_path,
            theta_x_limit,
            theta_y_limit,
            "EPU57 kick maps after current-strip correction",
            x_range_mm=arguments.x_range_mm,
            y_range_mm=arguments.y_range_mm,
            color_scale=arguments.color_scale,
            elevation_deg=plot_elevation_deg,
            azimuth_deg=plot_azimuth_deg,
            projection=arguments.projection,
            box_aspect_x=arguments.box_aspect_x,
            box_aspect_y=arguments.box_aspect_y,
            box_aspect_z=arguments.box_aspect_z,
            keep_open=arguments.show,
        )

        metrics += (
            f"current_symmetry_max_abs_mismatch_A="
            f"{current_symmetry_max_abs_mismatch_a:.12g}\n"
            f"current_symmetry_rms_mismatch_A="
            f"{current_symmetry_rms_mismatch_a:.12g}\n"
        )
        metrics_path.write_text(metrics, encoding="utf-8", newline="\n")

        plot_manifest = (
            f"release_id={RELEASE_ID}\n"
            f"acceptance_status={acceptance_status}\n"
            f"horizontal_weight={arguments.horizontal_weight:.12g}\n"
            f"vertical_weight={arguments.vertical_weight:.12g}\n"
            f"target_mode={arguments.target_mode}\n"
            f"normalized_horizontal_weight="
            f"{normalized_horizontal_weight:.12g}\n"
            f"normalized_vertical_weight="
            f"{normalized_vertical_weight:.12g}\n"
            f"strip_currents_one_plane_file="
            f"{prefix.name}_strip_currents_one_plane.png\n"
            f"strip_centre_half_span_mm={arguments.x_range_mm:.12g}\n"
            f"strip_pitch_mm={derived_pitch_m * 1.0e3:.12g}\n"
            f"strip_outer_half_width_actual_mm="
            f"{strip_outer_half_width_m * 1.0e3:.12g}\n"
            f"n_strips_per_row={n_strips_per_row}\n"
            f"raw_kickmaps_file={prefix.name}_raw_kickmaps_3d.png\n"
            f"corrected_kickmaps_file="
            f"{prefix.name}_corrected_kickmaps_3d.png\n"
            f"plot_elevation_deg={plot_elevation_deg:.12g}\n"
            f"plot_azimuth_deg={plot_azimuth_deg:.12g}\n"
        )
        plot_manifest_path.write_text(
            plot_manifest,
            encoding="utf-8",
            newline="\n",
        )

        if arguments.show:
            plt.show()
            for figure in (
                currents_figure,
                raw_kickmaps_3d_figure,
                corrected_kickmaps_3d_figure,
            ):
                if figure is not None:
                    plt.close(figure)

        publish_staged_outputs(staged_paths, final_paths)

    print(metrics, end="")
    print(plot_manifest, end="")


if __name__ == "__main__":
    main()
