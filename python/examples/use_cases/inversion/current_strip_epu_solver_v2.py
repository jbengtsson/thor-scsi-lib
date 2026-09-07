#!/usr/bin/env python3
"""
current_strip_epu_57_solver_13.py

Self-contained current-strip correction prototype for the supplied EPU57
RADIA kick map.

Key corrections relative to current_strip_epu_57_solver_2.py
-------------------------------------------------------------
1. The vertical kick is no longer sampled on y = 0, where it vanishes by
   symmetry.  Vertical diagnostics use the full 2D map and the x = 0 cut.
2. One physical strip-current vector is fitted to the horizontal mid-plane
   target and then applied to BOTH kick components.
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
import re
import sys
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

RELEASE_ID = "EPU57-SOLVER-13-21-STRIPS-PER-SHEET"
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
) -> np.ndarray:
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

    return result.x


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

    if arguments.bRho <= 0.0:
        raise ValueError("bRho must be positive")
    if not np.isfinite(arguments.kick_scale) or arguments.kick_scale <= 0.0:
        raise ValueError("kick-scale must be finite and positive")
    if not np.isfinite(arguments.x_range_mm) or arguments.x_range_mm <= 0.0:
        raise ValueError("x-range-mm must be finite and positive")
    if arguments.strips_per_plane < 1:
        raise ValueError("strips-per-plane must be positive")
    if arguments.y_range_mm is not None and arguments.y_range_mm <= 0.0:
        raise ValueError("y-range-mm must be positive")
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
        kick_map.length_m
        if arguments.length is None
        else arguments.length
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

    middle_y_index = int(np.argmin(np.abs(kick_map.y_m)))
    middle_x_index = int(np.argmin(np.abs(kick_map.x_m)))

    # The 2012 procedure fits one current distribution to the equivalent
    # first field integral on the horizontal mid-plane.
    fit_points = np.column_stack(
        [kick_map.x_m, np.zeros_like(kick_map.x_m)]
    )
    response_iy_midplane = response_matrix(
        strips,
        fit_points,
        component="Iy",
        quadrature_order=arguments.quadrature_order,
    )

    raw_theta_x_midplane_rad = (
        kick_map.theta_x_urad[middle_y_index, :] * 1.0e-6
    )

    # theta_x(strip) = -Iy/(B rho), so cancellation requires
    # Iy_target = (B rho) * theta_x(raw).
    target_iy_tm = arguments.bRho * raw_theta_x_midplane_rad

    currents_a = solve_currents(
        response=response_iy_midplane,
        target=target_iy_tm,
        ridge=arguments.ridge,
        current_limit_a=arguments.limit,
    )

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

    prefix = arguments.out_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)

    currents_path = prefix.with_name(prefix.name + "_strip_currents.csv")
    currents_plot_path = prefix.with_name(
        prefix.name + "_strip_currents_one_plane.png"
    )
    corrected_map_path = prefix.with_name(
        prefix.name + "_corrected_kickmap.dat"
    )
    metrics_path = prefix.with_name(prefix.name + "_metrics.txt")
    plot_manifest_path = prefix.with_name(
        prefix.name + "_plot_manifest.txt"
    )

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

    x_cut_m, raw_theta_x_cut_urad = crop_horizontal_range(
        kick_map.x_m,
        kick_map.theta_x_urad[middle_y_index, :],
        arguments.x_range_mm,
    )
    _, corrected_theta_x_cut_urad = crop_horizontal_range(
        kick_map.x_m,
        theta_x_corrected_urad[middle_y_index, :],
        arguments.x_range_mm,
    )
    plot_cut(
        x_cut_m * 1.0e3,
        raw_theta_x_cut_urad,
        corrected_theta_x_cut_urad,
        "Horizontal kick on the horizontal mid-plane (y = 0)",
        "x [mm]",
        "θx [µrad]",
        prefix.with_name(prefix.name + "_theta_x_midplane_cut.png"),
    )
    y_cut_m, raw_theta_y_cut_urad = crop_vertical_range(
        kick_map.y_m,
        kick_map.theta_y_urad[:, middle_x_index][:, np.newaxis],
        arguments.y_range_mm,
    )
    _, corrected_theta_y_cut_urad = crop_vertical_range(
        kick_map.y_m,
        theta_y_corrected_urad[:, middle_x_index][:, np.newaxis],
        arguments.y_range_mm,
    )
    plot_cut(
        y_cut_m * 1.0e3,
        raw_theta_y_cut_urad[:, 0],
        corrected_theta_y_cut_urad[:, 0],
        "Vertical kick on the vertical centreline (x = 0)",
        "y [mm]",
        "θy [µrad]",
        prefix.with_name(prefix.name + "_theta_y_centerline_cut.png"),
    )

    raw_kickmaps_3d_figure = plot_kickmaps_3d_pair(
        kick_map.x_m,
        kick_map.y_m,
        kick_map.theta_x_urad,
        kick_map.theta_y_urad,
        prefix.with_name(prefix.name + "_raw_kickmaps_3d.png"),
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
        prefix.with_name(prefix.name + "_corrected_kickmaps_3d.png"),
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

    if arguments.show:
        plt.show()
        for figure in (
            currents_figure,
            raw_kickmaps_3d_figure,
            corrected_kickmaps_3d_figure,
        ):
            if figure is not None:
                plt.close(figure)

    raw_theta_x_rms = float(
        np.sqrt(
            np.mean(
                kick_map.theta_x_urad[middle_y_index, :] ** 2
            )
        )
    )
    corrected_theta_x_rms = float(
        np.sqrt(
            np.mean(
                theta_x_corrected_urad[middle_y_index, :] ** 2
            )
        )
    )
    raw_theta_y_peak = float(
        np.max(np.abs(kick_map.theta_y_urad[:, middle_x_index]))
    )
    corrected_theta_y_peak = float(
        np.max(np.abs(theta_y_corrected_urad[:, middle_x_index]))
    )

    metrics = (
        f"release_id={RELEASE_ID}\n"
        f"input_file={arguments.kick}\n"
        f"grid_shape={kick_map.theta_x_urad.shape}\n"
        f"length_m={kick_map.length_m:.12g}\n"
        f"kick_scale={arguments.kick_scale:.12g}\n"
        f"color_scale={arguments.color_scale}\n"
        f"plot_elevation_deg={plot_elevation_deg:.12g}\n"
        f"plot_azimuth_deg={plot_azimuth_deg:.12g}\n"
        f"projection={arguments.projection}\n"
        f"box_aspect_x={arguments.box_aspect_x:.12g}\n"
        f"box_aspect_y={arguments.box_aspect_y:.12g}\n"
        f"box_aspect_z={arguments.box_aspect_z:.12g}\n"
        f"raw_theta_x_y0_min_urad="
        f"{np.min(kick_map.theta_x_urad[middle_y_index, :]):.12g}\n"
        f"raw_theta_x_y0_max_urad="
        f"{np.max(kick_map.theta_x_urad[middle_y_index, :]):.12g}\n"
        f"raw_theta_y_x0_min_urad="
        f"{np.min(kick_map.theta_y_urad[:, middle_x_index]):.12g}\n"
        f"raw_theta_y_x0_max_urad="
        f"{np.max(kick_map.theta_y_urad[:, middle_x_index]):.12g}\n"
        f"raw_theta_x_y0_rms_urad={raw_theta_x_rms:.12g}\n"
        f"corrected_theta_x_y0_rms_urad="
        f"{corrected_theta_x_rms:.12g}\n"
        f"raw_theta_y_x0_peak_urad={raw_theta_y_peak:.12g}\n"
        f"corrected_theta_y_x0_peak_urad="
        f"{corrected_theta_y_peak:.12g}\n"
        f"vertical_peak_reduction_factor="
        f"{raw_theta_y_peak / corrected_theta_y_peak:.12g}\n"
        f"x_range_mm={arguments.x_range_mm:.12g}\n"
        f"n_strips_per_row={n_strips_per_row}\n"
        f"strip_pitch_mm={derived_pitch_m * 1.0e3:.12g}\n"
        f"strip_width_mm={arguments.width * 1.0e3:.12g}\n"
        f"strip_outer_half_width_mm="
        f"{strip_outer_half_width_m * 1.0e3:.12g}\n"
        f"max_abs_current_A={np.max(np.abs(currents_a)):.12g}\n"
        f"current_symmetry_max_abs_mismatch_A="
        f"{current_symmetry_max_abs_mismatch_a:.12g}\n"
        f"current_symmetry_rms_mismatch_A="
        f"{current_symmetry_rms_mismatch_a:.12g}\n"
    )
    metrics_path.write_text(metrics, encoding="utf-8", newline="\n")

    plot_manifest = (
        f"release_id={RELEASE_ID}\n"
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

    print(metrics, end="")
    print(plot_manifest, end="")


if __name__ == "__main__":
    main()
