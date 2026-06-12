#!/usr/bin/env python3
"""
current_strip_epu_57_solver_3.py

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
5. Raw and corrected horizontal/vertical kick maps are saved as compact,
   aligned paired 3D surface plots; --show displays both figures.
6. Filled-contour 2D kick-map PNGs are not generated. The two 1D diagnostic
   cuts are retained.
7. --x-range-mm and --y-range-mm optionally limit the displayed ranges
   to ±X mm and ±Y mm without changing the fit, currents, or exported
   full kick map.
8. The paired 3D view follows the supplied gnuplot setting:
       set view 50, 348, 1, 1
9. The script has no dependency on the missing
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
import numpy as np
from scipy.optimize import lsq_linear

MU0 = 4.0e-7 * np.pi


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
    n_per_row: int,
    pitch_m: float,
    width_m: float,
    thickness_m: float,
    gap_m: float,
    length_m: float,
) -> list[Strip]:
    if n_per_row < 1:
        raise ValueError("n_per_row must be positive")
    if min(pitch_m, width_m, thickness_m, gap_m, length_m) <= 0.0:
        raise ValueError("All geometric dimensions must be positive")

    x_centres = (
        np.arange(n_per_row, dtype=float) - (n_per_row - 1) / 2.0
    ) * pitch_m

    # gap_m is the clear inner gap; the strip centres lie half a thickness
    # outside its two boundaries.
    y_centre = gap_m / 2.0 + thickness_m / 2.0

    return [
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

    figure = plt.figure(figsize=(8, 5))
    axes = figure.add_subplot(111)
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
    plt.close(figure)


def plot_cut(
    coordinate_mm: np.ndarray,
    raw_urad: np.ndarray,
    corrected_urad: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
    output_path: Path,
) -> None:
    figure = plt.figure(figsize=(8, 5))
    axes = figure.add_subplot(111)
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



def plot_kickmaps_3d(
    x_m: np.ndarray,
    y_m: np.ndarray,
    theta_x_urad: np.ndarray,
    theta_y_urad: np.ndarray,
    output_path: Path,
    theta_x_limit_urad: float,
    theta_y_limit_urad: float,
    x_range_mm: float | None = None,
    y_range_mm: float | None = None,
    keep_open: bool = False,
):
    """
    Save aligned paired 3D surfaces for theta_x and theta_y.

    The two panels use identical geometry and fixed axes positions. Raw and
    corrected figures share the same component-specific z and colour limits,
    allowing direct before/after comparison.
    """
    x_plot_m, theta_x_plot_urad = crop_horizontal_range(
        x_m, theta_x_urad, x_range_mm
    )
    _, theta_y_plot_urad = crop_horizontal_range(
        x_m, theta_y_urad, x_range_mm
    )
    y_plot_m, theta_x_plot_urad = crop_vertical_range(
        y_m, theta_x_plot_urad, y_range_mm
    )
    _, theta_y_plot_urad = crop_vertical_range(
        y_m, theta_y_plot_urad, y_range_mm
    )

    x_grid_mm, y_grid_mm = np.meshgrid(
        x_plot_m * 1.0e3,
        y_plot_m * 1.0e3,
    )

    # Explicit positions keep the two 3D axes and colour bars aligned,
    # independent of tick-label widths or the two different kick scales.
    figure = plt.figure(figsize=(10.8, 4.6))
    axes_positions = (
        (0.035, 0.12, 0.405, 0.80),
        (0.535, 0.12, 0.405, 0.80),
    )
    colorbar_positions = (
        (0.450, 0.31, 0.014, 0.38),
        (0.950, 0.31, 0.014, 0.38),
    )

    panels = (
        (
            theta_x_plot_urad * 1.0e-3,
            r"$\theta_x$ [mrad]",
            theta_x_limit_urad * 1.0e-3,
        ),
        (
            theta_y_plot_urad * 1.0e-3,
            r"$\theta_y$ [mrad]",
            theta_y_limit_urad * 1.0e-3,
        ),
    )

    for (
        values_mrad,
        panel_title,
        limit_mrad,
        axes_position,
        colorbar_position,
    ) in zip(
        (panel[0] for panel in panels),
        (panel[1] for panel in panels),
        (panel[2] for panel in panels),
        axes_positions,
        colorbar_positions,
    ):
        axes = figure.add_axes(axes_position, projection="3d")

        normalizer = matplotlib.colors.Normalize(
            vmin=-limit_mrad,
            vmax=limit_mrad,
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

        axes.set_title(panel_title, pad=1.0, fontsize=12)
        axes.set_xlabel("x [mm]", labelpad=2.0)
        axes.set_ylabel("y [mm]", labelpad=2.0)
        axes.set_zlabel("")

        axes.set_xlim(
            float(np.min(x_grid_mm)),
            float(np.max(x_grid_mm)),
        )
        axes.set_ylim(
            float(np.min(y_grid_mm)),
            float(np.max(y_grid_mm)),
        )
        axes.set_zlim(-limit_mrad, limit_mrad)

        # This compact orthographic view follows the supplied reference:
        # long x span, shallow y span, and matching panel baselines.
        axes.set_proj_type("ortho")
        # Gnuplot reference: set view 50, 348, 1, 1
        axes.view_init(elev=50.0, azim=348.0)
        axes.set_box_aspect((2.50, 0.80, 0.85))
        axes.tick_params(axis="both", labelsize=8, pad=0)

        for axis in (axes.xaxis, axes.yaxis, axes.zaxis):
            axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
            axis.pane.set_edgecolor((0.55, 0.55, 0.55, 0.55))

        color_axes = figure.add_axes(colorbar_position)
        colorbar = figure.colorbar(surface, cax=color_axes)
        colorbar.ax.tick_params(labelsize=8, pad=2)

    figure.savefig(output_path, dpi=180, bbox_inches="tight")

    if keep_open:
        return figure

    plt.close(figure)
    return None

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Current-strip correction prototype for an EPU kick map"
    )
    parser.add_argument("--kick", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, default=Path("epu57_fixed"))
    parser.add_argument("--bRho", type=float, default=10.0)
    parser.add_argument("--gap", type=float, default=0.016)
    parser.add_argument(
        "--length",
        type=float,
        default=None,
        help="Strip length in m; default is the kick-map undulator length",
    )
    parser.add_argument("--nstrips", type=int, default=48)
    parser.add_argument("--pitch", type=float, default=3.0e-3)
    parser.add_argument("--width", type=float, default=2.0e-3)
    parser.add_argument("--thickness", type=float, default=0.3e-3)
    parser.add_argument("--ridge", type=float, default=1.0e-5)
    parser.add_argument("--limit", type=float, default=6.0)
    parser.add_argument("--quadrature-order", type=int, default=4)
    parser.add_argument(
        "--x-range-mm",
        type=float,
        default=None,
        metavar="MM",
        help=(
            "limit displayed 3D kick maps and the horizontal cut to "
            "x = ±MM; the fit and exported kick map remain full-range"
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
        "--show",
        action="store_true",
        help=(
            "display the raw and corrected paired 3D kick-map figures "
            "after saving them"
        ),
    )
    arguments = parser.parse_args()

    if arguments.bRho <= 0.0:
        raise ValueError("bRho must be positive")
    if arguments.x_range_mm is not None and arguments.x_range_mm <= 0.0:
        raise ValueError("x-range-mm must be positive")
    if arguments.y_range_mm is not None and arguments.y_range_mm <= 0.0:
        raise ValueError("y-range-mm must be positive")

    kick_map = load_kick_map(arguments.kick)
    strip_length_m = (
        kick_map.length_m
        if arguments.length is None
        else arguments.length
    )

    strips = make_two_strip_rows(
        n_per_row=arguments.nstrips,
        pitch_m=arguments.pitch,
        width_m=arguments.width,
        thickness_m=arguments.thickness,
        gap_m=arguments.gap,
        length_m=strip_length_m,
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
    corrected_map_path = prefix.with_name(
        prefix.name + "_corrected_kickmap.dat"
    )
    metrics_path = prefix.with_name(prefix.name + "_metrics.txt")

    np.savetxt(
        currents_path,
        np.column_stack([np.arange(currents_a.size), currents_a]),
        delimiter=",",
        header="index,current_A",
        comments="",
        fmt=["%d", "%.10e"],
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

    raw_3d_figure = plot_kickmaps_3d(
        kick_map.x_m,
        kick_map.y_m,
        kick_map.theta_x_urad,
        kick_map.theta_y_urad,
        prefix.with_name(prefix.name + "_raw_kickmaps_3d.png"),
        theta_x_limit,
        theta_y_limit,
        x_range_mm=arguments.x_range_mm,
        y_range_mm=arguments.y_range_mm,
        keep_open=arguments.show,
    )
    corrected_3d_figure = plot_kickmaps_3d(
        kick_map.x_m,
        kick_map.y_m,
        theta_x_corrected_urad,
        theta_y_corrected_urad,
        prefix.with_name(prefix.name + "_corrected_kickmaps_3d.png"),
        theta_x_limit,
        theta_y_limit,
        x_range_mm=arguments.x_range_mm,
        y_range_mm=arguments.y_range_mm,
        keep_open=arguments.show,
    )

    if arguments.show:
        plt.show()
        if raw_3d_figure is not None:
            plt.close(raw_3d_figure)
        if corrected_3d_figure is not None:
            plt.close(corrected_3d_figure)

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
        f"input_file={arguments.kick}\n"
        f"grid_shape={kick_map.theta_x_urad.shape}\n"
        f"length_m={kick_map.length_m:.12g}\n"
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
        f"max_abs_current_A={np.max(np.abs(currents_a)):.12g}\n"
    )
    metrics_path.write_text(metrics, encoding="utf-8", newline="\n")

    print(metrics, end="")


if __name__ == "__main__":
    main()
