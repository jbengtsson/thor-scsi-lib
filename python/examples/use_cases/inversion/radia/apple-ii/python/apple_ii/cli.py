from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .analysis import (
    analyze_harmonics,
    evaluate_quality,
    harmonics_to_dict,
    polarization_metrics,
)
from .cli_common import (
    add_geometry_arguments,
    finite_float,
    integer_at_least,
    nonnegative_float,
    parameters_from_args,
    parameters_to_dict,
    quality_from_args,
)
from .field import field_integrals, sample_axis, write_csv
from .diagnostics import row_field_diagnostics
from .geometry import build_device
from ._util import require_radia


C_GEV_PER_TM = 0.299792458
ELECTRON_REST_ENERGY_EV = 510_998.95069


def _energy_to_brho(
    energy_ev: float,
    *,
    energy_kind: str,
    rest_energy_ev: float,
    charge_abs_e: float,
) -> dict[str, float | str]:
    """Convert relativistic beam energy to magnetic rigidity."""
    if not math.isfinite(energy_ev) or energy_ev <= 0:
        raise ValueError("beam energy must be finite and positive")
    if not math.isfinite(rest_energy_ev) or rest_energy_ev <= 0:
        raise ValueError("rest energy must be finite and positive")
    if not math.isfinite(charge_abs_e) or charge_abs_e <= 0:
        raise ValueError("absolute charge must be finite and positive")

    total_ev = energy_ev if energy_kind == "total" else energy_ev + rest_energy_ev
    if total_ev < rest_energy_ev:
        raise ValueError("total beam energy is below the particle rest energy")

    pc_ev = math.sqrt(max(0.0, total_ev * total_ev - rest_energy_ev * rest_energy_ev))
    if pc_ev <= 0:
        raise ValueError("beam momentum is zero")

    momentum_gev_c = pc_ev / 1.0e9
    brho_tm = momentum_gev_c / (C_GEV_PER_TM * charge_abs_e)
    return {
        "energy_input_ev": float(energy_ev),
        "energy_kind": energy_kind,
        "total_energy_ev": float(total_ev),
        "rest_energy_ev": float(rest_energy_ev),
        "momentum_eV_c": float(pc_ev),
        "momentum_GeV_c": float(momentum_gev_c),
        "charge_abs_e": float(charge_abs_e),
        "beam_rigidity_Tm": float(brho_tm),
        "beta": float(pc_ev / total_ev),
        "gamma": float(total_ev / rest_energy_ev),
    }


def _radia_field_xyz(radia, object_id: int, points_mm: list[list[float]]) -> np.ndarray:
    """Return Bx, By, Bz in tesla for a list of RADIA points in millimetres."""
    values = radia.Fld(object_id, "b", points_mm)
    array = np.asarray(values, dtype=float)

    if array.ndim == 1:
        if len(points_mm) == 1 and array.size == 3:
            array = array.reshape(1, 3)
        elif array.size == 3 * len(points_mm):
            array = array.reshape(len(points_mm), 3)

    if array.shape != (len(points_mm), 3):
        raise RuntimeError(
            "Unexpected RADIA field result shape "
            f"{array.shape}; expected ({len(points_mm)}, 3)"
        )
    return array


def _sample_transverse_field_grid(
    radia,
    object_id: int,
    *,
    x_min_mm: float,
    x_max_mm: float,
    x_points: int,
    y_min_mm: float,
    y_max_mm: float,
    y_points: int,
    z_min_mm: float,
    z_max_mm: float,
    z_points: int,
) -> pd.DataFrame:
    """Sample B(x,y,z) over a regular transverse grid."""
    xs = np.linspace(x_min_mm, x_max_mm, x_points)
    ys = np.linspace(y_min_mm, y_max_mm, y_points)
    zs = np.linspace(z_min_mm, z_max_mm, z_points)

    records: list[dict[str, float]] = []
    for y_mm in ys:
        for x_mm in xs:
            points = [[float(x_mm), float(y_mm), float(z_mm)] for z_mm in zs]
            fields = _radia_field_xyz(radia, object_id, points)
            for z_mm, (bx_t, by_t, bz_t) in zip(zs, fields):
                records.append(
                    {
                        "x_mm": float(x_mm),
                        "y_mm": float(y_mm),
                        "z_mm": float(z_mm),
                        "Bx_T": float(bx_t),
                        "By_T": float(by_t),
                        "Bz_T": float(bz_t),
                    }
                )
    return pd.DataFrame.from_records(records)


def _integrate_kicks(
    field_grid: pd.DataFrame,
    *,
    brho_tm: float,
    charge_sign: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Integrate fields and remove the on-axis first-integral steering term.

    The raw first field integrals are evaluated at every transverse point. The
    values at the grid point closest to (x, y) = (0, 0) are then subtracted
    from the complete map before conversion to angular kicks. With the default
    odd-numbered symmetric grid this reference point is exactly on axis.
    """
    rows: list[dict[str, float | int]] = []
    for (x_mm, y_mm), trace in field_grid.groupby(["x_mm", "y_mm"], sort=True):
        trace = trace.sort_values("z_mm")
        z_m = trace["z_mm"].to_numpy(float) * 1.0e-3
        if len(z_m) < 2 or not np.all(np.diff(z_m) > 0):
            raise ValueError(f"invalid z trace at x={x_mm:g} mm, y={y_mm:g} mm")

        int_bx_raw_tm = float(np.trapezoid(trace["Bx_T"].to_numpy(float), z_m))
        int_by_raw_tm = float(np.trapezoid(trace["By_T"].to_numpy(float), z_m))
        rows.append(
            {
                "x_mm": float(x_mm),
                "y_mm": float(y_mm),
                "x_m": float(x_mm) * 1.0e-3,
                "y_m": float(y_mm) * 1.0e-3,
                "int_Bx_raw_Tm": int_bx_raw_tm,
                "int_By_raw_Tm": int_by_raw_tm,
                "n_z": int(len(trace)),
                "z_min_mm": float(trace["z_mm"].iloc[0]),
                "z_max_mm": float(trace["z_mm"].iloc[-1]),
            }
        )

    result = pd.DataFrame(rows)
    if len(result) < 6:
        raise ValueError("at least six transverse grid points are required")

    radius2 = result["x_mm"].to_numpy(float) ** 2 + result["y_mm"].to_numpy(float) ** 2
    reference_index = int(np.argmin(radius2))
    reference = result.iloc[reference_index]

    x_ref_mm = float(reference["x_mm"])
    y_ref_mm = float(reference["y_mm"])
    x_step = float(np.min(np.diff(np.unique(result["x_mm"])))) if result["x_mm"].nunique() > 1 else 0.0
    y_step = float(np.min(np.diff(np.unique(result["y_mm"])))) if result["y_mm"].nunique() > 1 else 0.0
    tolerance_mm = max(abs(x_step), abs(y_step), 1.0) * 1.0e-9
    if abs(x_ref_mm) > tolerance_mm or abs(y_ref_mm) > tolerance_mm:
        raise ValueError(
            "transverse field grid does not contain the on-axis point (0, 0); "
            "use odd point counts and ranges that include zero"
        )

    int_bx_ref_tm = float(reference["int_Bx_raw_Tm"])
    int_by_ref_tm = float(reference["int_By_raw_Tm"])

    result["int_Bx_Tm"] = result["int_Bx_raw_Tm"] - int_bx_ref_tm
    result["int_By_Tm"] = result["int_By_raw_Tm"] - int_by_ref_tm
    result["theta_x_rad"] = -(charge_sign / brho_tm) * result["int_By_Tm"]
    result["theta_y_rad"] = +(charge_sign / brho_tm) * result["int_Bx_Tm"]

    subtraction = {
        "reference_x_mm": x_ref_mm,
        "reference_y_mm": y_ref_mm,
        "removed_int_Bx_Tm": int_bx_ref_tm,
        "removed_int_By_Tm": int_by_ref_tm,
        "removed_theta_x_rad": float(-(charge_sign / brho_tm) * int_by_ref_tm),
        "removed_theta_y_rad": float(+(charge_sign / brho_tm) * int_bx_ref_tm),
    }
    return result, subtraction


def _potential_matrix(x_m: np.ndarray, y_m: np.ndarray) -> np.ndarray:
    """Design matrix for kicks generated by a cubic scalar potential."""
    n = len(x_m)
    matrix = np.zeros((2 * n, 9), dtype=float)

    # W = c10*x + c01*y + c20*x² + c11*x*y + c02*y²
    #   + c30*x³ + c21*x²*y + c12*x*y² + c03*y³
    # theta_x = -dW/dx
    matrix[:n, 0] = -1.0
    matrix[:n, 2] = -2.0 * x_m
    matrix[:n, 3] = -y_m
    matrix[:n, 5] = -3.0 * x_m * x_m
    matrix[:n, 6] = -2.0 * x_m * y_m
    matrix[:n, 7] = -y_m * y_m

    # theta_y = -dW/dy
    matrix[n:, 1] = -1.0
    matrix[n:, 3] = -x_m
    matrix[n:, 4] = -2.0 * y_m
    matrix[n:, 6] = -x_m * x_m
    matrix[n:, 7] = -2.0 * x_m * y_m
    matrix[n:, 8] = -3.0 * y_m * y_m
    return matrix


def _fit_second_order_kick_map(kicks: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Fit a symplectic second-order kick map."""
    x_m = kicks["x_m"].to_numpy(float)
    y_m = kicks["y_m"].to_numpy(float)
    matrix = _potential_matrix(x_m, y_m)
    measured = np.concatenate(
        [
            kicks["theta_x_rad"].to_numpy(float),
            kicks["theta_y_rad"].to_numpy(float),
        ]
    )
    coefficients = np.linalg.lstsq(matrix, measured, rcond=None)[0]
    predicted = matrix @ coefficients
    n = len(kicks)

    fitted = kicks.copy()
    fitted["theta_x_fit_rad"] = predicted[:n]
    fitted["theta_y_fit_rad"] = predicted[n:]
    fitted["theta_x_residual_rad"] = (
        fitted["theta_x_rad"] - fitted["theta_x_fit_rad"]
    )
    fitted["theta_y_residual_rad"] = (
        fitted["theta_y_rad"] - fitted["theta_y_fit_rad"]
    )

    labels = ["c10", "c01", "c20", "c11", "c02", "c30", "c21", "c12", "c03"]
    residual_vector = np.hypot(
        fitted["theta_x_residual_rad"].to_numpy(float),
        fitted["theta_y_residual_rad"].to_numpy(float),
    )
    measured_vector = np.hypot(
        fitted["theta_x_rad"].to_numpy(float),
        fitted["theta_y_rad"].to_numpy(float),
    )
    scale = max(float(np.max(measured_vector)), np.finfo(float).eps)

    report = {
        "fit_mode": "cubic generating potential / second-order kicks",
        "symplectic_by_construction": True,
        "potential": (
            "W=c10*x+c01*y+c20*x^2+c11*x*y+c02*y^2+"
            "c30*x^3+c21*x^2*y+c12*x*y^2+c03*y^3"
        ),
        "kick_definition": "theta_x=-dW/dx; theta_y=-dW/dy",
        "potential_coefficients": dict(zip(labels, map(float, coefficients))),
        "statistics": {
            "n_transverse_points": int(n),
            "rms_theta_x_residual_rad": float(
                np.sqrt(np.mean(fitted["theta_x_residual_rad"].to_numpy(float) ** 2))
            ),
            "rms_theta_y_residual_rad": float(
                np.sqrt(np.mean(fitted["theta_y_residual_rad"].to_numpy(float) ** 2))
            ),
            "rms_vector_residual_rad": float(np.sqrt(np.mean(residual_vector**2))),
            "max_vector_residual_rad": float(np.max(residual_vector)),
            "max_measured_kick_rad": float(np.max(measured_vector)),
            "relative_rms_vector_residual": float(
                np.sqrt(np.mean(residual_vector**2)) / scale
            ),
        },
    }
    return fitted, report



def _write_on_axis_field_plot(
    field_csv: Path,
    output_png: Path,
    *,
    gap_mm: float,
    phase_mm: float,
    dpi: int,
    title: str,
) -> None:
    """Plot the modeled on-axis vertical field against longitudinal position."""
    frame = pd.read_csv(field_csv)
    required = {"z_mm", "By_T"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(
            "cannot generate field plot; missing CSV columns: "
            + ", ".join(missing)
        )

    z_mm = pd.to_numeric(frame["z_mm"], errors="coerce").to_numpy(float)
    by_t = pd.to_numeric(frame["By_T"], errors="coerce").to_numpy(float)
    valid = np.isfinite(z_mm) & np.isfinite(by_t)
    if np.count_nonzero(valid) < 2:
        raise ValueError("cannot generate field plot; fewer than two valid samples")

    z_mm = z_mm[valid]
    by_t = by_t[valid]
    order = np.argsort(z_mm)
    z_mm = z_mm[order]
    by_t = by_t[order]

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(z_mm, by_t, linewidth=1.15, label="Modeled")
    ax.axhline(0.0, linewidth=0.8)
    ax.axvline(0.0, linewidth=0.7, linestyle="--", alpha=0.7)
    ax.set_xlabel("Longitudinal position [mm]")
    ax.set_ylabel("Vertical field $B_y$ [T]")
    ax.set_title(title)
    ax.grid(True, linewidth=0.45, alpha=0.35)
    ax.legend(loc="best")
    ax.text(
        0.5,
        -0.22,
        f"On-axis vertical magnetic field at {gap_mm:g} mm gap "
        f"and {phase_mm:g} mm phase shift.",
        transform=ax.transAxes,
        ha="center",
        va="top",
    )
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _write_kick_map_plot(
    fitted: pd.DataFrame,
    output_png: Path,
    *,
    title: str,
    energy_ev: float,
    brho_tm: float,
    dpi: int,
) -> None:
    """Write a 3D fitted kick-magnitude surface with kick-direction arrows."""
    x_mm = fitted["x_mm"].to_numpy(float)
    y_mm = fitted["y_mm"].to_numpy(float)
    theta_x_urad = fitted["theta_x_fit_rad"].to_numpy(float) * 1.0e6
    theta_y_urad = fitted["theta_y_fit_rad"].to_numpy(float) * 1.0e6
    magnitude_urad = np.hypot(theta_x_urad, theta_y_urad)

    pivot = pd.DataFrame(
        {"x_mm": x_mm, "y_mm": y_mm, "kick": magnitude_urad}
    ).pivot(index="y_mm", columns="x_mm", values="kick")
    x_grid, y_grid = np.meshgrid(
        pivot.columns.to_numpy(float), pivot.index.to_numpy(float)
    )
    z_grid = pivot.to_numpy(float)

    fig = plt.figure(figsize=(10.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(
        x_grid,
        y_grid,
        z_grid,
        linewidth=0,
        antialiased=True,
        alpha=0.80,
    )

    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("|θ| [µrad]")
    ax.set_title(
        f"{title}\n"
        f"electron total energy={energy_ev:.6g} eV, Bρ={brho_tm:.6g} T·m"
    )
    ax.view_init(elev=28, azim=-55)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _add_kick_map_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("second-order kick map")
    group.add_argument(
        "--kickmap-prefix",
        default=None,
        help=(
            "output prefix for the shared transverse field grid, second-order kick map, "
            "JSON coefficients, and 3D PNG; omission disables generation"
        ),
    )
    group.add_argument(
        "--energy-ev",
        type=finite_float,
        default=None,
        help="electron beam energy in eV; required with --kickmap-prefix",
    )
    group.add_argument(
        "--energy-kind",
        choices=["total", "kinetic"],
        default="total",
        help="interpret --energy-ev as total or kinetic energy",
    )
    group.add_argument(
        "--charge-sign",
        type=int,
        choices=[-1, 1],
        default=-1,
        help="particle charge sign; default -1 for electrons",
    )
    group.add_argument(
        "--rest-energy-ev",
        type=finite_float,
        default=ELECTRON_REST_ENERGY_EV,
        help="particle rest energy m*c^2 in eV",
    )
    group.add_argument(
        "--charge-abs-e",
        type=finite_float,
        default=1.0,
        help="absolute particle charge in elementary-charge units",
    )
    field_group = parser.add_argument_group("transverse field map")
    field_group.add_argument(
        "--field-x-min-mm",
        type=finite_float,
        default=-10.0,
        help="minimum horizontal coordinate of the transverse field map",
    )
    field_group.add_argument(
        "--field-x-max-mm",
        type=finite_float,
        default=10.0,
        help="maximum horizontal coordinate of the transverse field map",
    )
    field_group.add_argument(
        "--field-x-points",
        type=integer_at_least(3),
        default=9,
        help="number of horizontal transverse field-map points",
    )
    field_group.add_argument(
        "--field-y-min-mm",
        type=finite_float,
        default=-5.0,
        help="minimum vertical coordinate of the transverse field map",
    )
    field_group.add_argument(
        "--field-y-max-mm",
        type=finite_float,
        default=5.0,
        help="maximum vertical coordinate of the transverse field map",
    )
    field_group.add_argument(
        "--field-y-points",
        type=integer_at_least(3),
        default=9,
        help="number of vertical transverse field-map points",
    )
    group.add_argument(
        "--kick-plot-dpi",
        type=integer_at_least(1),
        default=180,
    )
    group.add_argument(
        "--kick-plot-title",
        default="Second-order transverse kick map",
    )


def _generate_kick_map(
    args,
    *,
    radia,
    device,
    half_length_mm: float,
) -> dict:
    if args.energy_ev is None:
        raise ValueError("--energy-ev is required with --kickmap-prefix")
    if args.field_x_min_mm >= args.field_x_max_mm:
        raise ValueError("--field-x-min-mm must be smaller than --field-x-max-mm")
    if args.field_y_min_mm >= args.field_y_max_mm:
        raise ValueError("--field-y-min-mm must be smaller than --field-y-max-mm")

    beam = _energy_to_brho(
        args.energy_ev,
        energy_kind=args.energy_kind,
        rest_energy_ev=args.rest_energy_ev,
        charge_abs_e=args.charge_abs_e,
    )
    z_points = args.samples
    prefix = Path(args.kickmap_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    integrated_path = prefix.with_name(prefix.name + "_integrated.csv")
    plot_path = prefix.with_name(prefix.name + "_3d.png")

    field_grid = _sample_transverse_field_grid(
        radia,
        device.object_id,
        x_min_mm=args.field_x_min_mm,
        x_max_mm=args.field_x_max_mm,
        x_points=args.field_x_points,
        y_min_mm=args.field_y_min_mm,
        y_max_mm=args.field_y_max_mm,
        y_points=args.field_y_points,
        z_min_mm=-half_length_mm,
        z_max_mm=half_length_mm,
        z_points=z_points,
    )
    kicks, first_integral_subtraction = _integrate_kicks(
        field_grid,
        brho_tm=float(beam["beam_rigidity_Tm"]),
        charge_sign=args.charge_sign,
    )
    fitted, fit_report = _fit_second_order_kick_map(kicks)
    fitted.to_csv(integrated_path, index=False)

    report = {
        "scope": (
            "Second-order kick map fitted to the RADIA transverse field grid "
            "generated in the same device build"
        ),
        "beam": {
            **beam,
            "charge_sign": args.charge_sign,
        },
        "sampling": {
            "x_min_mm": args.field_x_min_mm,
            "x_max_mm": args.field_x_max_mm,
            "x_points": args.field_x_points,
            "y_min_mm": args.field_y_min_mm,
            "y_max_mm": args.field_y_max_mm,
            "y_points": args.field_y_points,
            "z_min_mm": -half_length_mm,
            "z_max_mm": half_length_mm,
            "z_points": z_points,
            "z_points_source": "--samples",
        },
        "physics": {
            "trajectory_direction": "+z",
            "theta_x": "-(charge_sign/B_rho)*corrected_integral(By dz)",
            "theta_y": "+(charge_sign/B_rho)*corrected_integral(Bx dz)",
            "first_integral_correction": (
                "subtract the on-axis raw first field integrals from every "
                "transverse point before kick conversion and fitting"
            ),
        },
        "first_integral_subtraction": first_integral_subtraction,
        "fit": fit_report,
        "outputs": {
            "integrated_csv": str(integrated_path),
            "plot_png": str(plot_path),
        },
        "plot": {
            "surface": "fitted kick magnitude only",
            "flow_arrows_included": False,
        },
    }
    _write_kick_map_plot(
        fitted,
        plot_path,
        title=args.kick_plot_title,
        energy_ev=args.energy_ev,
        brho_tm=float(beam["beam_rigidity_Tm"]),
        dpi=args.kick_plot_dpi,
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Build and analyze one APPLE-II operating point")
    add_geometry_arguments(parser)
    parser.add_argument("--phase-mm", type=finite_float, default=10.0)
    parser.add_argument("--samples", type=integer_at_least(8), default=1201)
    parser.add_argument("--central-periods", type=integer_at_least(1), default=6)
    parser.add_argument("--harmonics", type=integer_at_least(1), nargs="+", default=[1, 3, 5])
    parser.add_argument("--min-total-amplitude", type=nonnegative_float, default=1.0e-4)
    parser.add_argument("--max-relative-residual", type=nonnegative_float, default=0.10)
    parser.add_argument("--csv", default="apple2_field.csv")
    parser.add_argument("--json", default="apple2_analysis.json")
    parser.add_argument(
        "--diagnostics-json",
        default=None,
        help="optional geometry and per-row field diagnostic report",
    )
    parser.add_argument(
        "--field-plot",
        default="apple2_field.png",
        help="output PNG for the modeled on-axis vertical field",
    )
    parser.add_argument(
        "--field-plot-dpi",
        type=integer_at_least(1),
        default=180,
        help="resolution of the on-axis field plot",
    )
    parser.add_argument(
        "--field-plot-title",
        default="On-axis vertical magnetic field",
        help="title of the on-axis field plot",
    )
    parser.add_argument("--draw", action="store_true")
    _add_kick_map_arguments(parser)
    args = parser.parse_args()

    if 1 not in args.harmonics:
        parser.error("--harmonics must include 1")
    parameters = parameters_from_args(args, args.phase_mm)
    if args.central_periods > parameters.n_periods:
        parser.error("--central-periods must not exceed --periods")
    if args.kickmap_prefix and args.energy_ev is None:
        parser.error("--energy-ev is required with --kickmap-prefix")
    if args.energy_ev is not None and not args.kickmap_prefix:
        parser.error("--kickmap-prefix is required with --energy-ev")

    radia = require_radia()
    radia.UtiDelAll()
    device = build_device(parameters, radia_module=radia)
    half_length = 0.5 * (parameters.n_periods + 2) * parameters.period_mm

    rows = sample_axis(device, -half_length, half_length, args.samples)
    write_csv(rows, args.csv)
    _write_on_axis_field_plot(
        Path(args.csv),
        Path(args.field_plot),
        gap_mm=parameters.gap_mm,
        phase_mm=args.phase_mm,
        dpi=args.field_plot_dpi,
        title=args.field_plot_title,
    )
    harmonics = analyze_harmonics(
        rows,
        parameters.period_mm,
        harmonic_orders=args.harmonics,
        central_periods=args.central_periods,
    )
    bx1, by1 = harmonics[1]
    metrics = polarization_metrics(bx1, by1)
    quality = evaluate_quality(bx1, by1, quality_from_args(args))
    integrals = field_integrals(rows)
    report = {
        "parameters": parameters_to_dict(parameters),
        "harmonics": harmonics_to_dict(harmonics),
        "field_ellipse": metrics,
        "quality": asdict(quality),
        "field_integrals": asdict(integrals),
        "field_plot_png": str(args.field_plot),
    }

    if args.kickmap_prefix:
        report["kick_map"] = _generate_kick_map(
            args,
            radia=radia,
            device=device,
            half_length_mm=half_length,
        )

    Path(args.json).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if args.diagnostics_json:
        diagnostics = {
            "scope": (
                "Prototype geometry and per-row magnetic-field diagnostics; "
                "not a device-fidelity certification"
            ),
            "parameters": parameters_to_dict(parameters),
            "geometry": device.geometry_report(),
            "row_field_diagnostics": row_field_diagnostics(
                device,
                -half_length,
                half_length,
                args.samples,
                period_mm=parameters.period_mm,
                harmonic_orders=args.harmonics,
                central_periods=args.central_periods,
                combined_rows=rows,
                radia_module=radia,
            ),
        }
        Path(args.diagnostics_json).write_text(
            json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
        )

    print(f"Bx1={bx1.amplitude:.6g} T")
    print(f"By1={by1.amplitude:.6g} T")
    print(f"Bx1/By1={metrics['ratio']:.6g}")
    print(f"relative phase={metrics['phase_deg']:.6g} deg")
    print(f"magnetic circularity |s3|={metrics['circularity']:.6g}")
    print(f"quality gate={'PASS' if quality.passed else 'FAIL'}")
    print(f"field_plot_png={args.field_plot}")
    for reason in quality.reasons:
        print(f"  - {reason}")

    if args.kickmap_prefix:
        kick_report = report["kick_map"]
        print(
            "kick-map B*rho="
            f"{kick_report['beam']['beam_rigidity_Tm']:.12g} T*m"
        )
        correction = kick_report["first_integral_subtraction"]
        print(
            "removed first integrals: "
            f"IntBx={correction['removed_int_Bx_Tm']:.12g} T*m, "
            f"IntBy={correction['removed_int_By_Tm']:.12g} T*m"
        )
        for name, path in kick_report["outputs"].items():
            print(f"{name}={path}")

    if args.draw:
        radia.ObjDrwOpenGL(device.object_id)


if __name__ == "__main__":
    main()
