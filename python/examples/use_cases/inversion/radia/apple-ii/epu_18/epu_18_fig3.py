#!/usr/bin/env python3
"""Reproduce the WEP38 Figure-3-style K_rms versus gap tuning curve.

The script reuses the magnetic geometry and harmonic-fit helpers from a sibling
``epu_18_gen.py`` reference model. For each gap it rebuilds and solves the RADIA
model, samples the central periodic field, extracts the first-harmonic Bx/By
amplitudes, and evaluates

    Kx = 0.934 * Bx[T] * lambda_u[cm]
    Ky = 0.934 * By[T] * lambda_u[cm]
    K_rms = sqrt((Kx**2 + Ky**2) / 2)

The resonance wavelength uses the same convention as the reference model:

    lambda_r = lambda_u * (1 + K_rms**2) / (2 * gamma**2)

Outputs use the units shown in WEP38 Figure 3: gap in millimetres, K_rms
(dimensionless), and wavelength annotations in nanometres.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Sequence

import numpy as np

SCRIPT_REVISION = "wep38-fig3-krms-gap-v1"
DEFAULT_TARGET_WAVELENGTHS_NM = (5.8, 4.0, 3.5)
DEFAULT_TARGET_COLORS = ("red", "limegreen", "red")


def env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    value = float(raw) if raw is not None else default
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    return int(raw) if raw is not None else default


def env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    normalized = raw.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"{name} must be one of 0/1, true/false, yes/no, on/off")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate a WEP38 Figure-3-style K_rms versus gap curve.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("output_dir", nargs="?", type=Path, default=Path("runs/fig3"))
    parser.add_argument(
        "--model-script",
        type=Path,
        default=Path(__file__).resolve().with_name("epu_18_gen.py"),
        help="reference EPU-18 RADIA model to import",
    )
    parser.add_argument(
        "--mode",
        choices=("circular+", "circular-"),
        default="circular+",
        help="polarization mode used for the tuning curve",
    )
    parser.add_argument("--gap-min-mm", type=float, default=env_float("FIG3_GAP_MIN_MM", 1.5))
    parser.add_argument("--gap-max-mm", type=float, default=env_float("FIG3_GAP_MAX_MM", 5.0))
    parser.add_argument("--gap-points", type=int, default=env_int("FIG3_GAP_POINTS", 15))
    parser.add_argument("--period-mm", type=float, default=env_float("PERIOD_MM", 18.0))
    parser.add_argument("--periods", type=int, default=env_int("PERIODS", 110))
    parser.add_argument("--br-t", type=float, default=env_float("BR_T", 1.35))
    parser.add_argument("--energy-gev", type=float, default=env_float("ENERGY_GEV", 1.0))
    parser.add_argument(
        "--fit-periods",
        type=int,
        default=env_int("FIG3_FIT_PERIODS", 20),
        help="central periods used for first-harmonic fitting",
    )
    parser.add_argument(
        "--samples-per-period",
        type=int,
        default=env_int("FIG3_SAMPLES_PER_PERIOD", 64),
        help="on-axis longitudinal samples per fitted period",
    )
    parser.add_argument(
        "--target-wavelength-nm",
        type=float,
        nargs="+",
        default=list(DEFAULT_TARGET_WAVELENGTHS_NM),
        help="wavelengths annotated by dashed tuning guides",
    )
    parser.add_argument(
        "--target-colors",
        nargs="+",
        default=list(DEFAULT_TARGET_COLORS),
        help="matplotlib colors for target wavelength guides",
    )
    parser.add_argument(
        "--susceptibility",
        action=argparse.BooleanOptionalAction,
        default=env_bool("SUSCEPTIBILITY", True),
    )
    parser.add_argument(
        "--solve-precision",
        type=float,
        default=env_float("SOLVE_PRECISION", 1.0e-4),
    )
    parser.add_argument(
        "--solve-max-iterations",
        type=int,
        default=env_int("SOLVE_MAX_ITERATIONS", 1000),
    )
    parser.add_argument(
        "--radia-pythonpath",
        type=Path,
        default=(Path(os.environ["RADIA_PYTHONPATH"]) if os.environ.get("RADIA_PYTHONPATH") else None),
    )
    parser.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="reuse matching gap rows already present in the CSV",
    )
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.model_script.is_file():
        raise ValueError(f"model script not found: {args.model_script}")
    if not math.isfinite(args.gap_min_mm) or not math.isfinite(args.gap_max_mm):
        raise ValueError("gap limits must be finite")
    if args.gap_min_mm <= 0 or args.gap_max_mm <= args.gap_min_mm:
        raise ValueError("require 0 < gap_min_mm < gap_max_mm")
    if args.gap_points < 3:
        raise ValueError("gap_points must be at least 3")
    if args.period_mm <= 0 or args.periods < 1 or args.br_t <= 0 or args.energy_gev <= 0:
        raise ValueError("period, periods, Br and energy must be positive")
    if args.fit_periods < 1 or args.samples_per_period < 8:
        raise ValueError("fit_periods >= 1 and samples_per_period >= 8 are required")
    if args.solve_precision <= 0 or args.solve_max_iterations < 1:
        raise ValueError("invalid RADIA solve controls")
    if not args.target_wavelength_nm:
        raise ValueError("at least one target wavelength is required")
    if any((not math.isfinite(v) or v <= 0) for v in args.target_wavelength_nm):
        raise ValueError("target wavelengths must be finite and positive")
    if len(args.target_colors) not in {1, len(args.target_wavelength_nm)}:
        raise ValueError("target-colors must contain one color or one per wavelength")


def import_model_module(path: Path):
    spec = importlib.util.spec_from_file_location("epu_18_reference_model", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"unable to import model script: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    required = ("parse_args", "build_specs", "import_radia", "reshape_field", "first_harmonic")
    missing = [name for name in required if not hasattr(module, name)]
    if missing:
        raise RuntimeError(f"model script lacks required helpers: {', '.join(missing)}")
    return module


def make_base_config(module: Any, args: argparse.Namespace, output_dir: Path):
    argv = [
        args.mode,
        str(output_dir),
        "--period-mm", str(args.period_mm),
        "--periods", str(args.periods),
        "--gap-mm", "1.5",
        "--br-t", str(args.br_t),
        "--energy-gev", str(args.energy_gev),
        "--fit-periods", str(args.fit_periods),
        "--samples", str(args.fit_periods * args.samples_per_period + 1),
        "--no-kickmap",
        "--no-plot",
        "--solve-precision", str(args.solve_precision),
        "--solve-max-iterations", str(args.solve_max_iterations),
    ]
    argv.append("--susceptibility" if args.susceptibility else "--no-susceptibility")
    if args.radia_pythonpath is not None:
        argv.extend(["--radia-pythonpath", str(args.radia_pythonpath)])
    return module.parse_args(argv)


def compute_gap_point(module: Any, rad: Any, base_config: Any, gap_mm: float, args: argparse.Namespace) -> dict[str, Any]:
    config = replace(base_config, gap_mm=float(gap_mm))
    specs = module.build_specs(config)

    rad.UtiDelAll()
    objects = [
        rad.ObjThckPgn(
            spec.center_z_mm,
            spec.length_z_mm,
            [list(point) for point in spec.polygon_xy_mm],
            "z",
            list(spec.magnetization_t),
        )
        for spec in specs
    ]
    model = rad.ObjCnt(objects)

    material_handle = None
    solve_result = None
    if args.susceptibility:
        material_handle = int(rad.MatStd("NdFeB", args.br_t))
        for obj in objects:
            rad.MatApl(obj, material_handle)
        solve_result = rad.Solve(model, args.solve_precision, args.solve_max_iterations)

    fit_periods = min(args.fit_periods, args.periods)
    half_width_mm = 0.5 * fit_periods * args.period_mm
    sample_count = fit_periods * args.samples_per_period + 1
    z_mm = np.linspace(-half_width_mm, half_width_mm, sample_count)
    points = [[0.0, 0.0, float(z)] for z in z_mm]
    field = module.reshape_field(rad.Fld(model, "b", points), sample_count)
    if not np.all(np.isfinite(field)):
        raise RuntimeError(f"RADIA returned non-finite fields at gap={gap_mm:g} mm")

    bx_fit = module.first_harmonic(z_mm, field[:, 0], args.period_mm, half_width_mm)
    by_fit = module.first_harmonic(z_mm, field[:, 1], args.period_mm, half_width_mm)
    period_cm = args.period_mm / 10.0
    kx = 0.934 * bx_fit["amplitude_t"] * period_cm
    ky = 0.934 * by_fit["amplitude_t"] * period_cm
    krms = math.sqrt((kx * kx + ky * ky) / 2.0)

    rest_energy_gev = 0.00051099895
    gamma = args.energy_gev / rest_energy_gev
    resonance_nm = (
        args.period_mm * 1.0e-3 * (1.0 + krms * krms) / (2.0 * gamma * gamma) * 1.0e9
    )

    return {
        "gap_mm": float(gap_mm),
        "Bx_amplitude_T": float(bx_fit["amplitude_t"]),
        "By_amplitude_T": float(by_fit["amplitude_t"]),
        "Bx_fit_residual_rms_T": float(bx_fit["residual_rms_t"]),
        "By_fit_residual_rms_T": float(by_fit["residual_rms_t"]),
        "Kx": float(kx),
        "Ky": float(ky),
        "K_rms": float(krms),
        "resonance_wavelength_nm": float(resonance_nm),
        "within_reference_model_nominal_gap_range": bool(
            gap_mm <= float(getattr(module, "MAX_GAP_MM", gap_mm)) + 1.0e-12
        ),
        "radia_model_handle": int(model),
        "radia_material_handle": material_handle,
        "radia_solve_result": solve_result,
    }


def load_existing_rows(path: Path) -> dict[float, dict[str, Any]]:
    if not path.is_file():
        return {}
    rows: dict[float, dict[str, Any]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        for raw in reader:
            try:
                gap = float(raw["gap_mm"])
                rows[round(gap, 12)] = {
                    "gap_mm": gap,
                    "Bx_amplitude_T": float(raw["Bx_amplitude_T"]),
                    "By_amplitude_T": float(raw["By_amplitude_T"]),
                    "Bx_fit_residual_rms_T": float(raw["Bx_fit_residual_rms_T"]),
                    "By_fit_residual_rms_T": float(raw["By_fit_residual_rms_T"]),
                    "Kx": float(raw["Kx"]),
                    "Ky": float(raw["Ky"]),
                    "K_rms": float(raw["K_rms"]),
                    "resonance_wavelength_nm": float(raw["resonance_wavelength_nm"]),
                    "within_reference_model_nominal_gap_range": raw[
                        "within_reference_model_nominal_gap_range"
                    ].strip().lower() in {"1", "true", "yes"},
                }
            except (KeyError, TypeError, ValueError):
                return {}
    return rows


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    fieldnames = [
        "gap_mm",
        "Bx_amplitude_T",
        "By_amplitude_T",
        "Bx_fit_residual_rms_T",
        "By_fit_residual_rms_T",
        "Kx",
        "Ky",
        "K_rms",
        "resonance_wavelength_nm",
        "within_reference_model_nominal_gap_range",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row[name] for name in fieldnames})


def target_krms(wavelength_nm: float, period_mm: float, energy_gev: float) -> float:
    rest_energy_gev = 0.00051099895
    gamma = energy_gev / rest_energy_gev
    value = 2.0 * gamma * gamma * wavelength_nm * 1.0e-9 / (period_mm * 1.0e-3) - 1.0
    if value < 0:
        raise ValueError(
            f"target wavelength {wavelength_nm:g} nm is below the K=0 resonance"
        )
    return math.sqrt(value)


def interpolate_target_gap(gaps_mm: np.ndarray, krms: np.ndarray, target: float) -> float | None:
    order = np.argsort(gaps_mm)
    gaps = gaps_mm[order]
    kvals = krms[order]
    if target > float(np.max(kvals)) or target < float(np.min(kvals)):
        return None
    # K_rms decreases with gap. Reverse both arrays so np.interp sees an
    # increasing abscissa.
    return float(np.interp(target, kvals[::-1], gaps[::-1]))


def smooth_curve(gaps_mm: np.ndarray, krms: np.ndarray) -> tuple[np.ndarray, np.ndarray, str]:
    dense_gap = np.linspace(float(np.min(gaps_mm)), float(np.max(gaps_mm)), 500)
    try:
        from scipy.interpolate import PchipInterpolator  # type: ignore[import-not-found]
    except ImportError:
        dense_k = np.interp(dense_gap, gaps_mm, krms)
        return dense_gap, dense_k, "linear interpolation"
    interpolator = PchipInterpolator(gaps_mm, krms)
    return dense_gap, np.asarray(interpolator(dense_gap), dtype=float), "PCHIP"


def create_plot(
    output_dir: Path,
    gaps_mm: np.ndarray,
    krms: np.ndarray,
    target_records: Sequence[dict[str, Any]],
) -> tuple[str, str]:
    import matplotlib.pyplot as plt

    dense_gap, dense_k, interpolation_method = smooth_curve(gaps_mm, krms)
    figure, axis = plt.subplots(figsize=(7.2, 4.8))
    axis.plot(dense_gap, dense_k, linewidth=2.2, color="#d88900")
    axis.scatter(gaps_mm, krms, s=18, color="#d88900", zorder=3)

    gap_min = float(np.min(gaps_mm))
    for index, record in enumerate(target_records):
        target_gap = record["interpolated_gap_mm"]
        if target_gap is None:
            continue
        target_k = float(record["target_K_rms"])
        color = str(record["color"])
        axis.hlines(target_k, gap_min, target_gap, colors=color, linestyles="--", linewidth=1.1)
        axis.vlines(target_gap, 0.0, target_k, colors=color, linestyles="--", linewidth=1.1)
        label = rf"$\lambda={record['wavelength_nm']:g}\,\mathrm{{nm}}$"
        if index == 0:
            axis.text(gap_min + 0.06, max(0.08, 0.30 * target_k), label, color=color, fontsize=10)
        else:
            axis.text(target_gap + 0.04, target_k + 0.035, label, color=color, fontsize=10)

    axis.set_xlim(gap_min, float(np.max(gaps_mm)))
    axis.set_ylim(0.0, max(1.3, 1.06 * float(np.max(krms))))
    axis.set_xlabel("gap [mm]")
    axis.set_ylabel(r"$K_{\mathrm{rms}}$")
    axis.grid(False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(direction="out")
    figure.tight_layout()

    png_path = output_dir / "fig3_krms_vs_gap.png"
    pdf_path = output_dir / "fig3_krms_vs_gap.pdf"
    figure.savefig(png_path, dpi=220)
    figure.savefig(pdf_path)
    plt.close(figure)
    return str(png_path), interpolation_method


def script_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        validate_args(args)
        output_dir = args.output_dir.resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        model_script = args.model_script.resolve()
        module = import_model_module(model_script)
        gaps_mm = np.linspace(args.gap_min_mm, args.gap_max_mm, args.gap_points)

        colors = list(args.target_colors)
        if len(colors) == 1:
            colors *= len(args.target_wavelength_nm)
        target_records = [
            {
                "wavelength_nm": float(wavelength),
                "target_K_rms": target_krms(float(wavelength), args.period_mm, args.energy_gev),
                "color": colors[index],
            }
            for index, wavelength in enumerate(args.target_wavelength_nm)
        ]

        print("WEP38 Figure 3 — K_rms versus gap")
        print(f"  mode:                 {args.mode}")
        print(f"  gap scan:             {args.gap_min_mm:g} ... {args.gap_max_mm:g} mm, points={args.gap_points}")
        print(f"  period / periods:     {args.period_mm:g} mm / {args.periods}")
        print(f"  energy:               {args.energy_gev:g} GeV")
        print(f"  fit sampling:         {args.fit_periods} periods, {args.samples_per_period} points/period")
        print(f"  susceptibility solve: {args.susceptibility}")
        print(f"  model script:         {model_script}")
        print(f"  output directory:     {output_dir}")
        print(f"  dry run:              {args.dry_run}")

        if args.dry_run:
            payload = {
                "implementation": {"revision": SCRIPT_REVISION, "script_sha256": script_sha256()},
                "planned_gap_values_mm": [float(value) for value in gaps_mm],
                "target_wavelengths": target_records,
                "units": {"gap": "mm", "K_rms": "dimensionless", "wavelength": "nm"},
            }
            (output_dir / "fig3_manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
            print(json.dumps(payload, indent=2))
            return 0

        rad = module.import_radia(args.radia_pythonpath.resolve() if args.radia_pythonpath else None)
        base_config = make_base_config(module, args, output_dir)
        csv_path = output_dir / "fig3_krms_vs_gap.csv"
        existing = load_existing_rows(csv_path) if args.resume else {}
        results: list[dict[str, Any]] = []

        for index, gap_mm in enumerate(gaps_mm, start=1):
            key = round(float(gap_mm), 12)
            if key in existing:
                row = existing[key]
                print(f"  [{index:02d}/{len(gaps_mm):02d}] gap={gap_mm:.6g} mm — reused, K_rms={row['K_rms']:.6g}")
            else:
                row = compute_gap_point(module, rad, base_config, float(gap_mm), args)
                print(f"  [{index:02d}/{len(gaps_mm):02d}] gap={gap_mm:.6g} mm — K_rms={row['K_rms']:.6g}, lambda={row['resonance_wavelength_nm']:.6g} nm")
            results.append(row)
            write_csv(csv_path, results)

        results.sort(key=lambda row: row["gap_mm"])
        write_csv(csv_path, results)
        gaps = np.asarray([row["gap_mm"] for row in results], dtype=float)
        kvals = np.asarray([row["K_rms"] for row in results], dtype=float)
        for record in target_records:
            record["interpolated_gap_mm"] = interpolate_target_gap(gaps, kvals, record["target_K_rms"])

        plot_path = None
        interpolation_method = None
        if not args.no_plot:
            plot_path, interpolation_method = create_plot(output_dir, gaps, kvals, target_records)

        analysis = {
            "implementation": {
                "revision": SCRIPT_REVISION,
                "script_path": str(Path(__file__).resolve()),
                "script_sha256": script_sha256(),
                "model_script": str(model_script),
                "model_script_sha256": hashlib.sha256(model_script.read_bytes()).hexdigest(),
            },
            "mode": args.mode,
            "period_mm": args.period_mm,
            "periods": args.periods,
            "energy_gev": args.energy_gev,
            "remanence_T": args.br_t,
            "gap_scan": {
                "minimum_mm": args.gap_min_mm,
                "maximum_mm": args.gap_max_mm,
                "points": args.gap_points,
                "values_mm": [float(value) for value in gaps],
                "note": (
                    f"The imported reference model declares a nominal maximum gap of "
                    f"{getattr(module, 'MAX_GAP_MM', 'unspecified')} mm; larger plotted gaps are magnetic-model extrapolation."
                ),
            },
            "harmonic_fit": {
                "fit_periods": args.fit_periods,
                "samples_per_period": args.samples_per_period,
                "K_definition": "Kx=0.934*Bx[T]*period[cm], Ky=0.934*By[T]*period[cm], K_rms=sqrt((Kx^2+Ky^2)/2)",
            },
            "resonance": {
                "definition": "lambda_r=lambda_u*(1+K_rms^2)/(2*gamma^2)",
                "wavelength_units": "nm",
            },
            "target_wavelength_guides": target_records,
            "units": {"gap": "mm", "K_rms": "dimensionless", "wavelength": "nm"},
            "results": results,
            "outputs": {
                "csv": str(csv_path),
                "plot_png": plot_path,
                "plot_pdf": str(output_dir / "fig3_krms_vs_gap.pdf") if plot_path else None,
                "interpolation_method": interpolation_method,
            },
        }
        analysis_path = output_dir / "fig3_analysis.json"
        analysis_path.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(analysis, indent=2))
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
