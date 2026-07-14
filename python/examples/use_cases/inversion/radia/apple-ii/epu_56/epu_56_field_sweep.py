#!/usr/bin/env python3
"""Sweep APPLE-II magnetic gap with RADIA and plot vertical field versus gap.

This script is intended to live in ``apple-ii/epu_56/``.  When the checked-out
source tree is present, it prepends ``../python`` to ``sys.path`` so that the
local ``apple_ii`` package is used.  The RADIA extension must still be importable
by the selected Python interpreter, for example through the user's PYTHONPATH.

The default gap points match the sampling used for a Figure-6-style comparison:
20, 25, 30, 40, 50, 60, 70, 80, 90, and 100 mm.

For each gap the script:
  1. rebuilds the four-array device at the requested phase;
  2. samples the complete on-axis field;
  3. extracts central-region peak |By|;
  4. fits the first harmonic of Bx and By;
  5. writes CSV and JSON summaries;
  6. writes a PNG plot of vertical field versus gap.

The primary plotted quantity is the directly sampled central peak |By|.  No
empirical or fitted Daresbury values are inserted into the RADIA results.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Sequence


DEFAULT_GAPS_MM = (20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0)


def positive_float(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value <= 0.0:
        raise argparse.ArgumentTypeError("value must be finite and > 0")
    return value


def nonnegative_float(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value < 0.0:
        raise argparse.ArgumentTypeError("value must be finite and >= 0")
    return value


def integer_at_least(minimum: int):
    def convert(text: str) -> int:
        value = int(text)
        if value < minimum:
            raise argparse.ArgumentTypeError(f"value must be >= {minimum}")
        return value

    return convert


def build_parser(script_dir: Path) -> argparse.ArgumentParser:
    output_dir = script_dir / "runs" / "phase_0"

    parser = argparse.ArgumentParser(
        description="Sweep APPLE-II magnetic gap with RADIA at fixed phase"
    )

    gap_group = parser.add_mutually_exclusive_group()
    gap_group.add_argument(
        "--gaps",
        type=positive_float,
        nargs="+",
        metavar="MM",
        help=(
            "explicit gap values in mm; default: "
            + " ".join(f"{value:g}" for value in DEFAULT_GAPS_MM)
        ),
    )
    gap_group.add_argument(
        "--gap-range",
        type=positive_float,
        nargs=3,
        metavar=("MIN_MM", "MAX_MM", "STEP_MM"),
        help="inclusive uniform gap range",
    )

    parser.add_argument("--period", type=positive_float, default=56.0)
    parser.add_argument("--periods", type=integer_at_least(1), default=17)
    parser.add_argument("--width", type=positive_float, default=40.0)
    parser.add_argument("--height", type=positive_float, default=40.0)
    parser.add_argument("--remanence", type=positive_float, default=1.25)
    parser.add_argument(
        "--motion",
        choices=("elliptical", "opposite"),
        default="elliptical",
    )
    parser.add_argument("--phase-mm", type=float, default=0.0)
    parser.add_argument(
        "--inter-array-gap-x-mm",
        type=nonnegative_float,
        default=0.5,
    )
    parser.add_argument(
        "--end-block-fraction",
        type=nonnegative_float,
        default=0.5,
    )
    parser.add_argument(
        "--end-block-length-mm",
        type=positive_float,
        default=6.95,
    )
    parser.add_argument(
        "--end-magnetization-fraction",
        type=nonnegative_float,
        default=1.0,
    )
    parser.add_argument(
        "--end-clearance-mm",
        type=nonnegative_float,
        default=0.0,
    )
    parser.add_argument("--samples", type=integer_at_least(8), default=2401)
    parser.add_argument(
        "--central-periods",
        type=integer_at_least(1),
        default=12,
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=output_dir / "field_sweep.csv",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=output_dir / "field_sweep.json",
    )
    parser.add_argument(
        "--png",
        type=Path,
        default=output_dir / "field_sweep.png",
    )
    parser.add_argument(
        "--profiles-dir",
        type=Path,
        default=None,
        help="optional directory for one full on-axis field CSV per gap",
    )
    parser.add_argument(
        "--dpi",
        type=integer_at_least(72),
        default=300,
    )
    return parser


def resolve_gaps(args: argparse.Namespace, parser: argparse.ArgumentParser) -> list[float]:
    if args.gaps is not None:
        gaps = list(args.gaps)
    elif args.gap_range is not None:
        start, stop, step = args.gap_range
        if stop < start:
            parser.error("--gap-range MAX_MM must be >= MIN_MM")

        count = int(math.floor((stop - start) / step + 1.0e-12)) + 1
        gaps = [start + index * step for index in range(count)]

        tolerance = 1.0e-10 * max(1.0, abs(stop))
        if gaps[-1] < stop - tolerance:
            gaps.append(stop)
        elif abs(gaps[-1] - stop) <= tolerance:
            gaps[-1] = stop
    else:
        gaps = list(DEFAULT_GAPS_MM)

    if not gaps:
        parser.error("the gap sweep is empty")

    # RADIA calculations are deterministic; duplicate gaps only waste time.
    return sorted(set(float(value) for value in gaps))


def ensure_local_source_on_path(script_dir: Path) -> Path | None:
    source_dir = script_dir.parent / "python"
    package_dir = source_dir / "apple_ii"
    if package_dir.is_dir():
        sys.path.insert(0, str(source_dir))
        return source_dir
    return None


def import_runtime_dependencies() -> tuple[Any, Any, Any, Any, Any]:
    try:
        import numpy as np
    except ImportError as exc:
        raise SystemExit("NumPy is required: python3 -m pip install numpy") from exc

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "Matplotlib is required for the PNG output: "
            "python3 -m pip install matplotlib"
        ) from exc

    try:
        import radia
    except ImportError as exc:
        raise SystemExit(
            "RADIA is not importable. Add the directory containing radia.so "
            "to PYTHONPATH before running this script."
        ) from exc

    try:
        import apple_ii
        from apple_ii.geometry import Apple2Parameters, ArrayMotion, build_device
    except ImportError as exc:
        raise SystemExit(
            "The apple_ii package is not importable. Place this script in "
            "apple-ii/epu_56/ or add apple-ii/python to PYTHONPATH."
        ) from exc

    return np, plt, radia, apple_ii, (Apple2Parameters, ArrayMotion, build_device)


def first_harmonic(np: Any, z_mm: Any, values_t: Any, period_mm: float) -> dict[str, float]:
    wave_number = 2.0 * math.pi / period_mm
    matrix = np.column_stack(
        (
            np.cos(wave_number * z_mm),
            np.sin(wave_number * z_mm),
            np.ones_like(z_mm),
        )
    )
    coefficients, *_ = np.linalg.lstsq(matrix, values_t, rcond=None)
    cosine, sine, offset = (float(value) for value in coefficients)
    fitted = matrix @ coefficients

    amplitude = math.hypot(cosine, sine)
    phase_deg = math.degrees(math.atan2(-sine, cosine))
    residual_rms = float(np.sqrt(np.mean((values_t - fitted) ** 2)))
    centered_rms = float(np.sqrt(np.mean((values_t - np.mean(values_t)) ** 2)))
    if centered_rms == 0.0:
        relative_residual = 0.0 if residual_rms == 0.0 else math.inf
    else:
        relative_residual = residual_rms / centered_rms

    return {
        "amplitude_T": amplitude,
        "phase_deg": phase_deg,
        "offset_T": offset,
        "residual_rms_T": residual_rms,
        "relative_residual": relative_residual,
    }


def sample_axis(np: Any, radia: Any, object_id: int, z_mm: Any) -> tuple[Any, Any, Any]:
    bx = np.empty_like(z_mm, dtype=float)
    by = np.empty_like(z_mm, dtype=float)
    bz = np.empty_like(z_mm, dtype=float)

    for index, position_mm in enumerate(z_mm):
        field = radia.Fld(object_id, "b", [0.0, 0.0, float(position_mm)])
        if not isinstance(field, (list, tuple)) or len(field) < 3:
            raise RuntimeError(
                f"unexpected RADIA field result at z={position_mm:g} mm: {field!r}"
            )
        bx[index], by[index], bz[index] = (
            float(field[0]),
            float(field[1]),
            float(field[2]),
        )

    return bx, by, bz


def write_profile(
    path: Path,
    z_mm: Sequence[float],
    bx_t: Sequence[float],
    by_t: Sequence[float],
    bz_t: Sequence[float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("z_mm", "Bx_T", "By_T", "Bz_T"))
        writer.writerows(zip(z_mm, bx_t, by_t, bz_t))


def write_summary_csv(path: Path, records: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    os.replace(partial, path)


def write_json(path: Path, report: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    partial.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    os.replace(partial, path)


def write_plot(plt: Any, path: Path, records: Sequence[dict[str, Any]], dpi: int) -> None:
    gaps = [float(record["gap_mm"]) for record in records]
    peak_by = [float(record["peak_abs_By_T"]) for record in records]
    by1 = [float(record["By1_T"]) for record in records]

    path.parent.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(7.2, 5.2))
    axis = figure.add_subplot(1, 1, 1)
    axis.plot(gaps, peak_by, marker="x", label="central peak |By|")
    axis.plot(gaps, by1, marker="o", label="first-harmonic By1")
    axis.set_xlabel("Magnet Gap (mm)")
    axis.set_ylabel("Vertical Field (T)")
    axis.set_title("Prototype HU56 RADIA Gap Sweep\nFixed Phase")
    axis.grid(True, linewidth=0.35, alpha=0.35)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = build_parser(script_dir)
    args = parser.parse_args()

    if not math.isfinite(args.phase_mm):
        parser.error("--phase-mm must be finite")
    if args.central_periods > args.periods:
        parser.error("--central-periods must not exceed --periods")
    if args.end_block_fraction > 1.0:
        parser.error("--end-block-fraction must be <= 1")
    if args.end_magnetization_fraction > 1.0:
        parser.error("--end-magnetization-fraction must be <= 1")

    gaps = resolve_gaps(args, parser)
    source_dir = ensure_local_source_on_path(script_dir)
    np, plt, radia, apple_ii, geometry_api = import_runtime_dependencies()
    Apple2Parameters, ArrayMotion, build_device = geometry_api

    motion = (
        ArrayMotion.elliptical()
        if args.motion == "elliptical"
        else ArrayMotion.opposite_pair()
    )

    records: list[dict[str, Any]] = []
    half_scan_mm = 0.5 * (args.periods + 2) * args.period
    z_mm = np.linspace(-half_scan_mm, half_scan_mm, args.samples)
    central_half_length_mm = 0.5 * args.central_periods * args.period
    central_mask = np.abs(z_mm) <= central_half_length_mm + 1.0e-12

    print(f"apple_ii version: {apple_ii.__version__}")
    if source_dir is not None:
        print(f"apple_ii source: {source_dir}")
    print(f"gap points: {len(gaps)}")

    for index, gap_mm in enumerate(gaps, start=1):
        radia.UtiDelAll()

        parameters = Apple2Parameters(
            period_mm=args.period,
            n_periods=args.periods,
            gap_mm=gap_mm,
            row_width_x_mm=args.width,
            block_height_y_mm=args.height,
            remanence_t=args.remanence,
            phase_mm=args.phase_mm,
            motion=motion,
            inter_array_gap_x_mm=args.inter_array_gap_x_mm,
            end_block_fraction=args.end_block_fraction,
            end_block_length_mm=args.end_block_length_mm,
            end_magnetization_fraction=args.end_magnetization_fraction,
            end_clearance_mm=args.end_clearance_mm,
        )
        device = build_device(parameters, radia_module=radia)

        bx_t, by_t, bz_t = sample_axis(np, radia, device.object_id, z_mm)
        central_z = z_mm[central_mask]
        central_bx = bx_t[central_mask]
        central_by = by_t[central_mask]
        central_bz = bz_t[central_mask]

        bx1 = first_harmonic(np, central_z, central_bx, args.period)
        by1 = first_harmonic(np, central_z, central_by, args.period)

        record = {
            "gap_mm": gap_mm,
            "phase_mm": args.phase_mm,
            "peak_abs_Bx_T": float(np.max(np.abs(central_bx))),
            "peak_abs_By_T": float(np.max(np.abs(central_by))),
            "peak_abs_Bz_T": float(np.max(np.abs(central_bz))),
            "Bx1_T": bx1["amplitude_T"],
            "By1_T": by1["amplitude_T"],
            "Bx1_phase_deg": bx1["phase_deg"],
            "By1_phase_deg": by1["phase_deg"],
            "Bx1_relative_residual": bx1["relative_residual"],
            "By1_relative_residual": by1["relative_residual"],
            "total_B1_T": math.hypot(bx1["amplitude_T"], by1["amplitude_T"]),
        }
        records.append(record)

        print(
            f"[{index:02d}/{len(gaps):02d}] gap={gap_mm:8.3f} mm  "
            f"peak|By|={record['peak_abs_By_T']:.9g} T  "
            f"By1={record['By1_T']:.9g} T  "
            f"residual={record['By1_relative_residual']:.3g}"
        )

        if args.profiles_dir is not None:
            profile_name = f"field_gap_{gap_mm:g}mm.csv".replace(".", "p")
            write_profile(
                args.profiles_dir / profile_name,
                z_mm,
                bx_t,
                by_t,
                bz_t,
            )

    write_summary_csv(args.csv, records)

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "apple_ii_version": apple_ii.__version__,
        "quantity_plotted": "central sampled peak absolute vertical field",
        "parameters": {
            "period_mm": args.period,
            "n_periods": args.periods,
            "row_width_x_mm": args.width,
            "block_height_y_mm": args.height,
            "remanence_t": args.remanence,
            "motion": args.motion,
            "phase_mm": args.phase_mm,
            "inter_array_gap_x_mm": args.inter_array_gap_x_mm,
            "end_block_fraction": args.end_block_fraction,
            "end_block_length_mm": args.end_block_length_mm,
            "end_magnetization_fraction": args.end_magnetization_fraction,
            "end_clearance_mm": args.end_clearance_mm,
            "samples": args.samples,
            "central_periods": args.central_periods,
        },
        "gaps_mm": gaps,
        "records": records,
    }
    write_json(args.json, report)
    write_plot(plt, args.png, records, args.dpi)

    print(f"CSV:  {args.csv}")
    print(f"JSON: {args.json}")
    print(f"PNG:  {args.png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
