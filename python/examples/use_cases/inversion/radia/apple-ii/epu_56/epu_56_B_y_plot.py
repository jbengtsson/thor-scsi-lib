#!/usr/bin/env python3
"""Generate a HU56 on-axis vertical-field plot with RADIA.

The script is intended to live in ``apple-ii/epu_56/``.  It uses the checked-out
``apple_ii`` source package from ``../python`` when available and computes the
field directly with RADIA.

Default prototype parameters:
  period                         56 mm
  periods                        17
  magnetic gap                   21 mm
  block cross-section            40 x 40 mm
  remanence                      1.25 T
  horizontal inter-array gap     0.5 mm
  end-block length               6.95 mm
  end magnetization fraction     1.0
  row phase                      0 mm

Outputs:
  runs/fig5/on_axis_vertical_field.png
  runs/fig5/on_axis_vertical_field.csv
  runs/fig5/on_axis_vertical_field.json
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Sequence


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
    output_dir = script_dir / "runs" / "phase_0_refined_prototype"

    parser = argparse.ArgumentParser(
        description="Generate a HU56 on-axis vertical-field plot"
    )
    parser.add_argument("--period", type=positive_float, default=56.0)
    parser.add_argument("--periods", type=integer_at_least(1), default=17)
    parser.add_argument("--gap", type=positive_float, default=21.0)
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
        "--scan-padding-periods",
        type=nonnegative_float,
        default=1.0,
        help="extra longitudinal range on each end, in periods",
    )
    parser.add_argument(
        "--png",
        type=Path,
        default=output_dir / "on_axis_vertical_field.png",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=output_dir / "on_axis_vertical_field.csv",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=output_dir / "on_axis_vertical_field.json",
    )
    parser.add_argument("--dpi", type=integer_at_least(72), default=300)
    parser.add_argument(
        "--show",
        action="store_true",
        help="open an interactive Matplotlib window after saving",
    )
    return parser


def ensure_local_source_on_path(script_dir: Path) -> Path | None:
    source_dir = script_dir.parent / "python"
    if (source_dir / "apple_ii").is_dir():
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
            "Matplotlib is required: python3 -m pip install matplotlib"
        ) from exc

    try:
        import radia
    except ImportError as exc:
        raise SystemExit(
            "RADIA is not importable. Add the directory containing radia.so "
            "to PYTHONPATH."
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


def sample_on_axis(
    np: Any,
    radia: Any,
    object_id: int,
    z_mm: Sequence[float],
) -> tuple[Any, Any, Any]:
    bx = np.empty(len(z_mm), dtype=float)
    by = np.empty(len(z_mm), dtype=float)
    bz = np.empty(len(z_mm), dtype=float)

    for index, position_mm in enumerate(z_mm):
        field = radia.Fld(object_id, "b", [0.0, 0.0, float(position_mm)])
        if not isinstance(field, (list, tuple)) or len(field) < 3:
            raise RuntimeError(
                f"unexpected RADIA field at z={position_mm:g} mm: {field!r}"
            )
        bx[index] = float(field[0])
        by[index] = float(field[1])
        bz[index] = float(field[2])

    return bx, by, bz


def write_csv(
    path: Path,
    z_mm: Sequence[float],
    bx_t: Sequence[float],
    by_t: Sequence[float],
    bz_t: Sequence[float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("z_mm", "Bx_T", "By_T", "Bz_T"))
        writer.writerows(zip(z_mm, bx_t, by_t, bz_t))
    os.replace(temporary, path)


def write_json(path: Path, report: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = build_parser(script_dir)
    args = parser.parse_args()

    if not math.isfinite(args.phase_mm):
        parser.error("--phase-mm must be finite")
    if args.end_block_fraction > 1.0:
        parser.error("--end-block-fraction must be <= 1")
    if args.end_magnetization_fraction > 1.0:
        parser.error("--end-magnetization-fraction must be <= 1")

    source_dir = ensure_local_source_on_path(script_dir)
    np, plt, radia, apple_ii, geometry_api = import_runtime_dependencies()
    Apple2Parameters, ArrayMotion, build_device = geometry_api

    motion = (
        ArrayMotion.elliptical()
        if args.motion == "elliptical"
        else ArrayMotion.opposite_pair()
    )

    parameters = Apple2Parameters(
        period_mm=args.period,
        n_periods=args.periods,
        gap_mm=args.gap,
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

    radia.UtiDelAll()
    device = build_device(parameters, radia_module=radia)

    half_body_length_mm = 0.5 * args.periods * args.period
    half_scan_mm = half_body_length_mm + args.scan_padding_periods * args.period
    z_mm = np.linspace(-half_scan_mm, half_scan_mm, args.samples)

    bx_t, by_t, bz_t = sample_on_axis(
        np,
        radia,
        device.object_id,
        z_mm,
    )

    peak_abs_by_t = float(np.max(np.abs(by_t)))
    peak_abs_bx_t = float(np.max(np.abs(bx_t)))
    peak_abs_bz_t = float(np.max(np.abs(bz_t)))

    write_csv(args.csv, z_mm, bx_t, by_t, bz_t)

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "apple_ii_version": apple_ii.__version__,
        "apple_ii_source": str(source_dir) if source_dir is not None else None,
        "figure": "On-axis vertical magnetic field",
        "parameters": {
            "period_mm": args.period,
            "n_periods": args.periods,
            "gap_mm": args.gap,
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
            "scan_padding_periods": args.scan_padding_periods,
        },
        "results": {
            "peak_abs_Bx_T": peak_abs_bx_t,
            "peak_abs_By_T": peak_abs_by_t,
            "peak_abs_Bz_T": peak_abs_bz_t,
        },
        "outputs": {
            "png": str(args.png),
            "csv": str(args.csv),
            "json": str(args.json),
        },
    }
    write_json(args.json, report)

    args.png.parent.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(8.0, 5.2))
    axis = figure.add_subplot(1, 1, 1)
    axis.plot(z_mm, by_t, linewidth=1.15, label="RADIA prototype")
    axis.axhline(0.0, linewidth=0.7)
    axis.axvline(0.0, linewidth=0.7)
    axis.set_xlabel("Longitudinal position (mm)")
    axis.set_ylabel("Vertical field, $B_y$ (T)")
    axis.set_title(
        "Prototype HU56: On-Axis Vertical Magnetic Field\n"
        f"Gap {args.gap:g} mm, phase {args.phase_mm:g} mm"
    )
    axis.set_xlim(float(z_mm[0]), float(z_mm[-1]))
    vertical_limit = 1.08 * peak_abs_by_t if peak_abs_by_t > 0.0 else 1.0
    axis.set_ylim(-vertical_limit, vertical_limit)
    axis.grid(True, linewidth=0.35, alpha=0.35)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(args.png, dpi=args.dpi, bbox_inches="tight")

    if args.show:
        plt.show()
    else:
        plt.close(figure)

    print(f"apple_ii version: {apple_ii.__version__}")
    if source_dir is not None:
        print(f"apple_ii source: {source_dir}")
    print(f"peak |Bx|: {peak_abs_bx_t:.9g} T")
    print(f"peak |By|: {peak_abs_by_t:.9g} T")
    print(f"peak |Bz|: {peak_abs_bz_t:.9g} T")
    print(f"PNG:  {args.png}")
    print(f"CSV:  {args.csv}")
    print(f"JSON: {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
