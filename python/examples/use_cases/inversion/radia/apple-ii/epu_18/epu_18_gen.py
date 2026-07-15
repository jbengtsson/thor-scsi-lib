#!/usr/bin/env python3
"""WEP38-informed 18 mm AQUA-style APPLE-X RADIA reference model.

This is a parameterized magnetic reference model, not manufacturer CAD or a
released construction design. The FEL2022 WEP38 paper does not publish every
polygon, locking-tooth, mover datum, tolerance, susceptibility, or exact
magnetization-sequence construction detail. The model therefore:

* uses a symmetric chamfered-polygon approximation for each 18 x 18 mm block;
* omits the locking tooth and mechanical holder;
* treats the published gap as radial girder retraction from the 1.5 mm datum;
* uses ideal fixed magnetization without susceptibility or demagnetization;
* implements a four-block Halbach cycle with radial/longitudinal directions;
* uses a symmetry convention giving vertical linear field at zero phase and
  nominal circular field near a quarter-period diagonal-row shift.

Examples:
    python3 aqua_apple_x_18.py circular+
    python3 aqua_apple_x_18.py linear-v ./runs/aqua_x18_linear
    python3 aqua_apple_x_18.py circular+ --gap-mm 2.7
    python3 aqua_apple_x_18.py circular+ --periods 3 --dry-run
    python3 aqua_apple_x_18.py circular+ --radia-pythonpath /path/to/radia
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


MODE_PHASES: dict[str, float] = {
    "linear-v": 0.0,
    "circular+": 0.25,
    "circular-": -0.25,
    "linear-h": 0.5,
}

MODE_RUN_NAMES: dict[str, str] = {
    "linear-v": "aqua_apple_x18_linear_v",
    "circular+": "aqua_apple_x18_circular_plus",
    "circular-": "aqua_apple_x18_circular_minus",
    "linear-h": "aqua_apple_x18_linear_h",
}


@dataclass(frozen=True)
class MagnetSpec:
    girder: str
    quadrant: int
    role: str
    sequence_index: int
    center_z_mm: float
    length_z_mm: float
    magnetization_t: tuple[float, float, float]
    polygon_xy_mm: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class ModelConfig:
    mode: str
    output_dir: Path
    period_mm: float
    periods: int
    gap_mm: float
    br_t: float
    block_xy_mm: float
    central_hole_mm: float
    chamfer_mm: float
    samples: int
    fit_periods: int
    energy_gev: float
    dry_run: bool
    no_plot: bool
    radia_pythonpath: Path | None

    @property
    def phase_fraction(self) -> float:
        return MODE_PHASES[self.mode]

    @property
    def phase_mm(self) -> float:
        return self.phase_fraction * self.period_mm


def env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    value = float(raw) if raw is not None else default
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    return int(raw) if raw is not None else default


def env_bool(name: str, default: bool = False) -> bool:
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
        description=(
            "Generate and evaluate a WEP38-informed 18 mm AQUA-style "
            "APPLE-X RADIA reference model."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "mode",
        nargs="?",
        choices=tuple(MODE_PHASES),
        default="circular+",
        help="polarization/phase mode",
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        type=Path,
        help="output directory; defaults to runs/<mode-name> beside this script",
    )
    parser.add_argument("--period-mm", type=float, default=env_float("PERIOD_MM", 18.0))
    parser.add_argument("--periods", type=int, default=env_int("PERIODS", 110))
    parser.add_argument("--gap-mm", type=float, default=env_float("GAP_MM", 1.5))
    parser.add_argument("--br-t", type=float, default=env_float("BR_T", 1.35))
    parser.add_argument(
        "--block-xy-mm", type=float, default=env_float("BLOCK_XY_MM", 18.0)
    )
    parser.add_argument(
        "--central-hole-mm",
        type=float,
        default=env_float("CENTRAL_HOLE_MM", 5.5),
    )
    parser.add_argument(
        "--chamfer-mm", type=float, default=env_float("CHAMFER_MM", 1.5)
    )
    parser.add_argument("--samples", type=int, default=env_int("SAMPLES", 4001))
    parser.add_argument(
        "--fit-periods", type=int, default=env_int("FIT_PERIODS", 20)
    )
    parser.add_argument(
        "--energy-gev", type=float, default=env_float("ENERGY_GEV", 1.0)
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=env_bool("DRY_RUN", False),
        help="write geometry and manifest without importing or running RADIA",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="skip on-axis field PNG generation",
    )
    parser.add_argument(
        "--radia-pythonpath",
        type=Path,
        default=Path(os.environ["RADIA_PYTHONPATH"])
        if os.environ.get("RADIA_PYTHONPATH")
        else None,
        help="directory containing the RADIA Python binding",
    )
    return parser


def parse_args(argv: Sequence[str] | None = None) -> ModelConfig:
    parser = build_parser()
    args = parser.parse_args(argv)

    script_dir = Path(__file__).resolve().parent
    output_dir = (
        args.output_dir
        if args.output_dir is not None
        else script_dir / "runs" / MODE_RUN_NAMES[args.mode]
    )

    config = ModelConfig(
        mode=args.mode,
        output_dir=output_dir.resolve(),
        period_mm=float(args.period_mm),
        periods=int(args.periods),
        gap_mm=float(args.gap_mm),
        br_t=float(args.br_t),
        block_xy_mm=float(args.block_xy_mm),
        central_hole_mm=float(args.central_hole_mm),
        chamfer_mm=float(args.chamfer_mm),
        samples=int(args.samples),
        fit_periods=int(args.fit_periods),
        energy_gev=float(args.energy_gev),
        dry_run=bool(args.dry_run),
        no_plot=bool(args.no_plot),
        radia_pythonpath=args.radia_pythonpath.resolve()
        if args.radia_pythonpath is not None
        else None,
    )
    validate_config(config)
    return config


def validate_config(config: ModelConfig) -> None:
    finite_positive = {
        "period_mm": config.period_mm,
        "br_t": config.br_t,
        "block_xy_mm": config.block_xy_mm,
        "central_hole_mm": config.central_hole_mm,
        "energy_gev": config.energy_gev,
    }
    for name, value in finite_positive.items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(config.gap_mm):
        raise ValueError("gap_mm must be finite")
    if not math.isfinite(config.chamfer_mm) or config.chamfer_mm < 0:
        raise ValueError("chamfer_mm must be finite and non-negative")
    if config.periods < 1:
        raise ValueError("periods must be at least 1")
    if config.samples < 21:
        raise ValueError("samples must be at least 21")
    if config.fit_periods < 1:
        raise ValueError("fit_periods must be at least 1")


def chamfered_rectangle(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    chamfer: float,
) -> tuple[tuple[float, float], ...]:
    """Return a counter-clockwise chamfered rectangle."""
    width = x_max - x_min
    height = y_max - y_min
    if width <= 0 or height <= 0:
        raise ValueError("invalid rectangle dimensions")
    if not (0 <= chamfer < 0.5 * min(width, height)):
        raise ValueError("chamfer_mm must be less than half the block size")
    if chamfer == 0:
        return (
            (x_min, y_min),
            (x_max, y_min),
            (x_max, y_max),
            (x_min, y_max),
        )
    return (
        (x_min + chamfer, y_min),
        (x_max - chamfer, y_min),
        (x_max, y_min + chamfer),
        (x_max, y_max - chamfer),
        (x_max - chamfer, y_max),
        (x_min + chamfer, y_max),
        (x_min, y_max - chamfer),
        (x_min, y_min + chamfer),
    )


def quadrant_polygon(
    sx: int,
    sy: int,
    block_xy_mm: float,
    central_hole_mm: float,
    gap_mm: float,
    min_gap_mm: float,
    chamfer_mm: float,
) -> tuple[tuple[float, float], ...]:
    """Create one quadrant block cross-section under the inferred gap datum."""
    radial_retraction = gap_mm - min_gap_mm
    shift = radial_retraction / math.sqrt(2.0)
    inner_x = 0.5 * central_hole_mm + shift
    inner_y = 0.5 * central_hole_mm + shift

    if sx > 0:
        x_min, x_max = inner_x, inner_x + block_xy_mm
    else:
        x_min, x_max = -inner_x - block_xy_mm, -inner_x
    if sy > 0:
        y_min, y_max = inner_y, inner_y + block_xy_mm
    else:
        y_min, y_max = -inner_y - block_xy_mm, -inner_y

    return chamfered_rectangle(x_min, x_max, y_min, y_max, chamfer_mm)


def halbach_direction(
    sequence_index: int,
    radial_unit: np.ndarray,
    row_sign: int,
    br_t: float,
) -> tuple[float, float, float]:
    """Four-block cycle: +radial, +z, -radial, -z."""
    ez = np.array([0.0, 0.0, 1.0])
    cycle = (radial_unit, ez, -radial_unit, -ez)
    vector = row_sign * br_t * cycle[sequence_index % 4]
    return float(vector[0]), float(vector[1]), float(vector[2])


def build_specs(config: ModelConfig) -> list[MagnetSpec]:
    min_gap_mm = 1.5
    if not (min_gap_mm <= config.gap_mm <= 4.0):
        raise ValueError("gap_mm must be within the WEP38 range 1.5 ... 4.0 mm")

    standard_length = config.period_mm / 4.0
    body_blocks = config.periods * 4

    # Figure-1 ordering: G1 top-right, G2 top-left, G3 bottom-left,
    # G4 bottom-right. G1 and G3 form the shifted diagonal pair.
    quadrants = (
        ("G1", +1, +1, +1, True),
        ("G2", -1, +1, +1, False),
        ("G3", -1, -1, -1, True),
        ("G4", +1, -1, -1, False),
    )

    specs: list[MagnetSpec] = []
    for q_index, (girder, sx, sy, row_sign, shifted) in enumerate(
        quadrants, start=1
    ):
        radial = np.array([sx, sy, 0.0], dtype=float) / math.sqrt(2.0)
        polygon = quadrant_polygon(
            sx=sx,
            sy=sy,
            block_xy_mm=config.block_xy_mm,
            central_hole_mm=config.central_hole_mm,
            gap_mm=config.gap_mm,
            min_gap_mm=min_gap_mm,
            chamfer_mm=config.chamfer_mm,
        )
        z_shift = config.phase_mm if shifted else 0.0
        body_left_edge = -0.5 * body_blocks * standard_length + z_shift

        for k in range(body_blocks):
            specs.append(
                MagnetSpec(
                    girder=girder,
                    quadrant=q_index,
                    role="body",
                    sequence_index=k,
                    center_z_mm=body_left_edge + (k + 0.5) * standard_length,
                    length_z_mm=standard_length,
                    magnetization_t=halbach_direction(
                        k, radial, row_sign, config.br_t
                    ),
                    polygon_xy_mm=polygon,
                )
            )

        # WEP38 Figure 4: body outward = 3/4, 1/2, 1/4.
        end_fractions = (0.75, 0.50, 0.25)

        cursor = body_left_edge
        for step, fraction in enumerate(end_fractions, start=1):
            length = fraction * standard_length
            sequence_index = -step
            specs.append(
                MagnetSpec(
                    girder=girder,
                    quadrant=q_index,
                    role=f"end-left-{fraction:g}",
                    sequence_index=sequence_index,
                    center_z_mm=cursor - 0.5 * length,
                    length_z_mm=length,
                    magnetization_t=halbach_direction(
                        sequence_index, radial, row_sign, config.br_t
                    ),
                    polygon_xy_mm=polygon,
                )
            )
            cursor -= length

        cursor = body_left_edge + body_blocks * standard_length
        for step, fraction in enumerate(end_fractions):
            length = fraction * standard_length
            sequence_index = body_blocks + step
            specs.append(
                MagnetSpec(
                    girder=girder,
                    quadrant=q_index,
                    role=f"end-right-{fraction:g}",
                    sequence_index=sequence_index,
                    center_z_mm=cursor + 0.5 * length,
                    length_z_mm=length,
                    magnetization_t=halbach_direction(
                        sequence_index, radial, row_sign, config.br_t
                    ),
                    polygon_xy_mm=polygon,
                )
            )
            cursor += length

    return specs


def write_geometry_csv(path: Path, specs: Sequence[MagnetSpec]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "girder",
                "quadrant",
                "role",
                "sequence_index",
                "center_z_mm",
                "length_z_mm",
                "mx_t",
                "my_t",
                "mz_t",
                "polygon_xy_mm_json",
            ]
        )
        for spec in specs:
            writer.writerow(
                [
                    spec.girder,
                    spec.quadrant,
                    spec.role,
                    spec.sequence_index,
                    f"{spec.center_z_mm:.12g}",
                    f"{spec.length_z_mm:.12g}",
                    f"{spec.magnetization_t[0]:.12g}",
                    f"{spec.magnetization_t[1]:.12g}",
                    f"{spec.magnetization_t[2]:.12g}",
                    json.dumps(spec.polygon_xy_mm),
                ]
            )


def make_manifest(config: ModelConfig, specs: Sequence[MagnetSpec]) -> dict[str, object]:
    return {
        "model": "WEP38-informed AQUA APPLE-X reference model",
        "mode": config.mode,
        "period_mm": config.period_mm,
        "periods": config.periods,
        "body_length_mm": config.periods * config.period_mm,
        "phase_fraction": config.phase_fraction,
        "phase_mm": config.phase_mm,
        "shifted_girders": ["G1", "G3"],
        "gap_setting_mm": config.gap_mm,
        "gap_datum_assumption": (
            "5.5 mm central square aperture at the 1.5 mm setting; "
            "each girder retracts radially by gap_mm-1.5"
        ),
        "central_hole_mm": config.central_hole_mm,
        "block_transverse_envelope_mm": [
            config.block_xy_mm,
            config.block_xy_mm,
        ],
        "block_chamfer_mm_assumed": config.chamfer_mm,
        "remanence_t": config.br_t,
        "blocks_per_period": 4,
        "end_fractions_body_outward": [0.75, 0.50, 0.25],
        "total_magnets": len(specs),
        "energy_gev": config.energy_gev,
        "limitations": [
            "locking tooth and holder omitted",
            "exact WEP38 polygon dimensions not published; symmetric chamfer assumed",
            "gap/mover datum inferred rather than copied from mechanical drawings",
            "ideal fixed magnetization; susceptibility and demagnetization omitted",
            "no measured block errors, sorting, shimming, force compensation, or trim coils",
            "end fractions implemented without device-specific fine optimization",
        ],
        "dry_run": config.dry_run,
    }


def first_harmonic(
    z_mm: np.ndarray,
    values: np.ndarray,
    period_mm: float,
    half_width_mm: float,
) -> dict[str, float]:
    mask = np.abs(z_mm) <= half_width_mm
    if np.count_nonzero(mask) < 8:
        raise ValueError("too few central samples for harmonic fit")
    z = z_mm[mask]
    y = values[mask]
    omega = 2.0 * math.pi / period_mm
    design = np.column_stack(
        [np.cos(omega * z), np.sin(omega * z), np.ones_like(z)]
    )
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    cosine, sine, offset = coefficients
    fitted = design @ coefficients
    return {
        "amplitude_t": float(math.hypot(float(cosine), float(sine))),
        "phase_rad": float(math.atan2(float(sine), float(cosine))),
        "offset_t": float(offset),
        "residual_rms_t": float(np.sqrt(np.mean((y - fitted) ** 2))),
    }


def trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def import_radia(radia_pythonpath: Path | None):
    if radia_pythonpath is not None:
        sys.path.insert(0, str(radia_pythonpath))
    try:
        import radia as rad  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError(
            "RADIA Python binding is not importable. Use --radia-pythonpath "
            "or install/build the ESRF RADIA Python binding."
        ) from exc
    return rad


def run_radia(
    config: ModelConfig,
    specs: Sequence[MagnetSpec],
    manifest: dict[str, object],
) -> dict[str, object]:
    rad = import_radia(config.radia_pythonpath)
    rad.UtiDelAll()

    objects: list[int] = []
    for spec in specs:
        objects.append(
            rad.ObjThckPgn(
                spec.center_z_mm,
                spec.length_z_mm,
                [list(point) for point in spec.polygon_xy_mm],
                "z",
                list(spec.magnetization_t),
            )
        )
    model = rad.ObjCnt(objects)

    z_min = min(s.center_z_mm - 0.5 * s.length_z_mm for s in specs)
    z_max = max(s.center_z_mm + 0.5 * s.length_z_mm for s in specs)
    z_mm = np.linspace(z_min - config.period_mm, z_max + config.period_mm, config.samples)

    field = np.empty((config.samples, 3), dtype=float)
    for index, z_value in enumerate(z_mm):
        field[index, :] = np.asarray(
            rad.Fld(model, "b", [0.0, 0.0, float(z_value)]), dtype=float
        )
    if not np.all(np.isfinite(field)):
        raise RuntimeError("RADIA returned non-finite field values")

    with (config.output_dir / "on_axis_field.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["z_mm", "Bx_T", "By_T", "Bz_T"])
        for z_value, vector in zip(z_mm, field):
            writer.writerow(
                [
                    f"{z_value:.12g}",
                    f"{vector[0]:.12g}",
                    f"{vector[1]:.12g}",
                    f"{vector[2]:.12g}",
                ]
            )

    fit_half_width = 0.5 * min(config.fit_periods, config.periods) * config.period_mm
    bx_fit = first_harmonic(z_mm, field[:, 0], config.period_mm, fit_half_width)
    by_fit = first_harmonic(z_mm, field[:, 1], config.period_mm, fit_half_width)

    period_cm = config.period_mm / 10.0
    k_bx_peak = 0.934 * bx_fit["amplitude_t"] * period_cm
    k_by_peak = 0.934 * by_fit["amplitude_t"] * period_cm
    k_rms = math.sqrt(0.5 * (k_bx_peak**2 + k_by_peak**2))

    z_m = z_mm / 1000.0
    first_integral_bx_tm = trapezoid(field[:, 0], z_m)
    first_integral_by_tm = trapezoid(field[:, 1], z_m)

    electron_rest_gev = 0.00051099895
    gamma = config.energy_gev / electron_rest_gev
    resonance_m = (config.period_mm * 1.0e-3) / (2.0 * gamma**2) * (
        1.0 + k_rms**2
    )

    analysis: dict[str, object] = {
        **manifest,
        "dry_run": False,
        "radia_model_handle": int(model),
        "sample_count": config.samples,
        "sample_range_mm": [float(z_mm[0]), float(z_mm[-1])],
        "fit_central_periods": min(config.fit_periods, config.periods),
        "first_harmonic": {"Bx": bx_fit, "By": by_fit},
        "K_peak_from_Bx": k_bx_peak,
        "K_peak_from_By": k_by_peak,
        "K_rms_WEP38_convention": k_rms,
        "first_integral_T_m": {
            "Bx": first_integral_bx_tm,
            "By": first_integral_by_tm,
        },
        "first_integral_G_m": {
            "Bx": first_integral_bx_tm * 1.0e4,
            "By": first_integral_by_tm * 1.0e4,
        },
        "fundamental_resonance_nm": resonance_m * 1.0e9,
    }

    if not config.no_plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib is not installed; field plot skipped", file=sys.stderr)
        else:
            figure, axis = plt.subplots(figsize=(10, 5))
            axis.plot(z_mm, field[:, 0], label="Bx")
            axis.plot(z_mm, field[:, 1], label="By")
            axis.set_xlabel("z [mm]")
            axis.set_ylabel("B [T]")
            axis.set_title(
                f"AQUA-style APPLE-X 18 mm — {config.mode}, "
                f"gap={config.gap_mm:g} mm"
            )
            axis.grid(True, alpha=0.3)
            axis.legend()
            figure.tight_layout()
            figure.savefig(config.output_dir / "on_axis_field.png", dpi=180)
            plt.close(figure)

    return analysis


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def execute(config: ModelConfig) -> dict[str, object]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    specs = build_specs(config)
    write_geometry_csv(config.output_dir / "magnet_geometry.csv", specs)
    manifest = make_manifest(config, specs)
    write_json(config.output_dir / "manifest.json", manifest)

    if config.dry_run:
        return manifest

    analysis = run_radia(config, specs, manifest)
    write_json(config.output_dir / "analysis.json", analysis)
    return analysis


def print_summary(config: ModelConfig) -> None:
    print("AQUA-style APPLE-X 18 mm reference model")
    print(f"  mode:             {config.mode}")
    print(f"  phase fraction:   {config.phase_fraction:g} lambda_u")
    print(f"  period:           {config.period_mm:g} mm")
    print(f"  periods:          {config.periods}")
    print(f"  gap setting:      {config.gap_mm:g} mm")
    print(f"  output directory: {config.output_dir}")
    print(f"  dry run:          {config.dry_run}")


def main(argv: Sequence[str] | None = None) -> int:
    try:
        config = parse_args(argv)
        print_summary(config)
        result = execute(config)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
