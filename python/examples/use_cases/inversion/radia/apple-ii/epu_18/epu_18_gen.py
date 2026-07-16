#!/usr/bin/env python3
"""WEP38-informed 18 mm AQUA APPLE-X RADIA reference model.

Corrected geometry interpretation:

* ``--gap-mm`` is the clear horizontal/vertical separation between adjacent
  girder faces (1.5 ... 4 mm), as labelled in WEP38 Figure 1.
* At minimum gap each inner face is at +/-0.75 mm, not +/-2.75 mm.
* The beam-facing 45-degree corner cut is derived from the published 5.5 mm
  central-aperture across-flats dimension, interpreted as the inscribed-circle / across-flats
  dimension of the rotated-square opening at minimum gap.
* The opposite outer corner uses the 5 mm chamfer shown in Figure 1.
* RADIA's standard NdFeB susceptibility is included by default.

The locking tooth, holder, exact production polygon, measured block errors,
and the slight unpublished end-magnet optimization are not available in the
paper. This is therefore a public-information magnetic reference model, not
manufacturer CAD or a construction release.

When ``--kickmap`` is enabled, the model calls RADIA's native
``FldFocKickPer`` periodic focusing-potential and second-order kick function on
a regular Cartesian mesh. The third returned matrix is exported as the RADIA
focusing potential. The primary Figure-7-style plot masks all points outside
the valid diamond aperture and the one-mesh-step boundary exclusion; a second,
explicitly diagnostic plot retains the full native Cartesian matrix. The kick
magnitude is plotted separately and only inside ``--kick-fit-r-max-mm`` so it
is not confused with the focusing potential or dominated by aperture-edge
values. The primary kick-map CSV contains only points safely inside the
rotated-square magnetic aperture. A separate diagnostic CSV retains the full
rectangular RADIA mesh with explicit region labels; points outside the diamond
or within one transverse mesh step of the magnet boundary are never used for
physical extrema or the reduced fit. The reduced cubic-potential fit uses all
valid Cartesian points inside ``--kick-fit-r-max-mm``. First-integral values are
deliberately not written to CSV or JSON outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np


MIN_GAP_MM = 1.5
MAX_GAP_MM = 4.0
SCRIPT_REVISION = "cartesian-potential-v9"
MODE_PHASES = {
    "linear-v": 0.0,
    "circular+": 0.25,
    "circular-": -0.25,
    "linear-h": 0.5,
}
MODE_RUN_NAMES = {
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
    central_aperture_across_flats_mm: float
    outer_chamfer_mm: float
    samples: int
    fit_periods: int
    energy_gev: float
    susceptibility: bool
    solve_precision: float
    solve_max_iterations: int
    strict_reference: bool
    kickmap: bool
    field_x_min_mm: float
    field_x_max_mm: float
    field_x_points: int
    field_y_min_mm: float
    field_y_max_mm: float
    field_y_points: int
    kickmap_max_harmonic: int
    kickmap_points_per_period: int
    kickmap_derivative_step_mm: float
    kick_fit_r_max_mm: float
    dry_run: bool
    no_plot: bool
    radia_pythonpath: Path | None

    @property
    def phase_fraction(self) -> float:
        return MODE_PHASES[self.mode]

    @property
    def phase_mm(self) -> float:
        return self.phase_fraction * self.period_mm

    @property
    def inner_chamfer_mm(self) -> float:
        # For a first-quadrant line x+y=gap+c, the distance to the origin is
        # (gap+c)/sqrt(2). Set twice that distance equal to the aperture.
        return self.central_aperture_across_flats_mm / math.sqrt(2.0) - MIN_GAP_MM

    @property
    def aperture_vertex_mm(self) -> float:
        """Axis intercept of the rotated-square magnetic aperture."""
        return self.gap_mm + self.inner_chamfer_mm


def env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    value = float(raw) if raw is not None else default
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    return int(raw) if raw is not None else default


def env_optional_float(name: str) -> float | None:
    """Read an optional finite floating-point environment value."""
    raw = os.environ.get(name)
    if raw is None:
        return None
    value = float(raw)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value



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
        description="Generate the corrected WEP38-informed AQUA APPLE-X model.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "mode", nargs="?", choices=tuple(MODE_PHASES), default="circular+"
    )
    parser.add_argument("output_dir", nargs="?", type=Path)
    parser.add_argument("--period-mm", type=float, default=env_float("PERIOD_MM", 18.0))
    parser.add_argument("--periods", type=int, default=env_int("PERIODS", 110))
    parser.add_argument(
        "--gap-mm",
        type=float,
        default=env_float("GAP_MM", MIN_GAP_MM),
        help="clear separation between adjacent x/y girder faces",
    )
    parser.add_argument("--br-t", type=float, default=env_float("BR_T", 1.35))
    parser.add_argument(
        "--block-xy-mm", type=float, default=env_float("BLOCK_XY_MM", 18.0)
    )
    parser.add_argument(
        "--central-aperture-across-flats-mm",
        type=float,
        default=env_float("CENTRAL_APERTURE_ACROSS_FLATS_MM", 5.5),
        help=("across-flats dimension of the rotated-square magnetic aperture "
              "at the 1.5 mm minimum gap"),
    )
    parser.add_argument(
        "--outer-chamfer-mm",
        type=float,
        default=env_float("OUTER_CHAMFER_MM", 5.0),
        help="outer-corner chamfer shown in WEP38 Figure 1",
    )
    parser.add_argument("--chamfer-mm", type=float, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--samples", type=int, default=env_int("SAMPLES", 4001))
    parser.add_argument("--fit-periods", type=int, default=env_int("FIT_PERIODS", 20))
    parser.add_argument("--energy-gev", type=float, default=env_float("ENERGY_GEV", 1.0))
    parser.add_argument(
        "--susceptibility",
        action=argparse.BooleanOptionalAction,
        default=env_bool("SUSCEPTIBILITY", True),
        help="apply RADIA MatStd('NdFeB', Br) and solve the interaction",
    )
    parser.add_argument(
        "--solve-precision", type=float, default=env_float("SOLVE_PRECISION", 1.0e-4)
    )
    parser.add_argument(
        "--solve-max-iterations",
        type=int,
        default=env_int("SOLVE_MAX_ITERATIONS", 1000),
    )
    parser.add_argument(
        "--strict-reference",
        action="store_true",
        default=env_bool("STRICT_REFERENCE", False),
        help="return status 2 when nominal WEP38 checks fail",
    )
    parser.add_argument(
        "--kickmap",
        action=argparse.BooleanOptionalAction,
        default=env_bool("KICKMAP", False),
        help="compute RADIA's native periodic second-order kick map",
    )
    parser.add_argument(
        "--field-x-min-mm",
        type=float,
        default=env_optional_float("FIELD_X_MIN_MM"),
        help=(
            "minimum Cartesian x coordinate; default is the negative axis "
            "vertex of the rotated-square magnetic aperture"
        ),
    )
    parser.add_argument(
        "--field-x-max-mm",
        type=float,
        default=env_optional_float("FIELD_X_MAX_MM"),
        help=(
            "maximum Cartesian x coordinate; default is the positive axis "
            "vertex of the rotated-square magnetic aperture"
        ),
    )
    parser.add_argument(
        "--field-x-points",
        type=int,
        default=env_int("FIELD_X_POINTS", 101),
        help="number of Cartesian x samples; must be odd so x=0 is included",
    )
    parser.add_argument(
        "--field-y-min-mm",
        type=float,
        default=env_optional_float("FIELD_Y_MIN_MM"),
        help=(
            "minimum Cartesian y coordinate; default is the negative axis "
            "vertex of the rotated-square magnetic aperture"
        ),
    )
    parser.add_argument(
        "--field-y-max-mm",
        type=float,
        default=env_optional_float("FIELD_Y_MAX_MM"),
        help=(
            "maximum Cartesian y coordinate; default is the positive axis "
            "vertex of the rotated-square magnetic aperture"
        ),
    )
    parser.add_argument(
        "--field-y-points",
        type=int,
        default=env_int("FIELD_Y_POINTS", 101),
        help="number of Cartesian y samples; must be odd so y=0 is included",
    )
    parser.add_argument(
        "--kickmap-max-harmonic",
        type=int,
        default=env_int("KICKMAP_MAX_HARMONIC", 5),
        help="maximum field harmonic passed to RADIA FldFocKickPer",
    )
    parser.add_argument(
        "--kickmap-points-per-period",
        type=int,
        default=env_int("KICKMAP_POINTS_PER_PERIOD", 32),
        help="longitudinal points per period passed to RADIA FldFocKickPer",
    )
    parser.add_argument(
        "--kickmap-derivative-step-mm",
        type=float,
        default=env_float("KICKMAP_DERIVATIVE_STEP_MM", 0.0),
        help=(
            "transverse derivative step passed to RADIA; zero lets RADIA use "
            "the native rectangular mesh spacing"
        ),
    )
    parser.add_argument(
        "--kick-fit-r-max-mm",
        type=float,
        default=env_float("KICK_FIT_R_MAX_MM", 0.8),
        help=(
            "maximum radius used for the optional cubic-potential / quadratic-kick fit; "
            "the native tabulated map still spans the configured Cartesian bounds"
        ),
    )
    parser.add_argument(
        "--dry-run", action="store_true", default=env_bool("DRY_RUN", False)
    )
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument(
        "--radia-pythonpath",
        type=Path,
        default=(
            Path(os.environ["RADIA_PYTHONPATH"])
            if os.environ.get("RADIA_PYTHONPATH")
            else None
        ),
    )
    return parser


def parse_args(argv: Sequence[str] | None = None) -> ModelConfig:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.chamfer_mm is not None:
        parser.error(
            "--chamfer-mm was ambiguous and is rejected. Use "
            "--outer-chamfer-mm 5.0; the inner chamfer is derived from the "
            "5.5 mm aperture."
        )

    script_dir = Path(__file__).resolve().parent
    output_dir = (
        args.output_dir
        if args.output_dir is not None
        else script_dir / "runs" / MODE_RUN_NAMES[args.mode]
    )
    inner_chamfer_mm = (
        float(args.central_aperture_across_flats_mm) / math.sqrt(2.0) - MIN_GAP_MM
    )
    aperture_vertex_mm = float(args.gap_mm) + inner_chamfer_mm
    field_x_min_mm = (
        float(args.field_x_min_mm)
        if args.field_x_min_mm is not None
        else -aperture_vertex_mm
    )
    field_x_max_mm = (
        float(args.field_x_max_mm)
        if args.field_x_max_mm is not None
        else aperture_vertex_mm
    )
    field_y_min_mm = (
        float(args.field_y_min_mm)
        if args.field_y_min_mm is not None
        else -aperture_vertex_mm
    )
    field_y_max_mm = (
        float(args.field_y_max_mm)
        if args.field_y_max_mm is not None
        else aperture_vertex_mm
    )

    config = ModelConfig(
        mode=args.mode,
        output_dir=output_dir.resolve(),
        period_mm=float(args.period_mm),
        periods=int(args.periods),
        gap_mm=float(args.gap_mm),
        br_t=float(args.br_t),
        block_xy_mm=float(args.block_xy_mm),
        central_aperture_across_flats_mm=float(args.central_aperture_across_flats_mm),
        outer_chamfer_mm=float(args.outer_chamfer_mm),
        samples=int(args.samples),
        fit_periods=int(args.fit_periods),
        energy_gev=float(args.energy_gev),
        susceptibility=bool(args.susceptibility),
        solve_precision=float(args.solve_precision),
        solve_max_iterations=int(args.solve_max_iterations),
        strict_reference=bool(args.strict_reference),
        kickmap=bool(args.kickmap),
        field_x_min_mm=field_x_min_mm,
        field_x_max_mm=field_x_max_mm,
        field_x_points=int(args.field_x_points),
        field_y_min_mm=field_y_min_mm,
        field_y_max_mm=field_y_max_mm,
        field_y_points=int(args.field_y_points),
        kickmap_max_harmonic=int(args.kickmap_max_harmonic),
        kickmap_points_per_period=int(args.kickmap_points_per_period),
        kickmap_derivative_step_mm=float(args.kickmap_derivative_step_mm),
        kick_fit_r_max_mm=float(args.kick_fit_r_max_mm),
        dry_run=bool(args.dry_run),
        no_plot=bool(args.no_plot),
        radia_pythonpath=(
            args.radia_pythonpath.resolve()
            if args.radia_pythonpath is not None
            else None
        ),
    )
    validate_config(config)
    return config


def validate_config(config: ModelConfig) -> None:
    positive = {
        "period_mm": config.period_mm,
        "br_t": config.br_t,
        "block_xy_mm": config.block_xy_mm,
        "central_aperture_across_flats_mm": config.central_aperture_across_flats_mm,
        "energy_gev": config.energy_gev,
        "solve_precision": config.solve_precision,
    }
    for name, value in positive.items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not MIN_GAP_MM <= config.gap_mm <= MAX_GAP_MM:
        raise ValueError(f"gap_mm must be within {MIN_GAP_MM} ... {MAX_GAP_MM} mm")
    if not 0 <= config.outer_chamfer_mm < config.block_xy_mm:
        raise ValueError("outer_chamfer_mm must be in [0, block_xy_mm)")
    if not 0 < config.inner_chamfer_mm < config.block_xy_mm:
        raise ValueError("derived inner chamfer is outside the block envelope")
    if config.periods < 1 or config.samples < 21 or config.fit_periods < 1:
        raise ValueError("periods >= 1, samples >= 21 and fit_periods >= 1 are required")
    if config.solve_max_iterations < 1:
        raise ValueError("solve_max_iterations must be at least 1")
    if config.kickmap:
        bounds = {
            "field_x_min_mm": config.field_x_min_mm,
            "field_x_max_mm": config.field_x_max_mm,
            "field_y_min_mm": config.field_y_min_mm,
            "field_y_max_mm": config.field_y_max_mm,
        }
        for name, value in bounds.items():
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if config.field_x_min_mm >= config.field_x_max_mm:
            raise ValueError("field_x_min_mm must be smaller than field_x_max_mm")
        if config.field_y_min_mm >= config.field_y_max_mm:
            raise ValueError("field_y_min_mm must be smaller than field_y_max_mm")
        if config.field_x_points < 3 or config.field_x_points % 2 == 0:
            raise ValueError("field_x_points must be odd and at least 3")
        if config.field_y_points < 3 or config.field_y_points % 2 == 0:
            raise ValueError("field_y_points must be odd and at least 3")
        centering_tolerance_mm = 1.0e-12
        if abs(config.field_x_min_mm + config.field_x_max_mm) > centering_tolerance_mm:
            raise ValueError("the Cartesian x range must be centred on x=0")
        if abs(config.field_y_min_mm + config.field_y_max_mm) > centering_tolerance_mm:
            raise ValueError("the Cartesian y range must be centred on y=0")
        if config.kickmap_max_harmonic < 1:
            raise ValueError("kickmap_max_harmonic must be at least 1")
        if config.kickmap_points_per_period < 4:
            raise ValueError("kickmap_points_per_period must be at least 4")
        if (
            not math.isfinite(config.kickmap_derivative_step_mm)
            or config.kickmap_derivative_step_mm < 0
        ):
            raise ValueError("kickmap_derivative_step_mm must be finite and non-negative")
        if not math.isfinite(config.kick_fit_r_max_mm) or config.kick_fit_r_max_mm <= 0:
            raise ValueError("kick_fit_r_max_mm must be finite and positive")
        available_half_width_mm = min(
            abs(config.field_x_min_mm),
            abs(config.field_x_max_mm),
            abs(config.field_y_min_mm),
            abs(config.field_y_max_mm),
        )
        if config.kick_fit_r_max_mm > available_half_width_mm:
            raise ValueError(
                "kick_fit_r_max_mm must fit inside the Cartesian x/y ranges"
            )



def file_sha256(path: Path) -> str:
    """Return a SHA-256 digest for output provenance."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _matrix_on_native_grid(
    raw: Any,
    y_points: int,
    x_points: int,
    name: str,
) -> np.ndarray:
    """Normalize a RADIA matrix to [physical-y index, physical-x index]."""
    matrix = np.asarray(raw, dtype=float)
    expected_shape = (y_points, x_points)
    transposed_shape = (x_points, y_points)
    if matrix.shape == expected_shape:
        return matrix
    if matrix.shape == transposed_shape and transposed_shape != expected_shape:
        return matrix.T
    if matrix.size == x_points * y_points:
        return matrix.reshape(expected_shape)
    raise RuntimeError(f"unexpected {name} matrix shape: {matrix.shape}")


def polygon_area(points: Sequence[tuple[float, float]]) -> float:
    pairs = list(zip(points, list(points[1:]) + [points[0]]))
    return 0.5 * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in pairs)


def quadrant_polygon(
    sx: int,
    sy: int,
    block_xy_mm: float,
    gap_mm: float,
    inner_chamfer_mm: float,
    outer_chamfer_mm: float,
) -> tuple[tuple[float, float], ...]:
    """Create the Figure-1-informed magnet envelope in one quadrant."""
    inner = gap_mm / 2.0
    outer = inner + block_xy_mm
    ci = inner_chamfer_mm
    co = outer_chamfer_mm
    first_quadrant = [
        (inner + ci, inner),
        (outer, inner),
        (outer, outer - co),
        (outer - co, outer),
        (inner, outer),
        (inner, inner + ci),
    ]
    result = [(sx * x, sy * y) for x, y in first_quadrant]
    if polygon_area(result) < 0:
        result.reverse()
    return tuple(result)


def halbach_direction(
    sequence_index: int,
    radial_unit: np.ndarray,
    row_sign: int,
    br_t: float,
) -> tuple[float, float, float]:
    ez = np.array([0.0, 0.0, 1.0])
    cycle = (radial_unit, ez, -radial_unit, -ez)
    vector = row_sign * br_t * cycle[sequence_index % 4]
    return float(vector[0]), float(vector[1]), float(vector[2])


def build_specs(config: ModelConfig) -> list[MagnetSpec]:
    standard_length = config.period_mm / 4.0
    body_blocks = config.periods * 4
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
            sx,
            sy,
            config.block_xy_mm,
            config.gap_mm,
            config.inner_chamfer_mm,
            config.outer_chamfer_mm,
        )
        z_shift = config.phase_mm if shifted else 0.0
        body_left = -0.5 * body_blocks * standard_length + z_shift

        for k in range(body_blocks):
            specs.append(
                MagnetSpec(
                    girder,
                    q_index,
                    "body",
                    k,
                    body_left + (k + 0.5) * standard_length,
                    standard_length,
                    halbach_direction(k, radial, row_sign, config.br_t),
                    polygon,
                )
            )

        # Figure 4: from the body outward = 3/4, 1/2, 1/4.
        end_fractions = (0.75, 0.50, 0.25)
        cursor = body_left
        for step, fraction in enumerate(end_fractions, start=1):
            length = fraction * standard_length
            sequence_index = -step
            specs.append(
                MagnetSpec(
                    girder,
                    q_index,
                    f"end-left-{fraction:g}",
                    sequence_index,
                    cursor - 0.5 * length,
                    length,
                    halbach_direction(sequence_index, radial, row_sign, config.br_t),
                    polygon,
                )
            )
            cursor -= length

        cursor = body_left + body_blocks * standard_length
        for step, fraction in enumerate(end_fractions):
            length = fraction * standard_length
            sequence_index = body_blocks + step
            specs.append(
                MagnetSpec(
                    girder,
                    q_index,
                    f"end-right-{fraction:g}",
                    sequence_index,
                    cursor + 0.5 * length,
                    length,
                    halbach_direction(sequence_index, radial, row_sign, config.br_t),
                    polygon,
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
    aperture_from_polygon = math.sqrt(2.0) * (
        MIN_GAP_MM + config.inner_chamfer_mm
    )
    script_path = Path(__file__).resolve()
    total_grid_points = config.field_x_points * config.field_y_points
    return {
        "model": "WEP38-informed AQUA APPLE-X, corrected gap geometry",
        "implementation": {
            "revision": SCRIPT_REVISION,
            "script_path": str(script_path),
            "script_sha256": file_sha256(script_path),
            "kick_fit_sampling": (
                "all Cartesian FldFocKickPer points inside kick_fit_r_max_mm"
            ),
        },
        "mode": config.mode,
        "period_mm": config.period_mm,
        "periods": config.periods,
        "body_length_mm": config.periods * config.period_mm,
        "phase_fraction": config.phase_fraction,
        "phase_mm": config.phase_mm,
        "shifted_girders": ["G1", "G3"],
        "inter_girder_face_gap_mm": config.gap_mm,
        "gap_definition": "inner faces at +/-gap_mm/2 in x and y",
        "central_aperture_across_flats_mm_at_min_gap": config.central_aperture_across_flats_mm,
        "central_aperture_across_flats_from_polygon_mm": aperture_from_polygon,
        "rotated_square_aperture_vertex_mm": config.aperture_vertex_mm,
        "rotated_square_aperture_condition": (
            "abs(x_mm)+abs(y_mm) <= rotated_square_aperture_vertex_mm"
        ),
        "derived_inner_chamfer_mm": config.inner_chamfer_mm,
        "outer_chamfer_mm": config.outer_chamfer_mm,
        "block_transverse_envelope_mm": [config.block_xy_mm, config.block_xy_mm],
        "remanence_t": config.br_t,
        "blocks_per_period": 4,
        "end_fractions_body_outward": [0.75, 0.50, 0.25],
        "total_magnets": len(specs),
        "energy_gev": config.energy_gev,
        "ndfeb_susceptibility_enabled": config.susceptibility,
        "radia_solve_precision": config.solve_precision,
        "radia_solve_max_iterations": config.solve_max_iterations,
        "kickmap_requested": config.kickmap,
        "shared_cartesian_field_kickmap_grid": (
            {
                "coordinate_system": "regular Cartesian RADIA mesh",
                "x_min_mm": config.field_x_min_mm,
                "x_max_mm": config.field_x_max_mm,
                "x_points": config.field_x_points,
                "y_min_mm": config.field_y_min_mm,
                "y_max_mm": config.field_y_max_mm,
                "y_points": config.field_y_points,
                "total_points": total_grid_points,
                "radia_function": "FldFocKickPer",
                "physical_aperture_shape": "rotated square (diamond in x-y)",
                "physical_aperture_vertex_mm": config.aperture_vertex_mm,
                "kick_fit_r_max_mm": config.kick_fit_r_max_mm,
            }
            if config.kickmap
            else None
        ),
        "particle": {
            "species_assumed_by_RADIA_FldFocKickPer": "electron",
            "total_energy_gev": config.energy_gev,
        },
        "limitations": [
            "locking tooth, holder and pipe support omitted",
            "exact production polygon is unpublished",
            "unpublished slight end-magnet optimization is not reproduced",
            "no measured block errors, sorting, shimming or correction coils",
            "demagnetization margin and mechanical tolerances are not evaluated",
            "the full rectangular RADIA mesh is exported only as a diagnostic file with explicit region labels",
            "the primary kick map excludes points outside the diamond and one transverse mesh step near its boundary",
            "FldFocKickPer evaluates the periodic body; finite-device end fields are diagnosed separately",
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
        raise ValueError("too few samples for the central harmonic fit")
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


def cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    result = np.zeros_like(y, dtype=float)
    result[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return result


def import_radia(path: Path | None):
    if path is not None:
        sys.path.insert(0, str(path))
    try:
        import radia as rad  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError(
            "ESRF RADIA is not importable; use --radia-pythonpath"
        ) from exc
    return rad


def reshape_field(raw: Any, samples: int) -> np.ndarray:
    field = np.asarray(raw, dtype=float)
    if field.shape == (samples, 3):
        return field
    if field.size == samples * 3:
        return field.reshape(samples, 3)
    raise RuntimeError(f"unexpected RADIA field shape: {field.shape}")


def reference_checks(
    config: ModelConfig,
    k_rms: float,
    kx: float,
    ky: float,
    i2_gm2: dict[str, float],
) -> dict[str, object]:
    nominal = (
        config.mode in {"circular+", "circular-"}
        and math.isclose(config.period_mm, 18.0, abs_tol=1e-9)
        and math.isclose(config.gap_mm, 1.5, abs_tol=1e-9)
        and math.isclose(config.br_t, 1.35, abs_tol=1e-9)
        and math.isclose(config.block_xy_mm, 18.0, abs_tol=1e-9)
        and math.isclose(config.central_aperture_across_flats_mm, 5.5, abs_tol=1e-9)
    )
    if not nominal:
        return {"applicable": False}
    mismatch = abs(kx - ky) / max(abs(kx), abs(ky), 1e-15)
    checks = {
        "applicable": True,
        "K_rms_expected_band": [1.05, 1.30],
        "K_rms_pass": 1.05 <= k_rms <= 1.30,
        "circular_component_relative_mismatch": mismatch,
        "circular_symmetry_pass": mismatch <= 0.02,
        "second_integral_limit_G_m2": 0.05,
        "second_integral_pass": all(abs(v) < 0.05 for v in i2_gm2.values()),
        "note": "integral limits require the unpublished final end optimization",
    }
    checks["all_pass"] = all(
        checks[key]
        for key in (
            "K_rms_pass",
            "circular_symmetry_pass",
            "second_integral_pass",
        )
    )
    return checks






def fit_symplectic_second_order(
    x_m: np.ndarray,
    y_m: np.ndarray,
    theta_x_rad: np.ndarray,
    theta_y_rad: np.ndarray,
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    """Fit quadratic kicks as derivatives of one cubic generating potential.

    W = c10*x + c01*y + c20*x^2 + c11*x*y + c02*y^2
        + c30*x^3 + c21*x^2*y + c12*x*y^2 + c03*y^3

    theta_x = -dW/dx and theta_y = -dW/dy. Coordinates are normalized
    internally before least-squares fitting to avoid short-aperture conditioning
    problems; the reported coefficients are converted back to SI coordinates.
    """
    x_scale = max(float(np.max(np.abs(x_m))), 1.0e-12)
    y_scale = max(float(np.max(np.abs(y_m))), 1.0e-12)
    x_normalized = x_m / x_scale
    y_normalized = y_m / y_scale

    rows: list[list[float]] = []
    targets: list[float] = []
    for x_value, y_value, tx_value, ty_value in zip(
        x_normalized, y_normalized, theta_x_rad, theta_y_rad
    ):
        rows.append(
            [
                -1.0,
                0.0,
                -2.0 * x_value,
                -y_value,
                0.0,
                -3.0 * x_value * x_value,
                -2.0 * x_value * y_value,
                -y_value * y_value,
                0.0,
            ]
        )
        targets.append(x_scale * float(tx_value))
        rows.append(
            [
                0.0,
                -1.0,
                0.0,
                -x_value,
                -2.0 * y_value,
                0.0,
                -x_value * x_value,
                -2.0 * x_value * y_value,
                -3.0 * y_value * y_value,
            ]
        )
        targets.append(y_scale * float(ty_value))

    design = np.asarray(rows, dtype=float)
    target = np.asarray(targets, dtype=float)
    normalized_coefficients, _, rank, singular_values = np.linalg.lstsq(
        design, target, rcond=None
    )
    if rank < 9:
        raise RuntimeError(
            f"second-order symplectic fit is rank deficient: rank={rank}, expected 9"
        )

    a10, a01, a20, a11, a02, a30, a21, a12, a03 = normalized_coefficients
    physical_coefficients = {
        "c10": float(a10 / x_scale),
        "c01": float(a01 / y_scale),
        "c20_per_m": float(a20 / (x_scale * x_scale)),
        "c11_per_m": float(a11 / (x_scale * y_scale)),
        "c02_per_m": float(a02 / (y_scale * y_scale)),
        "c30_per_m2": float(a30 / (x_scale**3)),
        "c21_per_m2": float(a21 / (x_scale * x_scale * y_scale)),
        "c12_per_m2": float(a12 / (x_scale * y_scale * y_scale)),
        "c03_per_m2": float(a03 / (y_scale**3)),
    }

    x_value = x_normalized
    y_value = y_normalized
    fitted_x = -(
        a10
        + 2.0 * a20 * x_value
        + a11 * y_value
        + 3.0 * a30 * x_value * x_value
        + 2.0 * a21 * x_value * y_value
        + a12 * y_value * y_value
    ) / x_scale
    fitted_y = -(
        a01
        + a11 * x_value
        + 2.0 * a02 * y_value
        + a21 * x_value * x_value
        + 2.0 * a12 * x_value * y_value
        + 3.0 * a03 * y_value * y_value
    ) / y_scale

    residual_x = theta_x_rad - fitted_x
    residual_y = theta_y_rad - fitted_y
    residual_vector = np.concatenate([residual_x, residual_y])
    fit_record: dict[str, object] = {
        "model": "cubic generating potential; quadratic symplectic kicks",
        "coefficient_definition": (
            "W=c10*x+c01*y+c20*x^2+c11*x*y+c02*y^2+"
            "c30*x^3+c21*x^2*y+c12*x*y^2+c03*y^3; "
            "theta_x=-dW/dx, theta_y=-dW/dy; x,y in metres"
        ),
        "coefficients_SI": physical_coefficients,
        "normalization_scales_m": {"x": x_scale, "y": y_scale},
        "matrix_rank": int(rank),
        "singular_values": [float(value) for value in singular_values],
        "residual_rms_rad": float(np.sqrt(np.mean(residual_vector**2))),
        "residual_max_abs_rad": float(np.max(np.abs(residual_vector))),
        "residual_rms_urad": float(np.sqrt(np.mean(residual_vector**2)) * 1.0e6),
        "residual_max_abs_urad": float(np.max(np.abs(residual_vector)) * 1.0e6),
    }
    return fit_record, fitted_x, fitted_y



def evaluate_symplectic_second_order(
    coefficients: dict[str, float],
    x_m: np.ndarray,
    y_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the quadratic kicks derived from reported SI coefficients."""
    c10 = coefficients["c10"]
    c01 = coefficients["c01"]
    c20 = coefficients["c20_per_m"]
    c11 = coefficients["c11_per_m"]
    c02 = coefficients["c02_per_m"]
    c30 = coefficients["c30_per_m2"]
    c21 = coefficients["c21_per_m2"]
    c12 = coefficients["c12_per_m2"]
    c03 = coefficients["c03_per_m2"]
    fitted_x = -(
        c10 + 2.0 * c20 * x_m + c11 * y_m
        + 3.0 * c30 * x_m * x_m + 2.0 * c21 * x_m * y_m
        + c12 * y_m * y_m
    )
    fitted_y = -(
        c01 + c11 * x_m + 2.0 * c02 * y_m
        + c21 * x_m * x_m + 2.0 * c12 * x_m * y_m
        + 3.0 * c03 * y_m * y_m
    )
    return fitted_x, fitted_y

def locate_origin_index(
    r_mm: np.ndarray,
    *,
    tolerance_mm: float = 1.0e-9,
) -> int:
    """Return the retained-grid index corresponding to the magnetic axis.

    RADIA may return the nominal zero coordinate with a very small floating-point
    offset. Selecting the minimum-radius point and then applying a tolerance is
    therefore more reliable than testing coordinate values for exact equality.
    """
    radii = np.asarray(r_mm, dtype=float).reshape(-1)
    if radii.size == 0:
        raise RuntimeError("retained kick-map grid is empty")
    if not np.all(np.isfinite(radii)):
        raise RuntimeError("retained kick-map radii contain non-finite values")
    if tolerance_mm <= 0.0 or not math.isfinite(tolerance_mm):
        raise ValueError("origin tolerance must be finite and positive")

    origin_index = int(np.argmin(radii))
    minimum_radius_mm = float(radii[origin_index])
    if minimum_radius_mm > tolerance_mm:
        raise RuntimeError(
            "Retained kick-map grid does not contain the origin: "
            f"minimum radius is {minimum_radius_mm:.12g} mm"
        )
    return origin_index


def compute_second_order_kickmap(
    rad: Any,
    model: int,
    config: ModelConfig,
) -> dict[str, object]:
    """Compute the periodic second-order map on one Cartesian RADIA mesh.

    RADIA necessarily evaluates a full rectangular mesh. The primary exported
    map and all physical extrema use only a conservative interior diamond:

        abs(x) + abs(y) <= aperture_vertex_mm - boundary_margin_mm

    where ``boundary_margin_mm`` is one maximum transverse mesh step (or the
    explicit RADIA derivative step, when larger). The complete rectangular
    result is retained only in a clearly labelled diagnostic CSV.
    """
    x_points = config.field_x_points
    y_points = config.field_y_points
    total_points = x_points * y_points
    x_range_mm = config.field_x_max_mm - config.field_x_min_mm
    y_range_mm = config.field_y_max_mm - config.field_y_min_mm
    period_start_mm = -0.5 * config.period_mm
    precision = [
        config.kickmap_max_harmonic,
        config.kickmap_points_per_period,
        config.kickmap_derivative_step_mm,
        config.kickmap_derivative_step_mm,
    ]
    try:
        result = rad.FldFocKickPer(
            model,
            [0.0, 0.0, period_start_mm],
            [0.0, 0.0, 1.0],
            config.period_mm,
            config.periods,
            [1.0, 0.0, 0.0],
            x_range_mm,
            x_points,
            y_range_mm,
            y_points,
            "AQUA APPLE-X 18 mm periodic second-order Cartesian kick map",
            precision,
            "rad",
            config.energy_gev,
            # Compatibility workaround: the tested macOS Python binding
            # corrupts the heap with the formatted-text option "tab".
            "fix",
        )
    except (AttributeError, RuntimeError, TypeError) as exc:
        raise RuntimeError(
            "RADIA FldFocKickPer failed; verify that the installed ESRF RADIA "
            "Python binding includes the periodic focusing-kick function"
        ) from exc
    if not isinstance(result, (list, tuple)) or len(result) < 6:
        raise RuntimeError("unexpected RADIA FldFocKickPer return value")

    kick_first = _matrix_on_native_grid(
        result[0], y_points, x_points, "first kick"
    )
    kick_second = _matrix_on_native_grid(
        result[1], y_points, x_points, "second kick"
    )
    integrated_btr_squared = _matrix_on_native_grid(
        result[2], y_points, x_points, "integrated squared transverse field"
    )
    first_axis_mm = np.asarray(result[3], dtype=float).reshape(-1)
    second_axis_mm = np.asarray(result[4], dtype=float).reshape(-1)
    if first_axis_mm.size != x_points or second_axis_mm.size != y_points:
        raise RuntimeError("RADIA returned inconsistent transverse axes")

    # With n1=+x and ns=+z, RADIA's second transverse basis is n1 x ns=-y.
    x_axis_mm = first_axis_mm
    y_axis_mm = -second_axis_mm
    kick_x_grid_rad = kick_first
    kick_y_grid_rad = -kick_second

    x_order = np.argsort(x_axis_mm)
    y_order = np.argsort(y_axis_mm)
    x_axis_mm = x_axis_mm[x_order]
    y_axis_mm = y_axis_mm[y_order]
    endpoint_tolerance_mm = 1.0e-9
    if not (
        math.isclose(
            float(x_axis_mm[0]),
            config.field_x_min_mm,
            abs_tol=endpoint_tolerance_mm,
        )
        and math.isclose(
            float(x_axis_mm[-1]),
            config.field_x_max_mm,
            abs_tol=endpoint_tolerance_mm,
        )
        and math.isclose(
            float(y_axis_mm[0]),
            config.field_y_min_mm,
            abs_tol=endpoint_tolerance_mm,
        )
        and math.isclose(
            float(y_axis_mm[-1]),
            config.field_y_max_mm,
            abs_tol=endpoint_tolerance_mm,
        )
    ):
        raise RuntimeError(
            "RADIA returned transverse axes inconsistent with the requested "
            "centred Cartesian ranges"
        )
    kick_x_grid_rad = kick_x_grid_rad[np.ix_(y_order, x_order)]
    kick_y_grid_rad = kick_y_grid_rad[np.ix_(y_order, x_order)]
    integrated_btr_squared = integrated_btr_squared[np.ix_(y_order, x_order)]

    if not (
        np.all(np.isfinite(kick_x_grid_rad))
        and np.all(np.isfinite(kick_y_grid_rad))
        and np.all(np.isfinite(integrated_btr_squared))
    ):
        raise RuntimeError("RADIA FldFocKickPer returned non-finite values")

    xx_mm, yy_mm = np.meshgrid(x_axis_mm, y_axis_mm, indexing="xy")
    rows_grid, columns_grid = np.indices((y_points, x_points))
    x_mm = xx_mm.reshape(-1)
    y_mm = yy_mm.reshape(-1)
    rows = rows_grid.reshape(-1)
    columns = columns_grid.reshape(-1)
    r_mm = np.hypot(x_mm, y_mm)
    phi_rad = np.arctan2(y_mm, x_mm)
    theta_x_rad = kick_x_grid_rad.reshape(-1)
    theta_y_rad = kick_y_grid_rad.reshape(-1)
    btr_squared_native = integrated_btr_squared.reshape(-1)

    dx_mm = float(np.max(np.abs(np.diff(x_axis_mm))))
    dy_mm = float(np.max(np.abs(np.diff(y_axis_mm))))
    boundary_margin_mm = max(
        dx_mm,
        dy_mm,
        float(config.kickmap_derivative_step_mm),
    )
    geometric_aperture_mask_grid = (
        np.abs(xx_mm) + np.abs(yy_mm)
        <= config.aperture_vertex_mm + 1.0e-12
    )
    valid_aperture_limit_mm = config.aperture_vertex_mm - boundary_margin_mm
    if valid_aperture_limit_mm <= 0.0:
        raise RuntimeError(
            "Cartesian mesh is too coarse for the rotated-square aperture: "
            f"automatic boundary margin {boundary_margin_mm:.6g} mm leaves no "
            "valid interior"
        )
    valid_aperture_mask_grid = (
        np.abs(xx_mm) + np.abs(yy_mm)
        <= valid_aperture_limit_mm + 1.0e-12
    )
    boundary_exclusion_mask_grid = (
        geometric_aperture_mask_grid & ~valid_aperture_mask_grid
    )
    outside_aperture_mask_grid = ~geometric_aperture_mask_grid

    geometric_aperture_mask = geometric_aperture_mask_grid.reshape(-1)
    valid_aperture_mask = valid_aperture_mask_grid.reshape(-1)
    boundary_exclusion_mask = boundary_exclusion_mask_grid.reshape(-1)
    outside_aperture_mask = outside_aperture_mask_grid.reshape(-1)

    geometric_aperture_point_count = int(
        np.count_nonzero(geometric_aperture_mask)
    )
    valid_aperture_point_count = int(np.count_nonzero(valid_aperture_mask))
    boundary_exclusion_point_count = int(
        np.count_nonzero(boundary_exclusion_mask)
    )
    outside_aperture_point_count = int(np.count_nonzero(outside_aperture_mask))
    if valid_aperture_point_count == 0:
        raise RuntimeError("Cartesian grid contains no valid physical-aperture points")

    field_points = [[float(x), float(y), 0.0] for x, y in zip(x_mm, y_mm)]
    centre_field = reshape_field(rad.Fld(model, "b", field_points), total_points)
    if not np.all(np.isfinite(centre_field)):
        raise RuntimeError("RADIA returned non-finite field-grid values")

    fit_mask_grid = (
        xx_mm * xx_mm + yy_mm * yy_mm
        <= config.kick_fit_r_max_mm * config.kick_fit_r_max_mm + 1.0e-12
    ) & valid_aperture_mask_grid
    fit_mask = fit_mask_grid.reshape(-1)
    fit_point_count = int(np.count_nonzero(fit_mask))
    if fit_point_count < 9:
        raise RuntimeError(
            "too few valid Cartesian points inside kick_fit_r_max_mm for a "
            "full cubic-potential fit"
        )

    fit_record, fitted_fit_x_rad, fitted_fit_y_rad = fit_symplectic_second_order(
        x_mm[fit_mask] / 1000.0,
        y_mm[fit_mask] / 1000.0,
        theta_x_rad[fit_mask],
        theta_y_rad[fit_mask],
    )
    fit_record["fit_radius_mm"] = config.kick_fit_r_max_mm
    fit_record["fit_point_count"] = fit_point_count
    fit_record["sampling_grid"] = (
        "all valid regular Cartesian FldFocKickPer points inside fit_radius_mm"
    )
    fit_record["cartesian_grid_shape"] = [y_points, x_points]
    fit_record["evaluation_domain_note"] = (
        "coefficients are fitted only inside fit_radius_mm; fitted values "
        "outside that disk are extrapolations"
    )

    fitted_x_rad, fitted_y_rad = evaluate_symplectic_second_order(
        fit_record["coefficients_SI"],
        x_mm / 1000.0,
        y_mm / 1000.0,
    )
    residual_x_rad = theta_x_rad - fitted_x_rad
    residual_y_rad = theta_y_rad - fitted_y_rad
    measured_magnitude_urad = np.hypot(theta_x_rad, theta_y_rad) * 1.0e6
    fitted_magnitude_urad = np.hypot(fitted_x_rad, fitted_y_rad) * 1.0e6

    def point_region(index: int) -> str:
        if fit_mask[index]:
            return "fit_domain"
        if valid_aperture_mask[index]:
            return "valid_physical_aperture"
        if boundary_exclusion_mask[index]:
            return "boundary_exclusion"
        return "outside_geometric_aperture"

    common_header = [
        "point_index",
        "grid_row_y",
        "grid_column_x",
        "point_region",
        "inside_geometric_aperture",
        "inside_valid_physical_aperture",
        "inside_fit_domain",
        "x_mm",
        "y_mm",
        "r_mm",
        "phi_rad",
        "Bx_at_z0_T",
        "By_at_z0_T",
        "Bz_at_z0_T",
        "focusing_potential_1e2_T2_mm2",
        "kick_x_rad",
        "kick_y_rad",
        "kick_x_urad",
        "kick_y_urad",
        "fitted_kick_x_rad",
        "fitted_kick_y_rad",
        "fitted_kick_x_urad",
        "fitted_kick_y_urad",
        "residual_kick_x_rad",
        "residual_kick_y_rad",
        "residual_kick_x_urad",
        "residual_kick_y_urad",
        "kick_magnitude_urad",
        "fitted_kick_magnitude_urad",
    ]

    def map_row(index: int) -> list[object]:
        return [
            index,
            int(rows[index]),
            int(columns[index]),
            point_region(index),
            int(bool(geometric_aperture_mask[index])),
            int(bool(valid_aperture_mask[index])),
            int(bool(fit_mask[index])),
            f"{x_mm[index]:.12g}",
            f"{y_mm[index]:.12g}",
            f"{r_mm[index]:.12g}",
            f"{phi_rad[index]:.12g}",
            f"{centre_field[index, 0]:.12g}",
            f"{centre_field[index, 1]:.12g}",
            f"{centre_field[index, 2]:.12g}",
            f"{btr_squared_native[index]:.12g}",
            f"{theta_x_rad[index]:.12g}",
            f"{theta_y_rad[index]:.12g}",
            f"{theta_x_rad[index] * 1.0e6:.12g}",
            f"{theta_y_rad[index] * 1.0e6:.12g}",
            f"{fitted_x_rad[index]:.12g}",
            f"{fitted_y_rad[index]:.12g}",
            f"{fitted_x_rad[index] * 1.0e6:.12g}",
            f"{fitted_y_rad[index] * 1.0e6:.12g}",
            f"{residual_x_rad[index]:.12g}",
            f"{residual_y_rad[index]:.12g}",
            f"{residual_x_rad[index] * 1.0e6:.12g}",
            f"{residual_y_rad[index] * 1.0e6:.12g}",
            f"{measured_magnitude_urad[index]:.12g}",
            f"{fitted_magnitude_urad[index]:.12g}",
        ]

    # Primary physical map: outside and boundary-exclusion points are absent.
    output_path = config.output_dir / "kickmap_integrated.csv"
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(common_header)
        for index in np.flatnonzero(valid_aperture_mask):
            writer.writerow(map_row(int(index)))

    # Full rectangular RADIA result for debugging and reproducibility only.
    diagnostic_output_path = (
        config.output_dir / "kickmap_cartesian_diagnostic.csv"
    )
    with diagnostic_output_path.open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(common_header)
        for index in range(total_points):
            writer.writerow(map_row(index))

    # Dedicated full-mesh focusing-potential export. This is the third matrix
    # returned by FldFocKickPer and is the quantity comparable to WEP38 Fig. 7.
    focusing_potential_output_path = (
        config.output_dir / "focusing_potential_cartesian.csv"
    )
    with focusing_potential_output_path.open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "grid_row_y",
                "grid_column_x",
                "point_region",
                "inside_geometric_aperture",
                "inside_valid_physical_aperture",
                "x_mm",
                "y_mm",
                "focusing_potential_1e2_T2_mm2",
            ]
        )
        for index in range(total_points):
            writer.writerow(
                [
                    int(rows[index]),
                    int(columns[index]),
                    point_region(index),
                    int(bool(geometric_aperture_mask[index])),
                    int(bool(valid_aperture_mask[index])),
                    f"{x_mm[index]:.12g}",
                    f"{y_mm[index]:.12g}",
                    f"{btr_squared_native[index]:.12g}",
                ]
            )

    fit_output_path = config.output_dir / "kickmap_fit_cartesian.csv"
    fit_rows, fit_columns = np.nonzero(fit_mask_grid)
    fit_residual_x_rad = theta_x_rad[fit_mask] - fitted_fit_x_rad
    fit_residual_y_rad = theta_y_rad[fit_mask] - fitted_fit_y_rad
    with fit_output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "grid_row_y",
                "grid_column_x",
                "x_mm",
                "y_mm",
                "r_mm",
                "kick_x_rad",
                "kick_y_rad",
                "fitted_kick_x_rad",
                "fitted_kick_y_rad",
                "residual_kick_x_urad",
                "residual_kick_y_urad",
            ]
        )
        fit_x_mm = x_mm[fit_mask]
        fit_y_mm = y_mm[fit_mask]
        fit_theta_x_rad = theta_x_rad[fit_mask]
        fit_theta_y_rad = theta_y_rad[fit_mask]
        for index in range(fit_point_count):
            x_value = float(fit_x_mm[index])
            y_value = float(fit_y_mm[index])
            writer.writerow(
                [
                    int(fit_rows[index]),
                    int(fit_columns[index]),
                    f"{x_value:.12g}",
                    f"{y_value:.12g}",
                    f"{math.hypot(x_value, y_value):.12g}",
                    f"{fit_theta_x_rad[index]:.12g}",
                    f"{fit_theta_y_rad[index]:.12g}",
                    f"{fitted_fit_x_rad[index]:.12g}",
                    f"{fitted_fit_y_rad[index]:.12g}",
                    f"{fit_residual_x_rad[index] * 1.0e6:.12g}",
                    f"{fit_residual_y_rad[index] * 1.0e6:.12g}",
                ]
            )

    # Avoid confusing stale output from pre-v6 runs with current products.
    for stale_name in ("kickmap_3d.png",):
        stale_path = config.output_dir / stale_name
        try:
            stale_path.unlink()
        except FileNotFoundError:
            pass

    focusing_potential_plot_status = "not requested"
    focusing_potential_diagnostic_plot_status = "not requested"
    kick_fit_plot_status = "not requested"
    if not config.no_plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print(
                "Warning: matplotlib unavailable; focusing-potential and "
                "kick-map plots skipped",
                file=sys.stderr,
            )
            focusing_potential_plot_status = "matplotlib unavailable"
            focusing_potential_diagnostic_plot_status = "matplotlib unavailable"
            kick_fit_plot_status = "matplotlib unavailable"
        else:
            # WEP38 Figure 7 is a focusing-potential plot, not a plot of the
            # vector kick magnitude. The primary comparison product uses only
            # the valid beam-accessible diamond, excluding the numerical
            # boundary band. Masking rather than clipping prevents points in or
            # behind the magnets from generating artificial vertical walls.
            #
            # The third matrix returned by FldFocKickPer is displayed without
            # numerical rescaling. One plotted unit is labelled as
            # 10^2 T^2 mm^2, which is equivalent to 10^-4 T^2 m^2. This keeps
            # the Figure-7-style numerical scale near 1 ... 3 while expressing
            # the length unit consistently with the x/y axes in millimetres.
            focusing_potential_plot_scale_T2_mm2 = 1.0e2
            integrated_btr_squared_plot_units = np.asarray(
                integrated_btr_squared,
                dtype=float,
            )

            potential_for_plot = np.ma.masked_where(
                ~valid_aperture_mask_grid,
                integrated_btr_squared_plot_units,
            )
            valid_potential_values = integrated_btr_squared_plot_units[
                valid_aperture_mask_grid
            ]
            valid_potential_min = float(np.min(valid_potential_values))
            valid_potential_max = float(np.max(valid_potential_values))

            potential_figure = plt.figure(figsize=(9, 6))
            potential_axis = potential_figure.add_subplot(111, projection="3d")
            potential_surface = potential_axis.plot_surface(
                xx_mm,
                yy_mm,
                potential_for_plot,
                rstride=max(1, y_points // 80),
                cstride=max(1, x_points // 80),
                linewidth=0.2,
                antialiased=True,
                cmap="plasma",
            )
            potential_axis.set_xlim(float(x_axis_mm[0]), float(x_axis_mm[-1]))
            potential_axis.set_ylim(float(y_axis_mm[0]), float(y_axis_mm[-1]))
            potential_axis.set_zlim(valid_potential_min, valid_potential_max)
            potential_axis.set_xlabel("x [mm]")
            potential_axis.set_ylabel("y [mm]")
            potential_axis.set_zlabel(
                r"focusing potential [$10^{2} \mathrm{T^2 \, mm^2}$]"
            )
            potential_axis.set_title(
                f"RADIA focusing potential — {config.mode}, "
                f"gap={config.gap_mm:g} mm"
            )
            potential_axis.view_init(elev=27, azim=-58)
            potential_colorbar = potential_figure.colorbar(
                potential_surface,
                ax=potential_axis,
                shrink=0.72,
                pad=0.08,
                aspect=20,
            )
            potential_colorbar.set_label(
                r"focusing potential [$10^{2} \mathrm{T^2 \, mm^2}$]"
            )
            potential_figure.tight_layout()
            potential_figure.savefig(
                config.output_dir / "focusing_potential_3d.png",
                dpi=180,
            )
            plt.close(potential_figure)
            focusing_potential_plot_status = "written"

            # Preserve the full native matrix only as an explicitly diagnostic
            # plot. It can contain values outside the magnetic opening and must
            # not be interpreted as a beam-accessible focusing-potential map.
            diagnostic_potential_figure = plt.figure(figsize=(9, 6))
            diagnostic_potential_axis = diagnostic_potential_figure.add_subplot(
                111, projection="3d"
            )
            diagnostic_surface = diagnostic_potential_axis.plot_surface(
                xx_mm,
                yy_mm,
                integrated_btr_squared_plot_units,
                rstride=max(1, y_points // 80),
                cstride=max(1, x_points // 80),
                linewidth=0.2,
                antialiased=True,
                cmap="plasma",
            )
            diagnostic_potential_axis.set_xlabel("x [mm]")
            diagnostic_potential_axis.set_ylabel("y [mm]")
            diagnostic_potential_axis.set_zlabel(
                r"focusing potential [$10^{2} \mathrm{T^2 \, mm^2}$]"
            )
            diagnostic_potential_axis.set_title(
                f"Diagnostic full Cartesian focusing potential — "
                f"{config.mode}, gap={config.gap_mm:g} mm"
            )
            diagnostic_potential_axis.view_init(elev=27, azim=-58)
            diagnostic_colorbar = diagnostic_potential_figure.colorbar(
                diagnostic_surface,
                ax=diagnostic_potential_axis,
                shrink=0.72,
                pad=0.08,
                aspect=20,
            )
            diagnostic_colorbar.set_label(
                r"focusing potential [$10^{2} \mathrm{T^2 \, mm^2}$]"
            )
            diagnostic_potential_figure.tight_layout()
            diagnostic_potential_figure.savefig(
                config.output_dir
                / "focusing_potential_full_cartesian_diagnostic_3d.png",
                dpi=180,
            )
            plt.close(diagnostic_potential_figure)
            focusing_potential_diagnostic_plot_status = "written"

            # Keep the actual second-order kick visualization separate and
            # restrict it to the declared fit disk. This avoids aperture-edge
            # walls and makes clear that this is not the Figure-7 quantity.
            kick_figure = plt.figure(figsize=(9, 6))
            kick_axis = kick_figure.add_subplot(111, projection="3d")
            kick_surface = kick_axis.plot_trisurf(
                x_mm[fit_mask],
                y_mm[fit_mask],
                measured_magnitude_urad[fit_mask],
                linewidth=0.2,
                antialiased=True,
                cmap="viridis",
            )
            kick_axis.scatter(
                x_mm[fit_mask],
                y_mm[fit_mask],
                fitted_magnitude_urad[fit_mask],
                s=10,
                c="black",
                alpha=0.7,
                label="fitted values",
            )
            kick_axis.set_xlabel("x [mm]")
            kick_axis.set_ylabel("y [mm]")
            kick_axis.set_zlabel("|kick| [urad]")
            kick_axis.set_title(
                f"RADIA second-order kick magnitude — r <= "
                f"{config.kick_fit_r_max_mm:g} mm"
            )
            kick_colorbar = kick_figure.colorbar(
                kick_surface,
                ax=kick_axis,
                shrink=0.72,
                pad=0.08,
                aspect=20,
            )
            kick_colorbar.set_label("|kick| [urad]")
            kick_axis.legend(loc="upper left")
            kick_figure.tight_layout()
            kick_figure.savefig(
                config.output_dir / "kickmap_fit_domain_3d.png",
                dpi=180,
            )
            plt.close(kick_figure)
            kick_fit_plot_status = "written"

    origin_index = locate_origin_index(r_mm)
    if not valid_aperture_mask[origin_index]:
        raise RuntimeError("magnetic axis is not inside the valid aperture mask")
    on_axis_dynamic_kick = {
        "point_index": origin_index,
        "grid_row_y": int(rows[origin_index]),
        "grid_column_x": int(columns[origin_index]),
        "x_mm": float(x_mm[origin_index]),
        "y_mm": float(y_mm[origin_index]),
        "r_mm": float(r_mm[origin_index]),
        "kick_x_rad": float(theta_x_rad[origin_index]),
        "kick_y_rad": float(theta_y_rad[origin_index]),
        "kick_x_urad": float(theta_x_rad[origin_index] * 1.0e6),
        "kick_y_urad": float(theta_y_rad[origin_index] * 1.0e6),
    }

    valid_aperture_measured_max = float(
        np.max(measured_magnitude_urad[valid_aperture_mask])
    )
    valid_aperture_fitted_max = float(
        np.max(fitted_magnitude_urad[valid_aperture_mask])
    )
    fit_measured_max = float(np.max(measured_magnitude_urad[fit_mask]))
    fit_fitted_max = float(np.max(fitted_magnitude_urad[fit_mask]))
    formatted_output = result[5]

    return {
        "enabled": True,
        "method": (
            "RADIA FldFocKickPer periodic focusing-potential and "
            "second-order kick calculation"
        ),
        "radia_call": {
            "period_start_mm": period_start_mm,
            "longitudinal_direction": [0.0, 0.0, 1.0],
            "first_transverse_direction": [1.0, 0.0, 0.0],
            "second_transverse_direction": [0.0, -1.0, 0.0],
            "period_mm": config.period_mm,
            "full_periods": config.periods,
            "x_range_mm": x_range_mm,
            "x_points": x_points,
            "y_range_mm": y_range_mm,
            "y_points": y_points,
            "maximum_harmonic": config.kickmap_max_harmonic,
            "longitudinal_points_per_period": config.kickmap_points_per_period,
            "derivative_step_mm": config.kickmap_derivative_step_mm,
            "output_units": "rad",
            "formatted_output_format": "fix",
            "formatted_output_available": isinstance(formatted_output, str),
        },
        "shared_cartesian_field_kickmap_grid": {
            "coordinate_system": "regular Cartesian mesh",
            "x_min_mm": float(x_axis_mm[0]),
            "x_max_mm": float(x_axis_mm[-1]),
            "x_points": x_points,
            "x_spacing_mm": dx_mm,
            "y_min_mm": float(y_axis_mm[0]),
            "y_max_mm": float(y_axis_mm[-1]),
            "y_points": y_points,
            "y_spacing_mm": dy_mm,
            "total_points": total_points,
            "physical_aperture_shape": "rotated square (diamond in x-y)",
            "geometric_aperture_condition": (
                "abs(x_mm)+abs(y_mm) <= aperture_vertex_mm"
            ),
            "valid_physical_aperture_condition": (
                "abs(x_mm)+abs(y_mm) <= aperture_vertex_mm-boundary_margin_mm"
            ),
            "aperture_vertex_mm": config.aperture_vertex_mm,
            "boundary_margin_mm": boundary_margin_mm,
            "valid_aperture_limit_mm": valid_aperture_limit_mm,
            "points_inside_geometric_aperture": geometric_aperture_point_count,
            "points_inside_valid_physical_aperture": valid_aperture_point_count,
            "points_in_boundary_exclusion": boundary_exclusion_point_count,
            "points_outside_geometric_aperture": outside_aperture_point_count,
            "field_quantity_on_same_points": "B(x,y,z=0)",
            "kick_fit_r_max_mm": config.kick_fit_r_max_mm,
            "kick_fit_sampling_grid": (
                "all valid Cartesian points inside fit radius"
            ),
            "kick_fit_point_count": fit_point_count,
        },
        "radia_integrated_squared_transverse_field": {
            "description": (
                "third matrix returned by FldFocKickPer; exported and plotted "
                "as the RADIA focusing potential"
            ),
            "figure_7_comparison_quantity": True,
            "plot_domain": (
                "valid physical diamond with one-mesh-step boundary exclusion"
            ),
            "diagnostic_plot_domain": "full native Cartesian mesh",
            "minimum_on_full_cartesian_grid": float(
                np.min(btr_squared_native)
            ),
            "maximum_on_full_cartesian_grid": float(
                np.max(btr_squared_native)
            ),
            "minimum_inside_geometric_aperture": float(
                np.min(btr_squared_native[geometric_aperture_mask])
            ),
            "maximum_inside_geometric_aperture": float(
                np.max(btr_squared_native[geometric_aperture_mask])
            ),
            "minimum_inside_valid_physical_aperture": float(
                np.min(btr_squared_native[valid_aperture_mask])
            ),
            "maximum_inside_valid_physical_aperture": float(
                np.max(btr_squared_native[valid_aperture_mask])
            ),
            "minimum_inside_valid_physical_aperture_1e2_T2_mm2": float(
                np.min(btr_squared_native[valid_aperture_mask])
            ),
            "maximum_inside_valid_physical_aperture_1e2_T2_mm2": float(
                np.max(btr_squared_native[valid_aperture_mask])
            ),
            "definition": "<(integral Bx dz)^2 + (integral By dz)^2> over one undulator period",
            "plot_units": "1e2 T^2 mm^2",
            "equivalent_plot_units": "1e-4 T^2 m^2",
            "plot_scaling_from_RADIA_return": 1.0,
            "physical_value_per_plot_unit_T2_mm2": 1.0e2,
        },
        "particle": {
            "species_assumed_by_RADIA": "electron",
            "total_energy_gev": config.energy_gev,
        },
        "on_axis_dynamic_kick": on_axis_dynamic_kick,
        "fit": fit_record,
        "kick_extrema": {
            "fit_domain": {
                "r_max_mm": config.kick_fit_r_max_mm,
                "point_count": fit_point_count,
                "maximum_measured_kick_urad": fit_measured_max,
                "maximum_fitted_kick_urad": fit_fitted_max,
            },
            "valid_physical_aperture": {
                "shape": "rotated square with one-mesh-step boundary exclusion",
                "aperture_vertex_mm": config.aperture_vertex_mm,
                "boundary_margin_mm": boundary_margin_mm,
                "valid_aperture_limit_mm": valid_aperture_limit_mm,
                "point_count": valid_aperture_point_count,
                "maximum_measured_kick_urad": valid_aperture_measured_max,
                "maximum_fitted_kick_urad": valid_aperture_fitted_max,
                "fitted_value_note": (
                    "The polynomial is fitted only inside kick_fit_r_max_mm; "
                    "aperture-wide fitted extrema include extrapolation."
                ),
            },
        },
        # Backward-compatible aliases now refer only to the valid physical map.
        "maximum_measured_kick_urad": valid_aperture_measured_max,
        "maximum_fitted_kick_urad": valid_aperture_fitted_max,
        "diagnostics": {
            "full_cartesian_mesh_exported": True,
            "focusing_potential_and_kick_plots_are_separate": True,
            "focusing_potential_plot_quantity": (
                "third FldFocKickPer matrix (integrated squared transverse field)"
            ),
            "focusing_potential_plot_units": "1e2 T^2 mm^2",
            "focusing_potential_plot_rescaled": False,
            "focusing_potential_primary_plot_domain": (
                "valid diamond aperture with boundary exclusion"
            ),
            "focusing_potential_full_mesh_plot_is_diagnostic_only": True,
            "kick_plot_quantity": "hypot(kick_x, kick_y) inside fit radius",
            "outside_and_boundary_points_used_for_physical_extrema": False,
            "outside_and_boundary_points_used_for_fit": False,
            "diagnostic_csv_note": (
                "The diagnostic CSV contains raw rectangular-mesh values with "
                "point_region labels. It is not a beam-accessible kick map."
            ),
        },
        "outputs": {
            "kickmap_csv": str(output_path),
            "full_cartesian_diagnostic_csv": str(diagnostic_output_path),
            "focusing_potential_csv": str(focusing_potential_output_path),
            "fit_samples_csv": str(fit_output_path),
            "focusing_potential_plot": (
                str(config.output_dir / "focusing_potential_3d.png")
                if focusing_potential_plot_status == "written"
                else None
            ),
            "focusing_potential_plot_status": (
                focusing_potential_plot_status
            ),
            "focusing_potential_full_cartesian_diagnostic_plot": (
                str(
                    config.output_dir
                    / "focusing_potential_full_cartesian_diagnostic_3d.png"
                )
                if focusing_potential_diagnostic_plot_status == "written"
                else None
            ),
            "focusing_potential_full_cartesian_diagnostic_plot_status": (
                focusing_potential_diagnostic_plot_status
            ),
            "kickmap_fit_domain_plot": (
                str(config.output_dir / "kickmap_fit_domain_3d.png")
                if kick_fit_plot_status == "written"
                else None
            ),
            "kickmap_fit_domain_plot_status": kick_fit_plot_status,
            # Backward-compatible generic plot aliases now refer to the
            # fit-domain kick plot, never to the focusing-potential plot.
            "plot": (
                str(config.output_dir / "kickmap_fit_domain_3d.png")
                if kick_fit_plot_status == "written"
                else None
            ),
            "plot_status": kick_fit_plot_status,
        },
        "scope": (
            "periodic-body map from FldFocKickPer; finite entrance/exit end-field "
            "effects are not included in this periodic kick calculation"
        ),
    }

def run_radia(
    config: ModelConfig,
    specs: Sequence[MagnetSpec],
    manifest: dict[str, object],
) -> dict[str, object]:
    rad = import_radia(config.radia_pythonpath)
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

    material_handle: int | None = None
    solve_result: object | None = None
    if config.susceptibility:
        try:
            material_handle = int(rad.MatStd("NdFeB", config.br_t))
            for obj in objects:
                rad.MatApl(obj, material_handle)
            solve_result = rad.Solve(
                model, config.solve_precision, config.solve_max_iterations
            )
        except (AttributeError, RuntimeError, TypeError) as exc:
            raise RuntimeError(
                "RADIA susceptibility solve failed; verify the ESRF binding or "
                "use --no-susceptibility for a fixed-magnetization comparison"
            ) from exc

    z_min = min(s.center_z_mm - s.length_z_mm / 2.0 for s in specs)
    z_max = max(s.center_z_mm + s.length_z_mm / 2.0 for s in specs)
    z_mm = np.linspace(
        z_min - config.period_mm,
        z_max + config.period_mm,
        config.samples,
    )
    points = [[0.0, 0.0, float(z)] for z in z_mm]
    field = reshape_field(rad.Fld(model, "b", points), config.samples)
    if not np.all(np.isfinite(field)):
        raise RuntimeError("RADIA returned non-finite fields")

    with (config.output_dir / "on_axis_field.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["z_mm", "Bx_T", "By_T", "Bz_T"])
        for z, vector in zip(z_mm, field):
            writer.writerow([f"{z:.12g}", *(f"{v:.12g}" for v in vector)])

    half_width = config.period_mm * min(config.fit_periods, config.periods) / 2.0
    bx_fit = first_harmonic(z_mm, field[:, 0], config.period_mm, half_width)
    by_fit = first_harmonic(z_mm, field[:, 1], config.period_mm, half_width)
    period_cm = config.period_mm / 10.0
    kx = 0.934 * bx_fit["amplitude_t"] * period_cm
    ky = 0.934 * by_fit["amplitude_t"] * period_cm
    k_rms = math.sqrt((kx * kx + ky * ky) / 2.0)

    z_m = z_mm / 1000.0
    cum_bx = cumulative_trapezoid(field[:, 0], z_m)
    cum_by = cumulative_trapezoid(field[:, 1], z_m)
    i2_tm2 = {"Bx": trapezoid(cum_bx, z_m), "By": trapezoid(cum_by, z_m)}
    i2_gm2 = {key: value * 1.0e4 for key, value in i2_tm2.items()}

    gamma = config.energy_gev / 0.00051099895
    resonance_m = (config.period_mm * 1.0e-3) * (1.0 + k_rms**2) / (
        2.0 * gamma**2
    )
    checks = reference_checks(config, k_rms, kx, ky, i2_gm2)

    kickmap_analysis: dict[str, object] = {"enabled": False}
    if config.kickmap:
        kickmap_analysis = compute_second_order_kickmap(rad, model, config)

    analysis: dict[str, object] = {
        **manifest,
        "dry_run": False,
        "radia_model_handle": int(model),
        "radia_material_handle": material_handle,
        "radia_solve_result": solve_result,
        "sample_count": config.samples,
        "sample_range_mm": [float(z_mm[0]), float(z_mm[-1])],
        "fit_central_periods": min(config.fit_periods, config.periods),
        "first_harmonic": {"Bx": bx_fit, "By": by_fit},
        "K_peak_from_Bx": kx,
        "K_peak_from_By": ky,
        "K_rms_WEP38_convention": k_rms,
        "second_integral_T_m2": i2_tm2,
        "second_integral_G_m2": i2_gm2,
        "fundamental_resonance_nm": resonance_m * 1.0e9,
        "reference_checks": checks,
        "second_order_kickmap": kickmap_analysis,
    }

    if not config.no_plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib unavailable; plot skipped", file=sys.stderr)
        else:
            figure, axis = plt.subplots(figsize=(10, 5))
            axis.plot(z_mm, field[:, 0], label="Bx")
            axis.plot(z_mm, field[:, 1], label="By")
            axis.set_xlabel("z [mm]")
            axis.set_ylabel("B [T]")
            axis.set_title(
                f"AQUA APPLE-X 18 mm — {config.mode}, gap={config.gap_mm:g} mm"
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
    print("AQUA-style APPLE-X 18 mm — corrected WEP38 geometry")
    print(f"  implementation:       {SCRIPT_REVISION}")
    print(f"  mode:                 {config.mode}")
    print(f"  phase:                {config.phase_fraction:g} lambda_u")
    print(f"  period / periods:     {config.period_mm:g} mm / {config.periods}")
    print(f"  inter-girder gap:     {config.gap_mm:g} mm")
    print(
        "  aperture across flats:"
        f" {config.central_aperture_across_flats_mm:g} mm at min gap"
    )
    print(f"  derived inner cut:    {config.inner_chamfer_mm:.6g} mm")
    print(f"  aperture vertex:      +/-{config.aperture_vertex_mm:.6g} mm on x/y axes")
    print(f"  outer chamfer:        {config.outer_chamfer_mm:g} mm")
    print(f"  susceptibility solve: {config.susceptibility}")
    print(f"  second-order kickmap: {config.kickmap}")
    if config.kickmap:
        print(
            "  Cartesian grid x:     "
            f"[{config.field_x_min_mm:g}, {config.field_x_max_mm:g}] mm, "
            f"points={config.field_x_points}"
        )
        print(
            "  Cartesian grid y:     "
            f"[{config.field_y_min_mm:g}, {config.field_y_max_mm:g}] mm, "
            f"points={config.field_y_points}"
        )
        print(
            "  geometric aperture:   "
            f"|x|+|y| <= {config.aperture_vertex_mm:.6g} mm"
        )
        dx_mm = (config.field_x_max_mm - config.field_x_min_mm) / (
            config.field_x_points - 1
        )
        dy_mm = (config.field_y_max_mm - config.field_y_min_mm) / (
            config.field_y_points - 1
        )
        boundary_margin_mm = max(
            dx_mm,
            dy_mm,
            config.kickmap_derivative_step_mm,
        )
        print(
            "  valid physical map:   "
            f"|x|+|y| <= "
            f"{config.aperture_vertex_mm - boundary_margin_mm:.6g} mm "
            f"(margin={boundary_margin_mm:.6g} mm)"
        )
        print(
            "  RADIA native mesh:    "
            f"{config.field_x_points} x {config.field_y_points}, FldFocKickPer"
        )
        print(
            "  kick precision:       "
            f"harmonics={config.kickmap_max_harmonic}, "
            f"points/period={config.kickmap_points_per_period}"
        )
        print(f"  kick fit radius:      {config.kick_fit_r_max_mm:g} mm")
        print("  kick fit sampling:    complete Cartesian fit disk")
        print(f"  particle:             electron, E={config.energy_gev:g} GeV")
    print(f"  output directory:     {config.output_dir}")
    print(f"  dry run:              {config.dry_run}")


def main(argv: Sequence[str] | None = None) -> int:
    try:
        config = parse_args(argv)
        print_summary(config)
        result = execute(config)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    checks = result.get("reference_checks")
    if (
        config.strict_reference
        and isinstance(checks, dict)
        and checks.get("applicable") is True
        and checks.get("all_pass") is not True
    ):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
