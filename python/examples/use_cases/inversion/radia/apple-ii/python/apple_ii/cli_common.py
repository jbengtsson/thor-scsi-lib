from __future__ import annotations

import argparse
from dataclasses import asdict
import math

from .analysis import QualityCriteria
from .geometry import Apple2Parameters, ArrayMotion


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be a finite number > 0")
    return number


def nonnegative_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("must be a finite number >= 0")
    return number


def finite_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise argparse.ArgumentTypeError("must be finite")
    return number


def integer_at_least(minimum: int):
    def parse(value: str) -> int:
        number = int(value)
        if number < minimum:
            raise argparse.ArgumentTypeError(f"must be an integer >= {minimum}")
        return number
    return parse


def add_geometry_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--period", type=positive_float, default=40.0)
    parser.add_argument("--periods", type=integer_at_least(1), default=10)
    parser.add_argument("--gap", type=positive_float, default=12.0)
    parser.add_argument("--width", type=positive_float, default=30.0)
    parser.add_argument("--height", type=positive_float, default=12.0)
    parser.add_argument("--remanence", type=positive_float, default=1.2)
    parser.add_argument("--motion", choices=("elliptical", "opposite"), default="elliptical")
    parser.add_argument("--end-block-fraction", type=nonnegative_float, default=0.5)
    parser.add_argument("--end-magnetization-fraction", type=nonnegative_float, default=0.5)
    parser.add_argument("--end-clearance-mm", type=nonnegative_float, default=0.0)


def parameters_from_args(args: argparse.Namespace, phase_mm: float) -> Apple2Parameters:
    motion = ArrayMotion.elliptical() if args.motion == "elliptical" else ArrayMotion.opposite_pair()
    return Apple2Parameters(
        period_mm=args.period,
        n_periods=args.periods,
        gap_mm=args.gap,
        row_width_x_mm=args.width,
        block_height_y_mm=args.height,
        remanence_t=args.remanence,
        phase_mm=phase_mm,
        motion=motion,
        end_block_fraction=args.end_block_fraction,
        end_magnetization_fraction=args.end_magnetization_fraction,
        end_clearance_mm=args.end_clearance_mm,
    )


def quality_from_args(args: argparse.Namespace) -> QualityCriteria:
    return QualityCriteria(args.min_total_amplitude, args.max_relative_residual)


def parameters_to_dict(parameters: Apple2Parameters) -> dict[str, object]:
    result = asdict(parameters)
    return result
