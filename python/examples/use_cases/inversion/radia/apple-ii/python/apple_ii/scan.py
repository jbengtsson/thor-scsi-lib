from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path

from ._util import require_radia
from .cli_common import (
    add_geometry_arguments,
    finite_float,
    integer_at_least,
    nonnegative_float,
    parameters_from_args,
    parameters_to_dict,
    quality_from_args,
)
from .geometry import build_device
from .optimize import evaluate_phase, optimize_phase


def main() -> None:
    parser = argparse.ArgumentParser(description="Scan and refine APPLE-II row phase")
    add_geometry_arguments(parser)
    parser.add_argument("--phase-min", type=finite_float, default=-20.0)
    parser.add_argument("--phase-max", type=finite_float, default=20.0)
    parser.add_argument("--phase-steps", type=integer_at_least(2), default=41)
    parser.add_argument("--refinement-levels", type=integer_at_least(0), default=3)
    parser.add_argument("--refinement-steps", type=integer_at_least(3), default=11)
    parser.add_argument("--samples", type=integer_at_least(8), default=801)
    parser.add_argument("--central-periods", type=integer_at_least(1), default=6)
    parser.add_argument("--harmonics", type=integer_at_least(1), nargs="+", default=[1, 3, 5])
    parser.add_argument("--min-total-amplitude", type=nonnegative_float, default=1.0e-4)
    parser.add_argument("--max-relative-residual", type=nonnegative_float, default=0.10)
    parser.add_argument("--csv", default="apple2_phase_scan.csv")
    parser.add_argument("--json", default="apple2_phase_optimization.json")
    args = parser.parse_args()

    if args.phase_max <= args.phase_min:
        parser.error("--phase-max must be greater than --phase-min")
    if 1 not in args.harmonics:
        parser.error("--harmonics must include 1")
    parameters = parameters_from_args(args, 0.0)
    if args.central_periods > parameters.n_periods:
        parser.error("--central-periods must not exceed --periods")

    radia = require_radia()
    radia.UtiDelAll()
    device = build_device(parameters, radia_module=radia)
    half_length = 0.5 * (parameters.n_periods + 2) * parameters.period_mm
    criteria = quality_from_args(args)

    def evaluator(phase_mm: float):
        evaluation = evaluate_phase(
            device,
            phase_mm,
            z_min_mm=-half_length,
            z_max_mm=half_length,
            n_points=args.samples,
            central_periods=args.central_periods,
            harmonic_orders=args.harmonics,
            quality_criteria=criteria,
        )
        status = "PASS" if evaluation.quality.passed else "FAIL"
        print(
            f"{phase_mm:10.6f} mm  |s3|={evaluation.metrics['circularity']:.6f}  "
            f"B1={evaluation.quality.total_amplitude_t:.6g} T  {status}"
        )
        return evaluation

    result = optimize_phase(
        evaluator,
        args.phase_min,
        args.phase_max,
        coarse_steps=args.phase_steps,
        refinement_levels=args.refinement_levels,
        refinement_steps=args.refinement_steps,
    )

    all_evaluations = {evaluation.phase_mm: evaluation for evaluation in result.coarse_scan}
    all_evaluations.update(
        {evaluation.phase_mm: evaluation for evaluation in result.refinement_evaluations}
    )
    ordered = [all_evaluations[phase] for phase in sorted(all_evaluations)]
    records = [{"gap_mm": parameters.gap_mm, **evaluation.flat_record()} for evaluation in ordered]

    csv_path = Path(args.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    report = {
        "parameters": parameters_to_dict(parameters),
        "bounds_mm": result.bounds_mm,
        "quality_criteria": asdict(criteria),
        "best": {"gap_mm": parameters.gap_mm, **result.best.flat_record()},
        "coarse_count": len(result.coarse_scan),
        "refinement_evaluation_count": len(result.refinement_evaluations),
    }

    json_path = Path(args.json)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("best bounded point:")
    print(json.dumps(result.best.flat_record(), indent=2))


if __name__ == "__main__":
    main()
