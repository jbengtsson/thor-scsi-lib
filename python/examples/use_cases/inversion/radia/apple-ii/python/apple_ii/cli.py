from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

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
    parser.add_argument("--draw", action="store_true")
    args = parser.parse_args()

    if 1 not in args.harmonics:
        parser.error("--harmonics must include 1")
    parameters = parameters_from_args(args, args.phase_mm)
    if args.central_periods > parameters.n_periods:
        parser.error("--central-periods must not exceed --periods")

    radia = require_radia()
    radia.UtiDelAll()
    device = build_device(parameters, radia_module=radia)
    half_length = 0.5 * (parameters.n_periods + 2) * parameters.period_mm
    rows = sample_axis(device, -half_length, half_length, args.samples)
    write_csv(rows, args.csv)
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
    }
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
    for reason in quality.reasons:
        print(f"  - {reason}")
    if args.draw:
        radia.ObjDrwOpenGL(device.object_id)


if __name__ == "__main__":
    main()
