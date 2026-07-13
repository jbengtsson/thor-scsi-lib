from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Iterable, Mapping


REQUIRED_COLUMNS = ("gap_mm", "phase_mm", "Bx1_T", "By1_T", "phase_deg", "circularity")


@dataclass(frozen=True)
class ReferencePoint:
    gap_mm: float
    phase_mm: float
    Bx1_T: float
    By1_T: float
    phase_deg: float
    circularity: float


@dataclass(frozen=True)
class ComparisonPoint:
    gap_mm: float
    phase_mm: float
    bx_error_t: float
    by_error_t: float
    phase_error_deg: float
    circularity_error: float


def _wrapped_phase_error_deg(value: float, reference: float) -> float:
    return (value - reference + 180.0) % 360.0 - 180.0


def load_points(path: str | Path) -> list[ReferencePoint]:
    with Path(path).open(newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        missing = [column for column in REQUIRED_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"missing required CSV columns: {', '.join(missing)}")
        points: list[ReferencePoint] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                values = {column: float(row[column]) for column in REQUIRED_COLUMNS}
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid numeric value at CSV line {line_number}") from exc
            if not all(math.isfinite(value) for value in values.values()):
                raise ValueError(f"non-finite value at CSV line {line_number}")
            points.append(ReferencePoint(**values))
    if not points:
        raise ValueError("reference CSV contains no data rows")
    return points


def compare_points(
    calculated: Iterable[ReferencePoint],
    reference: Iterable[ReferencePoint],
    *,
    coordinate_tolerance_mm: float = 1.0e-9,
) -> tuple[list[ComparisonPoint], dict[str, float | int]]:
    if coordinate_tolerance_mm < 0 or not math.isfinite(coordinate_tolerance_mm):
        raise ValueError("coordinate_tolerance_mm must be finite and >= 0")
    calculated_points = list(calculated)
    reference_points = list(reference)
    comparisons: list[ComparisonPoint] = []
    used: set[int] = set()
    for ref in reference_points:
        matches = [
            (index, point)
            for index, point in enumerate(calculated_points)
            if index not in used
            and abs(point.gap_mm - ref.gap_mm) <= coordinate_tolerance_mm
            and abs(point.phase_mm - ref.phase_mm) <= coordinate_tolerance_mm
        ]
        if len(matches) != 1:
            raise ValueError(
                f"expected exactly one calculated point at gap={ref.gap_mm}, "
                f"phase={ref.phase_mm}; found {len(matches)}"
            )
        index, point = matches[0]
        used.add(index)
        comparisons.append(
            ComparisonPoint(
                gap_mm=ref.gap_mm,
                phase_mm=ref.phase_mm,
                bx_error_t=point.Bx1_T - ref.Bx1_T,
                by_error_t=point.By1_T - ref.By1_T,
                phase_error_deg=_wrapped_phase_error_deg(point.phase_deg, ref.phase_deg),
                circularity_error=point.circularity - ref.circularity,
            )
        )

    def rmse(attribute: str) -> float:
        return math.sqrt(
            sum(getattr(point, attribute) ** 2 for point in comparisons) / len(comparisons)
        )

    summary: dict[str, float | int] = {
        "count": len(comparisons),
        "Bx_RMSE_T": rmse("bx_error_t"),
        "By_RMSE_T": rmse("by_error_t"),
        "phase_RMSE_deg": rmse("phase_error_deg"),
        "circularity_RMSE": rmse("circularity_error"),
    }
    return comparisons, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare calculated APPLE-II metrics to reference data")
    parser.add_argument("calculated_csv")
    parser.add_argument("reference_csv")
    parser.add_argument("--coordinate-tolerance-mm", type=float, default=1.0e-9)
    parser.add_argument("--json", default="apple2_reference_comparison.json")
    args = parser.parse_args()
    comparisons, summary = compare_points(
        load_points(args.calculated_csv),
        load_points(args.reference_csv),
        coordinate_tolerance_mm=args.coordinate_tolerance_mm,
    )
    report = {
        "summary": summary,
        "points": [asdict(point) for point in comparisons],
    }
    Path(args.json).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
