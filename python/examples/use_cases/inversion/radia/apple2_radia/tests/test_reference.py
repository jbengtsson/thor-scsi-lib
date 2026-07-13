from pathlib import Path

import pytest

from apple2.reference import ReferencePoint, compare_points, load_points


def test_reference_comparison_wraps_phase_and_reports_rmse():
    reference = [ReferencePoint(12.0, 10.0, 0.16, 0.29, 179.0, 0.85)]
    calculated = [ReferencePoint(12.0, 10.0, 0.17, 0.28, -179.0, 0.80)]
    points, summary = compare_points(calculated, reference)
    assert points[0].phase_error_deg == pytest.approx(2.0)
    assert summary["Bx_RMSE_T"] == pytest.approx(0.01)
    assert summary["circularity_RMSE"] == pytest.approx(0.05)


def test_load_points_rejects_missing_columns(tmp_path: Path):
    path = tmp_path / "bad.csv"
    path.write_text("gap_mm,phase_mm\n12,10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing required"):
        load_points(path)


def test_compare_requires_unique_coordinate_match():
    point = ReferencePoint(12.0, 10.0, 0.16, 0.29, 90.0, 0.85)
    with pytest.raises(ValueError, match="exactly one"):
        compare_points([], [point])
