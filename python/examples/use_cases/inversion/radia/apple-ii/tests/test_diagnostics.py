from __future__ import annotations

import math

from apple_ii.diagnostics import row_field_diagnostics
from apple_ii.field import sample_axis
from apple_ii.geometry import Apple2Parameters, build_device
from fake_radia import FakeRadia


class FieldFakeRadia(FakeRadia):
    def __init__(self):
        super().__init__()
        self.row_index = {}
        self.shifts = {}

    def ObjCnt(self, object_ids):
        object_id = super().ObjCnt(object_ids)
        children = self.containers[object_id]
        if children and all(child in self.blocks for child in children):
            self.row_index[object_id] = len(self.row_index)
        return object_id

    def TrfOrnt(self, object_id, transform):
        super().TrfOrnt(object_id, transform)
        self.shifts[object_id] = self.shifts.get(object_id, 0.0) + float(transform[1][2])
        return object_id

    def Fld(self, object_id, component, point):
        z_mm = float(point[2])
        if object_id in self.row_index:
            phase = 2.0 * math.pi * (z_mm - self.shifts.get(object_id, 0.0)) / 40.0
            return [0.0, 0.1 * math.cos(phase), 0.0]
        if object_id in self.containers:
            fields = [self.Fld(child, component, point) for child in self.containers[object_id]]
            return [sum(field[index] for field in fields) for index in range(3)]
        return [0.0, 0.0, 0.0]


def test_row_diagnostics_reports_four_rows_and_superposition():
    radia = FieldFakeRadia()
    device = build_device(
        Apple2Parameters(n_periods=2, end_block_fraction=0.0),
        radia_module=radia,
    )
    combined = sample_axis(device, -80.0, 80.0, 161)
    report = row_field_diagnostics(
        device,
        -80.0,
        80.0,
        161,
        period_mm=40.0,
        harmonic_orders=(1, 3, 5),
        central_periods=2,
        combined_rows=combined,
        radia_module=radia,
    )

    assert sorted(report["rows"]) == ["LL", "LR", "UL", "UR"]
    assert report["superposition"]["checked"] is True
    assert max(report["superposition"]["max_abs_error_T"].values()) < 1.0e-12
