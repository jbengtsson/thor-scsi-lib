from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from types import ModuleType
from typing import Mapping, Sequence

from ._util import finite_number, nonnegative_number, positive_number, require_radia


@dataclass(frozen=True)
class ArrayMotion:
    """Longitudinal row displacement coefficients multiplied by ``phase_mm``."""

    ul: float = 0.0
    ur: float = 1.0
    ll: float = 1.0
    lr: float = 0.0

    @classmethod
    def elliptical(cls) -> "ArrayMotion":
        return cls()

    @classmethod
    def opposite_pair(cls) -> "ArrayMotion":
        return cls(-0.5, 0.5, 0.5, -0.5)

    def validate(self) -> None:
        for name, value in self.as_dict().items():
            finite_number(f"motion.{name}", value)

    def as_dict(self) -> dict[str, float]:
        return {"UL": self.ul, "UR": self.ur, "LL": self.ll, "LR": self.lr}


@dataclass(frozen=True)
class Apple2Parameters:
    period_mm: float = 40.0
    n_periods: int = 10
    gap_mm: float = 12.0
    row_width_x_mm: float = 30.0
    block_height_y_mm: float = 12.0
    remanence_t: float = 1.20
    blocks_per_period: int = 4
    phase_mm: float = 0.0
    motion: ArrayMotion = field(default_factory=ArrayMotion.elliptical)
    inter_array_gap_x_mm: float = 0.0
    end_block_fraction: float = 0.50
    end_block_length_mm: float | None = None
    end_magnetization_fraction: float = 0.50
    end_clearance_mm: float = 0.0

    @property
    def block_length_z_mm(self) -> float:
        return self.period_mm / self.blocks_per_period

    @property
    def end_block_length_z_mm(self) -> float:
        if self.end_block_length_mm is not None:
            return self.end_block_length_mm
        return self.block_length_z_mm * self.end_block_fraction

    @property
    def row_center_x_mm(self) -> float:
        """Absolute x coordinate of each row centre.

        ``inter_array_gap_x_mm`` is the horizontal face-to-face gap between the
        left and right arrays. A value of zero reproduces the original model.
        """
        return 0.5 * (self.row_width_x_mm + self.inter_array_gap_x_mm)

    def validate(self) -> None:
        positive_number("period_mm", self.period_mm)
        positive_number("gap_mm", self.gap_mm)
        positive_number("row_width_x_mm", self.row_width_x_mm)
        positive_number("block_height_y_mm", self.block_height_y_mm)
        positive_number("remanence_t", self.remanence_t)
        finite_number("phase_mm", self.phase_mm)
        nonnegative_number("inter_array_gap_x_mm", self.inter_array_gap_x_mm)
        nonnegative_number("end_block_fraction", self.end_block_fraction)
        if self.end_block_length_mm is not None:
            positive_number("end_block_length_mm", self.end_block_length_mm)
        nonnegative_number("end_magnetization_fraction", self.end_magnetization_fraction)
        nonnegative_number("end_clearance_mm", self.end_clearance_mm)
        if not isinstance(self.n_periods, int) or isinstance(self.n_periods, bool) or self.n_periods < 1:
            raise ValueError("n_periods must be an integer >= 1")
        if self.blocks_per_period != 4:
            raise ValueError("the current Halbach sequence requires blocks_per_period == 4")
        if self.end_block_fraction > 1:
            raise ValueError("end_block_fraction must be <= 1")
        if self.end_magnetization_fraction > 1:
            raise ValueError("end_magnetization_fraction must be <= 1")
        self.motion.validate()


@dataclass(frozen=True)
class BlockRecord:
    """Exact geometry and magnetization supplied to RADIA for one block."""

    object_id: int
    kind: str
    sequence_index: int
    center_mm: tuple[float, float, float]
    dimensions_mm: tuple[float, float, float]
    magnetization_t: tuple[float, float, float]

    def as_dict(self, current_shift_mm: float = 0.0) -> dict[str, object]:
        current_center = (
            self.center_mm[0],
            self.center_mm[1],
            self.center_mm[2] + current_shift_mm,
        )
        return {
            "object_id": self.object_id,
            "kind": self.kind,
            "sequence_index": self.sequence_index,
            "center_mm_zero_phase": list(self.center_mm),
            "center_mm_current": list(current_center),
            "dimensions_mm": list(self.dimensions_mm),
            "magnetization_t": list(self.magnetization_t),
        }


@dataclass
class RowHandle:
    name: str
    object_id: int
    motion_coefficient: float
    blocks: tuple[BlockRecord, ...] = ()
    current_shift_mm: float = 0.0

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "object_id": self.object_id,
            "motion_coefficient": self.motion_coefficient,
            "current_shift_mm": self.current_shift_mm,
            "block_count": len(self.blocks),
            "blocks": [block.as_dict(self.current_shift_mm) for block in self.blocks],
        }


@dataclass
class Apple2Device:
    """Persistent RADIA device and row handles supporting in-place phase motion."""

    object_id: int
    rows: dict[str, RowHandle]
    parameters: Apple2Parameters
    phase_mm: float
    _radia: object = field(repr=False)

    def set_phase(self, phase_mm: float) -> None:
        """Translate each row by only the increment needed for the requested phase."""
        phase_mm = finite_number("phase_mm", phase_mm)
        for row in self.rows.values():
            target = row.motion_coefficient * phase_mm
            delta = target - row.current_shift_mm
            if delta != 0.0:
                transform = self._radia.TrfTrsl([0.0, 0.0, delta])
                self._radia.TrfOrnt(row.object_id, transform)
                row.current_shift_mm = target
        self.phase_mm = phase_mm
        self.parameters = replace(self.parameters, phase_mm=phase_mm)

    def row_ids(self) -> Mapping[str, int]:
        return {name: row.object_id for name, row in self.rows.items()}

    def geometry_report(self) -> dict[str, object]:
        return {
            "device_object_id": self.object_id,
            "phase_mm": self.phase_mm,
            "row_center_x_mm": self.parameters.row_center_x_mm,
            "inter_array_gap_x_mm": self.parameters.inter_array_gap_x_mm,
            "resolved_end_block_length_mm": self.parameters.end_block_length_z_mm,
            "rows": {name: row.as_dict() for name, row in self.rows.items()},
        }


def object_id(device: Apple2Device | int) -> int:
    return device.object_id if isinstance(device, Apple2Device) else int(device)


def _mag(i: int, inward: Sequence[float], longitudinal_sign: float, br: float) -> list[float]:
    angle = 0.5 * math.pi * i
    normal = br * math.cos(angle)
    longitudinal = br * longitudinal_sign * math.sin(angle)
    return [normal * inward[0], normal * inward[1], longitudinal]


def _create_block(
    *,
    radia,
    center: Sequence[float],
    dimensions: Sequence[float],
    magnetization: Sequence[float],
    color: Sequence[float],
    kind: str,
    sequence_index: int,
) -> BlockRecord:
    block_id = int(radia.ObjRecMag(list(center), list(dimensions), list(magnetization)))
    radia.ObjDrwAtr(block_id, list(color), 0.001)
    return BlockRecord(
        object_id=block_id,
        kind=kind,
        sequence_index=sequence_index,
        center_mm=tuple(float(value) for value in center),
        dimensions_mm=tuple(float(value) for value in dimensions),
        magnetization_t=tuple(float(value) for value in magnetization),
    )


def _row(
    p: Apple2Parameters,
    x_mm: float,
    y_mm: float,
    inward: Sequence[float],
    longitudinal_sign: float,
    color: Sequence[float],
    radia,
) -> tuple[int, tuple[BlockRecord, ...]]:
    n_blocks = p.n_periods * p.blocks_per_period
    dz = p.block_length_z_mm
    first_center = -0.5 * (n_blocks - 1) * dz
    records: list[BlockRecord] = []

    for i in range(n_blocks):
        records.append(
            _create_block(
                radia=radia,
                center=[x_mm, y_mm, first_center + i * dz],
                dimensions=[p.row_width_x_mm, p.block_height_y_mm, dz],
                magnetization=_mag(i, inward, longitudinal_sign, p.remanence_t),
                color=color,
                kind="body",
                sequence_index=i,
            )
        )

    if (p.end_block_length_mm is not None or p.end_block_fraction > 0) and p.end_magnetization_fraction > 0:
        end_dz = p.end_block_length_z_mm
        clearance = p.end_clearance_mm
        lower_center = first_center - 0.5 * (dz + end_dz) - clearance
        last_center = first_center + (n_blocks - 1) * dz
        upper_center = last_center + 0.5 * (dz + end_dz) + clearance
        end_specs = ((-1, lower_center), (n_blocks, upper_center))
        for sequence_index, z_mm in end_specs:
            records.append(
                _create_block(
                    radia=radia,
                    center=[x_mm, y_mm, z_mm],
                    dimensions=[p.row_width_x_mm, p.block_height_y_mm, end_dz],
                    magnetization=_mag(
                        sequence_index,
                        inward,
                        longitudinal_sign,
                        p.remanence_t * p.end_magnetization_fraction,
                    ),
                    color=color,
                    kind="end",
                    sequence_index=sequence_index,
                )
            )

    row_id = int(radia.ObjCnt([record.object_id for record in records]))
    return row_id, tuple(records)


def build_device(
    parameters: Apple2Parameters,
    *,
    radia_module: ModuleType | object | None = None,
) -> Apple2Device:
    """Build four zero-phase rows once, then apply the requested phase in place."""
    parameters.validate()
    radia = require_radia(radia_module)
    x = parameters.row_center_x_mm
    y = 0.5 * (parameters.gap_mm + parameters.block_height_y_mm)
    # Mid-plane magnetic symmetry for the phase-zero planar mode:
    # transverse magnetization is the same in all four rows, while the
    # longitudinal sequence reverses between the upper and lower pairs.
    specs = (
        ("UL", -x, +y, (0.0, +1.0), +1.0, (0.85, 0.25, 0.25)),
        ("UR", +x, +y, (0.0, +1.0), +1.0, (0.25, 0.35, 0.90)),
        ("LL", -x, -y, (0.0, +1.0), -1.0, (0.25, 0.35, 0.90)),
        ("LR", +x, -y, (0.0, +1.0), -1.0, (0.85, 0.25, 0.25)),
    )
    coefficients = parameters.motion.as_dict()
    rows: dict[str, RowHandle] = {}
    for name, x_mm, y_mm, inward, sign, color in specs:
        row_id, blocks = _row(parameters, x_mm, y_mm, inward, sign, color, radia)
        rows[name] = RowHandle(name, row_id, coefficients[name], blocks)

    device_id = int(radia.ObjCnt([row.object_id for row in rows.values()]))
    device = Apple2Device(device_id, rows, parameters, 0.0, radia)
    device.set_phase(parameters.phase_mm)
    return device
