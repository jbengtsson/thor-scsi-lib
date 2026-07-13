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
    end_block_fraction: float = 0.50
    end_magnetization_fraction: float = 0.50
    end_clearance_mm: float = 0.0

    @property
    def block_length_z_mm(self) -> float:
        return self.period_mm / self.blocks_per_period

    @property
    def end_block_length_z_mm(self) -> float:
        return self.block_length_z_mm * self.end_block_fraction

    def validate(self) -> None:
        positive_number("period_mm", self.period_mm)
        positive_number("gap_mm", self.gap_mm)
        positive_number("row_width_x_mm", self.row_width_x_mm)
        positive_number("block_height_y_mm", self.block_height_y_mm)
        positive_number("remanence_t", self.remanence_t)
        finite_number("phase_mm", self.phase_mm)
        nonnegative_number("end_block_fraction", self.end_block_fraction)
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


@dataclass
class RowHandle:
    name: str
    object_id: int
    motion_coefficient: float
    current_shift_mm: float = 0.0


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


def object_id(device: Apple2Device | int) -> int:
    return device.object_id if isinstance(device, Apple2Device) else int(device)


def _mag(i: int, inward: Sequence[float], longitudinal_sign: float, br: float) -> list[float]:
    angle = 0.5 * math.pi * i
    normal = br * math.cos(angle)
    longitudinal = br * longitudinal_sign * math.sin(angle)
    return [normal * inward[0], normal * inward[1], longitudinal]


def _row(
    p: Apple2Parameters,
    x_mm: float,
    y_mm: float,
    inward: Sequence[float],
    longitudinal_sign: float,
    color: Sequence[float],
    radia,
) -> int:
    n_blocks = p.n_periods * p.blocks_per_period
    dz = p.block_length_z_mm
    first_center = -0.5 * (n_blocks - 1) * dz
    blocks: list[int] = []

    for i in range(n_blocks):
        block = radia.ObjRecMag(
            [x_mm, y_mm, first_center + i * dz],
            [p.row_width_x_mm, p.block_height_y_mm, dz],
            _mag(i, inward, longitudinal_sign, p.remanence_t),
        )
        radia.ObjDrwAtr(block, list(color), 0.001)
        blocks.append(block)

    if p.end_block_fraction > 0 and p.end_magnetization_fraction > 0:
        end_dz = p.end_block_length_z_mm
        clearance = p.end_clearance_mm
        lower_center = first_center - 0.5 * (dz + end_dz) - clearance
        last_center = first_center + (n_blocks - 1) * dz
        upper_center = last_center + 0.5 * (dz + end_dz) + clearance
        end_specs = ((-1, lower_center), (n_blocks, upper_center))
        for sequence_index, z_mm in end_specs:
            block = radia.ObjRecMag(
                [x_mm, y_mm, z_mm],
                [p.row_width_x_mm, p.block_height_y_mm, end_dz],
                _mag(
                    sequence_index,
                    inward,
                    longitudinal_sign,
                    p.remanence_t * p.end_magnetization_fraction,
                ),
            )
            radia.ObjDrwAtr(block, list(color), 0.001)
            blocks.append(block)

    return int(radia.ObjCnt(blocks))


def build_device(
    parameters: Apple2Parameters,
    *,
    radia_module: ModuleType | object | None = None,
) -> Apple2Device:
    """Build four zero-phase rows once, then apply the requested phase in place."""
    parameters.validate()
    radia = require_radia(radia_module)
    x = 0.5 * parameters.row_width_x_mm
    y = 0.5 * (parameters.gap_mm + parameters.block_height_y_mm)
    specs = (
        ("UL", -x, +y, (0.0, -1.0), +1.0, (0.85, 0.25, 0.25)),
        ("UR", +x, +y, (0.0, -1.0), +1.0, (0.25, 0.35, 0.90)),
        ("LL", -x, -y, (0.0, +1.0), -1.0, (0.25, 0.35, 0.90)),
        ("LR", +x, -y, (0.0, +1.0), -1.0, (0.85, 0.25, 0.25)),
    )
    coefficients = parameters.motion.as_dict()
    rows: dict[str, RowHandle] = {}
    for name, x_mm, y_mm, inward, sign, color in specs:
        row_id = _row(parameters, x_mm, y_mm, inward, sign, color, radia)
        rows[name] = RowHandle(name, row_id, coefficients[name])

    device_id = int(radia.ObjCnt([row.object_id for row in rows.values()]))
    device = Apple2Device(device_id, rows, parameters, 0.0, radia)
    device.set_phase(parameters.phase_mm)
    return device
