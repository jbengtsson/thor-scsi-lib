#!/usr/bin/env python3
"""RADIA model of an APPLE-II undulator in circular-polarization mode.

Coordinate convention
---------------------
    x : horizontal transverse direction [mm]
    y : vertical transverse direction [mm]
    z : electron-beam / undulator axis [mm]

The device consists of four pure-permanent-magnet Halbach arrays. Circular
polarization is obtained by shifting one diagonal pair longitudinally by
+/- period/4 relative to the other diagonal pair. The sign selects helicity.

This is a compact starting model, not a manufacturing design. End corrections,
magnet chamfers, measured magnet errors, support material, and trajectory / field
integral optimization are intentionally omitted.

Requires the RADIA Python module, commonly imported as ``radia``.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import math
from typing import Iterable, Sequence

try:
    import radia as rad
except ImportError as exc:  # pragma: no cover - depends on external RADIA install
    raise SystemExit(
        "RADIA is not installed or not visible to this Python interpreter. "
        "Install/build RADIA and ensure `import radia` works."
    ) from exc


@dataclass(frozen=True)
class Apple2Parameters:
    period_mm: float = 40.0
    n_periods: int = 10
    gap_mm: float = 12.0
    block_width_x_mm: float = 30.0
    block_height_y_mm: float = 12.0
    remanence_t: float = 1.20
    blocks_per_period: int = 4
    helicity: int = +1  # +1 or -1
    phase_fraction: float = 0.25  # circular mode: +/- lambda_u / 4

    @property
    def block_length_z_mm(self) -> float:
        return self.period_mm / self.blocks_per_period

    @property
    def phase_shift_mm(self) -> float:
        return self.helicity * self.phase_fraction * self.period_mm

    def validate(self) -> None:
        if self.period_mm <= 0 or self.gap_mm <= 0:
            raise ValueError("period_mm and gap_mm must be positive")
        if self.n_periods < 1:
            raise ValueError("n_periods must be at least 1")
        if self.blocks_per_period != 4:
            raise ValueError("This compact Halbach model requires 4 blocks per period")
        if self.helicity not in (-1, +1):
            raise ValueError("helicity must be +1 or -1")
        if self.block_width_x_mm <= 0 or self.block_height_y_mm <= 0:
            raise ValueError("magnet block dimensions must be positive")
        if self.remanence_t <= 0:
            raise ValueError("remanence_t must be positive")


def _halbach_magnetization(
    block_index: int,
    inward_axis: Sequence[float],
    longitudinal_sign: float,
    remanence_t: float,
) -> list[float]:
    """Return a four-block Halbach magnetization vector.

    The magnetization rotates by 90 degrees per block in the plane spanned by
    the array's inward normal and the longitudinal z axis.
    """
    angle = 0.5 * math.pi * block_index
    inward = math.cos(angle) * remanence_t
    longitudinal = math.sin(angle) * remanence_t * longitudinal_sign
    return [
        inward * inward_axis[0],
        inward * inward_axis[1],
        longitudinal,
    ]


def _make_array(
    p: Apple2Parameters,
    center_xy: Sequence[float],
    inward_axis: Sequence[float],
    phase_z_mm: float,
    longitudinal_sign: float,
    color: Sequence[float],
) -> int:
    """Create one finite Halbach array and return its RADIA container handle."""
    n_blocks = p.n_periods * p.blocks_per_period
    dz = p.block_length_z_mm
    z0 = -0.5 * (n_blocks - 1) * dz + phase_z_mm
    blocks: list[int] = []

    for i in range(n_blocks):
        center = [center_xy[0], center_xy[1], z0 + i * dz]
        dims = [p.block_width_x_mm, p.block_height_y_mm, dz]
        magnetization = _halbach_magnetization(
            i, inward_axis, longitudinal_sign, p.remanence_t
        )
        block = rad.ObjRecMag(center, dims, magnetization)
        rad.ObjDrwAtr(block, list(color), 0.001)
        blocks.append(block)

    return rad.ObjCnt(blocks)


def build_apple2(p: Apple2Parameters) -> int:
    """Build a four-array APPLE-II device in nominal circular mode.

    Diagonal pair A (upper-left and lower-right) is unshifted.
    Diagonal pair B (upper-right and lower-left) is shifted by +/- lambda_u/4.
    """
    p.validate()
    y_center = 0.5 * p.gap_mm + 0.5 * p.block_height_y_mm
    x_offset = 0.5 * p.block_width_x_mm

    # Array positions are separated left/right to expose the four-array APPLE-II
    # topology. The beam travels along z through x=y=0.
    specs = [
        # center (x,y), inward normal, z shift, Halbach rotation, RGB
        #
        # Facing upper/lower Halbach rows require opposite longitudinal
        # rotation senses.  With equal senses their on-axis fields cancel.
        # The mechanical APPLE-II phase is applied to one diagonal pair.
        ((-x_offset, +y_center), (0.0, -1.0), 0.0, +1.0, (0.85, 0.25, 0.25)),
        ((+x_offset, -y_center), (0.0, +1.0), 0.0, -1.0, (0.85, 0.25, 0.25)),
        ((+x_offset, +y_center), (0.0, -1.0), p.phase_shift_mm, +1.0, (0.25, 0.35, 0.90)),
        ((-x_offset, -y_center), (0.0, +1.0), p.phase_shift_mm, -1.0, (0.25, 0.35, 0.90)),
    ]

    arrays = [
        _make_array(p, center, inward, shift, rotation, color)
        for center, inward, shift, rotation, color in specs
    ]
    return rad.ObjCnt(arrays)


def sample_on_axis_field(
    device: int, z_min_mm: float, z_max_mm: float, n_points: int
) -> list[tuple[float, float, float]]:
    """Return [(z, Bx, By), ...] on the magnetic axis."""
    if n_points < 2:
        raise ValueError("n_points must be at least 2")
    # RADIA's Fld() evaluates the field at one point.  Some RADIA builds do
    # not support the multi-point calling form, so sample the axis explicitly.
    dz = (z_max_mm - z_min_mm) / (n_points - 1)
    rows: list[tuple[float, float, float]] = []
    for i in range(n_points):
        z = z_min_mm + i * dz
        field = rad.Fld(device, "b", [0.0, 0.0, z])
        if not isinstance(field, (list, tuple)) or len(field) < 2:
            raise RuntimeError(
                f"Unexpected RADIA field result at z={z:g} mm: {field!r}"
            )
        rows.append((z, float(field[0]), float(field[1])))
    return rows


def write_csv(rows: Iterable[tuple[float, float, float]], path: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as stream:
        stream.write("z_mm,Bx_T,By_T\n")
        for z, bx, by in rows:
            stream.write(f"{z:.9g},{bx:.9g},{by:.9g}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--period", type=float, default=40.0, help="period [mm]")
    parser.add_argument("--periods", type=int, default=10, help="number of periods")
    parser.add_argument("--gap", type=float, default=12.0, help="magnetic gap [mm]")
    parser.add_argument("--remanence", type=float, default=1.20, help="remanence [T]")
    parser.add_argument(
        "--helicity", type=int, choices=(-1, +1), default=+1,
        help="circular-polarization handedness selector",
    )
    parser.add_argument(
        "--phase-fraction", type=float, default=0.25,
        help="longitudinal shift as a fraction of one period; 0 gives planar baseline",
    )
    parser.add_argument("--samples", type=int, default=801, help="axis samples")
    parser.add_argument("--csv", default="apple2_circular_field.csv", help="output CSV")
    parser.add_argument("--draw", action="store_true", help="open RADIA 3D viewer")
    args = parser.parse_args()

    params = Apple2Parameters(
        period_mm=args.period,
        n_periods=args.periods,
        gap_mm=args.gap,
        remanence_t=args.remanence,
        helicity=args.helicity,
        phase_fraction=args.phase_fraction,
    )

    rad.UtiDelAll()
    device = build_apple2(params)

    half_scan = 0.5 * (params.n_periods + 2) * params.period_mm
    rows = sample_on_axis_field(device, -half_scan, +half_scan, args.samples)
    write_csv(rows, args.csv)

    peak_bx = max(abs(row[1]) for row in rows)
    peak_by = max(abs(row[2]) for row in rows)
    print(f"APPLE-II RADIA object: {device}")
    print(f"phase shift: {params.phase_shift_mm:.6g} mm")
    print(f"peak |Bx| on sampled axis: {peak_bx:.6g} T")
    print(f"peak |By| on sampled axis: {peak_by:.6g} T")
    if max(peak_bx, peak_by) < 1.0e-6:
        raise RuntimeError(
            "The assembled arrays produced numerical-zero on-axis field. "
            "Run once with --phase-fraction 0 to diagnose the planar baseline."
        )
    print(f"field data written to: {args.csv}")

    if args.draw:
        rad.ObjDrwOpenGL(device)


if __name__ == "__main__":
    main()
