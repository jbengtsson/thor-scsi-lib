"""Check on closed orbit computation

Current checks:
    to be used as beam based alignment simulation back bone
"""
from dataclasses import dataclass
from pathlib import Path
import os

import numpy as np
import pytest
import thor_scsi.lib as tslib

import gtpsa
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.closed_orbit import compute_closed_orbit

path = (
    Path(os.environ["HOME"])
    / "Devel"
    / "gitlab"
    / "dt4acc"
    / "lattices/b2_stduser_beamports_blm_tracy_corr.lat"
)


def create_acc():
    return accelerator_from_config(path)


def test10_closed_orbit_for_zero():
    calc_config = tslib.ConfigType()
    acc = create_acc()

    # compute for a start at zero
    result = compute_closed_orbit(acc, calc_config, delta=0e0)
    assert result.found_closed_orbit == True
    x0 = result.x0
    assert x0.x == pytest.approx(0, abs=1e-9)
    assert x0.px == pytest.approx(0, abs=1e-9)
    assert x0.y == pytest.approx(0, abs=1e-9)
    assert x0.py == pytest.approx(0, abs=1e-9)
    assert x0.delta == pytest.approx(0, abs=1e-9)
    # should be only 2 degrees of freedom?
    ct_save = x0.ct
    assert ct_save == pytest.approx(0, abs=1e-4)
    # assert x0.ct == pytest.approx(0,abs=1e-9)

    result = compute_closed_orbit(acc, calc_config, x0=x0, delta=None)
    assert result.found_closed_orbit == True
    x0 = result.x0
    assert x0.x == pytest.approx(0, abs=1e-9)
    assert x0.px == pytest.approx(0, abs=1e-9)
    assert x0.y == pytest.approx(0, abs=1e-9)
    assert x0.py == pytest.approx(0, abs=1e-9)
    assert x0.delta == pytest.approx(0, abs=1e-9)
    assert x0.ct == pytest.approx(ct_save, abs=1e-9)


def compensate_quadrupole_offset_with_dipole_kick(quad):
    dx = quad.get_dx()
    assert quad.get_main_multipole_number() == 2
    muls = quad.get_multipoles()
    k = muls.get_multipole(2)
    b = k * dx
    bref = muls.get_multipole(1)
    new_b = bref + b
    muls.set_multipole(1, new_b)
    assert quad.get_multipoles().get_multipole(1).real == pytest.approx(new_b, 1e-12)


def test20_closed_orbit_for_moved_quad_with_dipole_kick():
    """quadrupole compensated by dipole kick"""
    calc_config = tslib.ConfigType()
    acc = create_acc()

    quad_name = "q1m1d1r"
    quad = acc.find(quad_name, 0)

    dx = 1e-4
    compensate_quadrupole_offset_with_dipole_kick(quad)

    result = compute_closed_orbit(acc, calc_config, delta=0e0)
    assert result.found_closed_orbit == True
    x0 = result.x0
    assert x0.x == pytest.approx(0, abs=1e-9)
    assert x0.px == pytest.approx(0, abs=1e-9)
    assert x0.y == pytest.approx(0, abs=1e-9)
    assert x0.py == pytest.approx(0, abs=1e-9)
    assert x0.delta == pytest.approx(0, abs=1e-9)
    assert x0.ct == pytest.approx(0, abs=1e-9)

    dx = 2e-4
    result = compute_closed_orbit(acc, calc_config, delta=0e0)
    assert result.found_closed_orbit == True
    x0 = result.x0
    assert x0.x == pytest.approx(0, abs=1e-9)
    assert x0.px == pytest.approx(0, abs=1e-9)
    assert x0.y == pytest.approx(0, abs=1e-9)
    assert x0.py == pytest.approx(0, abs=1e-9)
    assert x0.delta == pytest.approx(0, abs=1e-9)
    assert x0.ct == pytest.approx(0, abs=1e-9)


@dataclass
class DistortedOrbit:
    # offset of quad
    dy: float
    # fixed point in x
    y: float
    # fixed point in x
    py: float


def test30_closed_orbit_for_moved_quad():
    """Check that xo, px scale linearly with moved quad

    Vertical plane: avoid dispersion effects if so
    """
    calc_config = tslib.ConfigType()
    acc = create_acc()

    quad_name = "q2m2d2r"
    quad = acc.find(quad_name, 0)

    ref = None
    dy_vals = [1e-4] + (np.arange(-3, 4) * 1e-4).tolist()
    for cnt, dy in enumerate(dy_vals):
        quad.set_dy(dy)

        result = compute_closed_orbit(acc, calc_config, delta=0e0, eps=1e-9)
        assert result.found_closed_orbit == True
        x0 = result.x0
        print(f"{cnt=}, {dy=}, {x0=}")
        if cnt == 0:
            assert x0.y == pytest.approx(0e0, abs=1e-3)
            assert x0.py == pytest.approx(0e0, abs=1e-3)
            ref = DistortedOrbit(dy=dy, y=x0.y, py=x0.py)
        else:
            scale = dy / ref.dy
            print(f"{ref=} {scale=} ")

            eps_rel = 1e-4 * abs(scale)
            assert x0.y == pytest.approx(ref.y * scale, abs=1e-9, rel=eps_rel)
            assert x0.py == pytest.approx(ref.py * scale, abs=1e-9, rel=eps_rel)

            if abs(scale) > 1e-9:
                print(
                    f"y/y_ref/scale = {x0.y/ref.y/scale}, py/py_ref/scale = {x0.py/ref.py/scale}"
                )

        assert x0.x == pytest.approx(0, abs=1e-5)
        assert x0.px == pytest.approx(0, abs=3e-6)
        assert x0.delta == pytest.approx(0, abs=1e-9)
        assert x0.ct == pytest.approx(0, abs=1e-4)
