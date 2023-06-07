import numpy as np

import gtpsa
import os.path
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.linear_optics import instrument_with_standard_observers
import thor_scsi.lib as tslib
from pathlib import Path
import pytest
import matplotlib.pyplot as plt

path = (
    Path(os.environ["HOME"])
    / "Devel"
    / "gitlab"
    / "dt4acc"
    / "lattices/b2_stduser_beamports_blm_tracy_corr.lat"
)


def create_acc():
    return accelerator_from_config(path)


default_mapping = gtpsa._gtpsa.default_mapping()


def reverse_direction(
    ps: gtpsa.ss_vect_double, copy: bool = True
) -> gtpsa.ss_vect_tpsa:
    """Reverse direction of phase space

    Todo:
        where to put it into
    """
    if copy:
        nps = ps.copy()
    else:
        nps = ps

    nps.px = -ps.px
    nps.py = -ps.py
    # why in ct ... extra time for tardy electrons
    nps.ct = -ps.ct

    print(
        f"tracked started recevied  {ps}\nstarting backward with    {nps}",
    )

    return nps


def compare_backward_tracked_phase_space_to_forward(
    fwd_ps: gtpsa.ss_vect_double, bwd_ps: gtpsa.ss_vect_double
) -> None:
    """Compute difference between forward and backward phase space

    Args:
        fwd_ps: phase space that was used at the beginning to track forward
        bwd_ps: phase space that was recieved after tracking back the electrons

    """
    dps = bwd_ps - fwd_ps
    # double minus: one for returned direction
    #               one for comparison
    dps.px = bwd_ps.px + fwd_ps.px
    dps.py = bwd_ps.py + fwd_ps.py
    print(
        f"tracking started forward with  {fwd_ps}"
        f"backward tracking received     {bwd_ps}"
        f"difference reverting direction {dps}"
    )
    # fmt:off
    assert dps.x     == pytest.approx(0, abs=1e-12)
    assert dps.y     == pytest.approx(0, abs=1e-12)
    assert dps.px    == pytest.approx(0, abs=1e-12)
    assert dps.py    == pytest.approx(0, abs=1e-12)
    assert dps.delta == pytest.approx(0, abs=1e-12)
    # Why this deviation in ct
    assert dps.ct    == pytest.approx(0, abs=1e-8)
    # fmt:on


def test010_backtrack_drift():
    calc_config = tslib.ConfigType()

    acc = create_acc()
    drift = acc[1]
    assert isinstance(drift, tslib.Drift)

    x = 1e-3
    y = -1e-3
    px = -2e-4
    py = -0.3e-4
    ps_ref = gtpsa.ss_vect_double([x, px, y, py, 0, 0])
    ps = ps_ref.copy()

    drift.propagate(calc_config, ps)
    rps = reverse_direction(ps).copy()
    drift.propagate(calc_config, rps)

    compare_backward_tracked_phase_space_to_forward(ps_ref, rps)


def test020_backtrack_quadrupole():
    calc_config = tslib.ConfigType()

    acc = create_acc()
    quad = acc.find("q4m2d1r", 0)
    assert isinstance(quad, tslib.Quadrupole)

    x = 1e-3
    y = -1e-3
    px = -2e-4
    py = -0.3e-4
    ps_ref = gtpsa.ss_vect_double([x, px, y, py, 0, 0])
    ps = ps_ref.copy()

    quad.propagate(calc_config, ps)
    rps = reverse_direction(ps).copy()
    quad.propagate(calc_config, rps)

    compare_backward_tracked_phase_space_to_forward(ps_ref, rps)


def test030_backtrack_sextupole():
    calc_config = tslib.ConfigType()

    acc = create_acc()
    sext = acc.find("s3m2d1rr", 0)
    assert isinstance(sext, tslib.Sextupole)

    x = 1e-3
    y = -1e-3
    px = -2e-4
    py = -0.3e-4
    ps_ref = gtpsa.ss_vect_double([x, px, y, py, 0, 0])
    ps = ps_ref.copy()

    sext.propagate(calc_config, ps)
    rps = reverse_direction(ps)
    sext.propagate(calc_config, rps)

    compare_backward_tracked_phase_space_to_forward(ps_ref, rps)


def test040_backtrack_dipole():
    calc_config = tslib.ConfigType()

    acc = create_acc()
    dip = acc.find("bm2d1r12", 0)
    assert isinstance(dip, tslib.Bending)

    # without momentum
    x = 1e-3
    y = -1e-3
    px = -2e-4
    py = -0.3e-4
    ps_ref = gtpsa.ss_vect_double([x, px, y, py, 0, 0])
    ps = ps_ref.copy()

    dip.propagate(calc_config, ps)
    rps = reverse_direction(ps).copy()
    dip.propagate(calc_config, rps)

    compare_backward_tracked_phase_space_to_forward(ps_ref, rps)

    # off momentum
    x = 1e-3
    y = -1e-3
    px = -2e-4
    py = -0.3e-4
    delta = 1e-3
    ps_ref = gtpsa.ss_vect_double([x, px, y, py, 0, 0])
    ps_ref.delta = delta
    ps = ps_ref.copy()

    dip.propagate(calc_config, ps)
    rps = reverse_direction(ps).copy()
    dip.propagate(calc_config, rps)

    compare_backward_tracked_phase_space_to_forward(ps_ref, rps)

def test100_backtrack():
    calc_config = tslib.ConfigType()

    acc = create_acc()
    observers = instrument_with_standard_observers(acc, mapping=default_mapping)

    x = 1e-3
    y = -1e-3
    px = -2e-4
    py = -0.3e-4
    ps_ref = gtpsa.ss_vect_double([x, px, y, py, 0, 0])
    ps = ps_ref.copy()

    n_elems = 160
    n_elems = 2
    pscale = 1000e0

    # Forward tracking and plot
    print(f"Start tracking with {ps=}")
    next_element_index = acc.propagate(calc_config, ps, 1, n_elems)
    last_index = next_element_index - 1
    last_elem = acc[last_index]
    print(
        f"Last propagated element (returned index {last_index}): {last_elem.name}, {last_elem.index}"
    )

    data = np.array([ob.get_phase_space().iloc for ob in observers][:n_elems]) * pscale
    _, axes = plt.subplots(2, 2, sharex=True)
    ax_x, ax_y = axes[0]
    ax_x_diff, ax_y_diff = axes[1]
    (line,) = ax_x.plot(data[:, 0], "-")
    ax_y.plot(data[:, 2], "-", color=line.get_color())
    ax_x.set_ylabel("x [mm]")
    ax_y.set_ylabel("y [mm]")
    ax_x_diff.set_ylabel("dx [mm]")
    ax_y_diff.set_ylabel("dy [mm]")

    ref_data = data

    # Backward tracking
    rps = reverse_direction(ps)
    [ob.reset() for ob in observers]

    # fmt:off
    print(f"Result of tracking:         {ps}"
          f"starting back tracking with {rps}")
    # fmt:on
    next_element_index = acc.propagate(calc_config, rps, last_index, -n_elems)
    last_index = next_element_index + 1
    print(
        f"Last propagated element (returned index {last_index}): {last_elem.name}, {last_elem.index}"
    )
    last_elem = acc[last_index]
    print(f"Backward tracking {next_element_index=} : {rps}")
    compare_backward_tracked_phase_space_to_forward(ps_ref, rps)

    data = np.array([ob.get_phase_space().iloc for ob in observers][:n_elems]) * pscale
    diff_data = data - ref_data
    (line,) = ax_x.plot(data[:, 0], "-")
    ax_y.plot(data[:, 2], "-", color=line.get_color())
    ax_x_diff.plot(diff_data[:, 0], "-", color=line.get_color())
    ax_y_diff.plot(diff_data[:, 2], "-", color=line.get_color())

    return
    # fmt: off
    assert dps.x  == pytest.approx(0, abs=1e-5)
    assert dps.px == pytest.approx(0, abs=5e-6)
    assert dps.y  == pytest.approx(0, abs=1e-6)
    assert dps.py == pytest.approx(0, abs=1e-6)
    # fmt: on


if __name__ == "__main__":
    plt.ion()
    try:
        test100_backtrack()
    except:
        raise
    else:
        plt.ioff()
        plt.show()
