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

def test_backtrack():
    calc_config = tslib.ConfigType()

    acc = create_acc()
    observers = instrument_with_standard_observers(acc, mapping=default_mapping)

    x = 1e-3
    y = -1e-3
    px = -2e-4
    py = -.3e-4
    ps_ref = gtpsa.ss_vect_double([x, px, y, py, 0, 0])
    ps = ps_ref.copy()

    n_elems = 160
    n_elems = 1000
    pscale = 1000e0

    # Forward tracking and plot
    next_element_index = acc.propagate(calc_config, ps, 0, n_elems)
    data = np.array([ob.get_phase_space().iloc for ob in observers][:n_elems]) * pscale
    _, axes = plt.subplots(2,2, sharex=True)
    ax_x, ax_y = axes[0]
    ax_x_diff, ax_y_diff = axes[1]
    line, = ax_x.plot(data[:, 0], '-')
    ax_y.plot(data[:, 2], '-', color=line.get_color())
    ax_x.set_ylabel("x [mm]")
    ax_y.set_ylabel("y [mm]")
    ax_x_diff.set_ylabel("dx [mm]")
    ax_y_diff.set_ylabel("dy [mm]")
    ref_data = data
    [ob.reset() for ob in observers]
    print(f"Forward tracking {next_element_index=} : {ps-ps_ref}")

    # Backward tracking
    ps2 = ps.copy()
    ps.px = -ps2.px
    ps.py = -ps2.py
    print(f"Result of tracking: {ps2} : starting back tracking with {ps}")
    next_element_index = acc.propagate(calc_config, ps, n_elems-1, -n_elems)
    data = np.array([ob.get_phase_space().iloc for ob in observers][:n_elems]) * pscale
    diff_data = data - ref_data
    line, = ax_x.plot(data[:, 0], '-')
    ax_y.plot(data[:, 2], '-', color=line.get_color())
    ax_x_diff.plot(diff_data[:, 0], '-', color=line.get_color())
    ax_y_diff.plot(diff_data[:, 2], '-', color=line.get_color())

    dps = ps-ps_ref
    print(f"Backward tracking {next_element_index=} : {dps}")

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
        test_backtrack()
    except:
        raise
    else:
        plt.ioff()
        plt.show()