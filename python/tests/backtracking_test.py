import gtpsa
import os.path
from thor_scsi.factory import accelerator_from_config
import thor_scsi.lib as tslib
from pathlib import Path
import pytest

path = (
    Path(os.environ["HOME"])
    / "Devel"
    / "gitlab"
    / "dt4acc"
    / "lattices/b2_stduser_beamports_blm_tracy_corr.lat"
)


def create_acc():
    return accelerator_from_config(path)


def test_backtrack():
    calc_config = tslib.ConfigType()

    acc = create_acc()
    ps_ref = gtpsa.ss_vect_double(0e0)
    ps_ref.x = 1e-3
    ps_ref.x = -1e-3

    ps = ps_ref.copy()

    n_elems = 160
    next_element_index = acc.propagate(calc_config, ps, 0, n_elems)
    print(f"Forward tracking {next_element_index=} : {ps-ps_ref}")
    ps2 = ps.copy()
    ps.px = -ps2.px
    ps.py = -ps2.py
    next_element_index = acc.propagate(calc_config, ps, n_elems - 1, -n_elems)

    dps = ps-ps_ref
    print(f"Backward tracking {next_element_index=} : {dps}")

    # fmt: off
    assert dps.x  == pytest.approx(0, abs=1e-5)
    assert dps.px == pytest.approx(0, abs=5e-6)
    assert dps.y  == pytest.approx(0, abs=1e-6)
    assert dps.py == pytest.approx(0, abs=1e-6)
    # fmt: on

