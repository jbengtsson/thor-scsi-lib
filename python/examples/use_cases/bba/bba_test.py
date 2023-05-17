import logging
import sys
from dataclasses import dataclass
from pathlib import Path
import os
from typing import Sequence

import pytest

import gtpsa
import numpy as np
import scipy.optimize
from sklearn.linear_model import Ridge
import xarray as xr
import thor_scsi.lib as tslib
from bact_math_utils.distorted_orbit import closed_orbit_distortion
from bact_analysis.transverse.distorted_orbit import closed_orbit_distortion as co_dist
import matplotlib.pyplot as plt

import thor_scsi.utils.knobs
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.accelerator import (
    instrument_with_standard_observers,
    extract_orbit_from_standard_observers,
)
from thor_scsi.utils.extract_info import accelerator_info
from thor_scsi.utils.linear_optics import compute_Twiss_along_lattice

logger = logging.getLogger("thor-scsi")

t_file = (
    Path(os.environ["HOME"])
    / "Devel"
    / "gitlab"
    / "dt4acc"
    / "lattices"
    / "b2_stduser_beamports_blm_tracy_corr.lat"
)


named_index_d = dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6, dx=7, dy=8)
named_index = gtpsa.IndexMapping(named_index_d)


@dataclass
class AcceleratorForTest:
    acc: tslib.Accelerator
    acc_info: xr.Dataset
    quad: tslib.Quadrupole
    observers: Sequence[tslib.StandardObserver]
    bpm_observers: Sequence[tslib.StandardObserver]
    bpm: Sequence[tslib.BPM]
    observed_elements: Sequence[tslib.ElemType]


def create_accelerator():
    """create a new one for each test

    Returns: AcceleratorForTest

    """
    print(f"Reading lattice file {t_file}")
    acc = accelerator_from_config(t_file)
    acc_info = accelerator_info(acc)

    # Required to see what happens in bpms
    quad = acc.find("q4m1d1r", 0)
    # quad = acc.find("q1m1d1r", 0)
    # quad = acc.find("q3m2d8r", 0)
    # and the elements before and after
    before_quad = acc[quad.index - 1]
    after_quad = acc[quad.index + 1]
    bpm_elems = [elem for elem in acc if isinstance(elem, tslib.BPM)]
    observed_elements = bpm_elems + [before_quad, quad, after_quad]
    observers = instrument_with_standard_observers(observed_elements, mapping=named_index)
    bpm_observers = observers[: len(bpm_elems)]
    return AcceleratorForTest(
        acc=acc,
        acc_info=acc_info,
        quad=quad,
        observers=observers,
        bpm_observers=bpm_observers,
        bpm=bpm_elems,
        observed_elements=observed_elements,
    )


def func_find_fixed_point(ps_start, acc, calc_config):
    ps = gtpsa.ss_vect_double(0e0)
    tmp = ps.iloc
    for cnt, val in enumerate(ps_start):
        tmp[cnt] = val
    acc.propagate(calc_config, ps)
    offset = np.array(ps.iloc) - ps_start
    assert offset is not None
    return offset


def find_fixed_point(acc: tslib.Accelerator, *, calc_config=tslib.ConfigType()):
    x0 = np.zeros(6)
    result = scipy.optimize.root(
        func_find_fixed_point, x0, tol=1e-13, args=(acc, calc_config)
    )
    return result


def dipole_due_to_quadrupole_gradient_change(
    quadrupole: tslib.Quadrupole, strength_change: float, offset: complex
) -> complex:
    """
    Args:
        strength change: expected to be precentage of main

    """
    gradient = quadrupole.get_main_multipole_strength()
    length = quadrupole.get_length()

    added_dipole_kick = strength_change * gradient * offset * length
    return added_dipole_kick


desc = gtpsa.desc(6, 2, 3, 1)

def quadrupole_response(dc: AcceleratorForTest):
    # Find the quadrupole in question
    t_idx = dc.quad.index
    t_quad = dc.acc[t_idx]
    # quadrupole with knobs

    def convert(elem):
        if elem.index == t_idx:
            elem = thor_scsi.utils.knobs.convert_if_quadrupole(elem)
            return thor_scsi.utils.knobs.make_magnet_knobbable(elem, po=1, desc=desc, named_index=named_index, offset=True)
        return elem
    elems = [convert(elem) for elem in dc.acc]

    acc_tpsa = tslib.AcceleratorTpsa(elems)
    observers = [tslib.StandardObserver(named_index) for elem in acc_tpsa]
    for elem, ob in zip(acc_tpsa, observers):
        elem.set_observer(ob)

    calc_config = tslib.ConfigType()
    ps = gtpsa.ss_vect_tpsa(desc, 2, 6, named_index)
    ps.set_identity()
    # Propagate it once to get the effect of the quad
    acc_tpsa.propagate(calc_config, ps)
    # Propagate it a second time to store it in the bpms ...
    acc_tpsa.propagate(calc_config, ps)
    # a few more to weed out the effect
    for cnt in range(0):
        acc_tpsa.propagate(calc_config, ps)

    for elem in acc_tpsa:
        assert(elem.get_observer().has_truncated_power_series_a())
    return acc_tpsa, observers, ps

# @pytest.mark.skip
def test01_quadrupole_moved_steerer_compensation():
    dc = create_accelerator()
    desc = gtpsa.desc(6, 2, 3, 1)

    # acc = accelerator_from_config(t_file)
    twiss = compute_Twiss_along_lattice(
        2, dc.acc, desc=desc, mapping=named_index
    )

    acc_tpsa, observer_tpsa, ps_orig = quadrupole_response(dc)
    calc_config = tslib.ConfigType()
    ps = ps_orig.copy()
    ps_ref = ps.copy()
    if True:
        # why that extra iteration?
        for i in range(1):
            acc_tpsa.propagate(calc_config, ps)
        # should not be necessary
        for ob in observer_tpsa:
            ob.reset()
        # need to run it again ... to get the effect "into the observers"
        acc_tpsa.propagate(calc_config, ps)

    def extract_forecast_from_observer(elem):
        observer = elem.get_observer()
        # observer = observer_tpsa[elem.index]
        try:
            ps = observer.get_truncated_power_series_a() #.copy()
        except:
            logger.error(f"While evaluating element {elem}")
            raise
        res = [ps.x.get(dx=1), ps.y.get(dy=1), observer.get_observed_name(), observer.get_observed_index()]
        return res

    forecast_displacement = np.array(
        [extract_forecast_from_observer(elem) for elem in acc_tpsa], dtype=object
    )

    def track_to_element(ps, start, n_elems):
        # print(f"{named_index=}")
        # print("ps.x.get_mapping()", ps.x.get_mapping())
        assert(n_elems == 1)
        next_elem = acc_tpsa.propagate(calc_config, ps, start, n_elems)
        assert(next_elem == start + n_elems)
        elem = acc_tpsa[start]
        # print(f"{elem=} {ps}")
        # print("ps.x.get_mapping()", ps.x.get_mapping())
        return [ps.x.get(dx=1), ps.y.get(dy=1), elem.name, elem.index]


    if True:
        ps = gtpsa.ss_vect_tpsa(desc, 2, 6, named_index)
        ps.set_identity()
        for i in range(3):
            acc_tpsa.propagate(calc_config, ps)
        forecast_displacement2 = np.array([track_to_element(ps, start, 1) for start in range(len(acc_tpsa))], dtype=np.object)
        forecast_displacement2_s = np.add.accumulate([elem.get_length() for elem in acc_tpsa])
    t_nu_x = twiss.twiss.sel(par="dnu", plane="x").cumsum()  # / (2*np.pi)
    print(t_nu_x)
    pi2 = np.pi * 2
    co = closed_orbit_distortion(
        twiss.twiss.sel(par="beta", plane="x").values,
        t_nu_x.values * pi2,
        tune=float(t_nu_x.isel(index=-1)),
        theta_i=1e-5,
        beta_i=twiss.twiss.sel(par="beta", plane="x", index=dc.quad.index).values,
        mu_i=t_nu_x.sel(index=dc.quad.index).values * pi2,
    )

    phase_advance_plot = False
    twiss_plot = False
    closed_orbit_plot = False

    if phase_advance_plot:
        fig, ax = plt.subplots(1, 1)
        ax.plot(twiss.s, t_nu_x)

    if twiss_plot:
        fig, axes = plt.subplots(2, 1, sharex=True)
        ax_nu, ax_s = axes
        ax_nu.plot(t_nu_x, twiss.twiss.sel(par="beta", plane="x"), "-")
        ax_nu.plot(t_nu_x, twiss.twiss.sel(par="beta", plane="y"), "-")
        ax_s.plot(twiss.s, twiss.twiss.sel(par="beta", plane="x"), "-")
        ax_s.plot(twiss.s, twiss.twiss.sel(par="beta", plane="y"), "-")

    if closed_orbit_plot:
        fig, axes = plt.subplots(2, 1, sharex=True)
        ax_nu, ax_s = axes
        ax_nu.plot(t_nu_x, co, "-")
        ax_s.plot(twiss.s, co, "-")
        # return

        # ax.plot(nu_x, observed_displacement[:, 0] * 1e8, "bx")
        # ax.plot(nu_x, displaced_quad_effect.sel(dep="dx", plane="x").values * 1, "g.")
        # ax.plot(nu_x, displaced_quad_effect.sel(dep="dKdx", plane="x").values * 1 * 4, "r.")
        fig.savefig("closed_orbit_distortion_x.pdf")


    # chose an offset
    dx = 1e-4
    # dx = 0e0
    quad = dc.quad
    # set the offset
    quad.set_dx(dx)
    # add the quadrupole a steerer component compensationg the offset
    # then the beam should go through undisturbed
    offset = quad.get_dx() + quad.get_dy() * 1j
    quad_muls = quad.get_multipoles()
    dipole_from_feed_down = quad_muls.get_multipole(2) * offset
    quad_muls.set_multipole(1, dipole_from_feed_down)

    assert quad.get_dx() == pytest.approx(dx, 1e-12)
    observers = instrument_with_standard_observers(
        dc.observed_elements, mapping=named_index
    )
    ps = gtpsa.ss_vect_double(0e0)
    acc = dc.acc
    acc.propagate(calc_config, ps)

    if False:
        # Check that the steerer compensates the offset effect
        print(type(ps))
        print(dir(ps))
        assert ps.x == pytest.approx(0e0, abs=1e-12)
        assert ps.px == pytest.approx(0e0, abs=1e-12)
        assert ps.y == pytest.approx(0e0, abs=1e-12)
        assert ps.py == pytest.approx(0e0, abs=1e-12)
        assert ps.ct == pytest.approx(0e0, abs=1e-12)
        assert ps.delta == pytest.approx(0e0, abs=1e-12)

        # check that all observers reported 0
        for bpm in dc.bpm:
            ps = bpm.get_observer().get_phase_space()
            # Check that the steerer compensates the offset effect
            # of the quadrupole for all bpms ...
            assert ps.x == pytest.approx(0e0, abs=1e-12)
            assert ps.px == pytest.approx(0e0, abs=1e-12)
            assert ps.y == pytest.approx(0e0, abs=1e-12)
            assert ps.py == pytest.approx(0e0, abs=1e-12)
            assert ps.ct == pytest.approx(0e0, abs=1e-12)
            assert ps.delta == pytest.approx(0e0, abs=1e-12)

    # Change quadrupole powering strength ... does not affect
    # the steerer: typically in the percent range
    strength_change = 0.01 + 1
    gradient = quad_muls.get_multipole(2)
    gradient_changed = gradient * strength_change
    quad_muls.set_multipole(2, gradient_changed)
    assert quad_muls.get_multipole(2) == pytest.approx(gradient_changed, 1e-12)

    # Find the fixed point of the orbit ...
    r = find_fixed_point(acc)
    ps = gtpsa.ss_vect_double(r.x)
    print("Distored fixed point\n", r)
    # and check that it is one within acceptable limits ...
    acc.propagate(calc_config, ps)
    assert ps.x == pytest.approx(r.x[0], abs=1e-12)
    assert ps.px == pytest.approx(r.x[1], abs=1e-12)
    assert ps.y == pytest.approx(r.x[2], abs=1e-12)
    assert ps.py == pytest.approx(r.x[3], abs=1e-12)
    assert ps.delta == pytest.approx(r.x[4], abs=1e-12)
    assert ps.ct == pytest.approx(r.x[5], abs=1e-12)

    # one extra propagation ... just in case the observers have not seen it
    acc.propagate(calc_config, ps)

    # get observer values and check to forecast
    dc.acc

    if True:
        # extract displacement from observers
        def extract_orbit_from_observer(elem):
            observer = elem.get_observer()
            ps = observer.get_phase_space()
            return [ps.x, ps.y, observer.get_observed_name(), observer.get_observed_index()]

        observed_displacement = np.array(
            [extract_orbit_from_observer(bpm) for bpm in dc.bpm], dtype=object
        )

        # print(observed_displacement[:, 3])
        s = dc.acc_info.s.sel(index=observed_displacement[:, 3])
        # check that vertical displacements are all close to zero
        assert np.sum(np.absolute(observed_displacement[:, 1])) == pytest.approx(0e0, 1e-8)

    bpm_indices = observed_displacement[:, 3]
    twiss_bpm = twiss.sel(index=bpm_indices)
    t_nu_x_bpm = t_nu_x.sel(index=bpm_indices)
    theta_guess = 1e-5
    co_bpm = closed_orbit_distortion(
        twiss_bpm.twiss.sel(par="beta", plane="x").values,
        t_nu_x_bpm.values * pi2,
        tune=float(t_nu_x.isel(index=-1)),
        theta_i=theta_guess,
        beta_i=twiss.twiss.sel(par="beta", plane="x", index=dc.quad.index).values,
        mu_i=t_nu_x.sel(index=dc.quad.index).values * pi2,
    )

    fig, axes = plt.subplots(2, 1)
    ax_nu, ax_s = axes

    nu_x = [t_nu_x.sel(index=idx) for idx in observed_displacement[:, 3]]

    gradient = quad.get_main_multipole_strength().real
    txt = f"""
    Gradient         {gradient:6.3f} T/m
    offset           {offset*1e6:6.3f} um
    strength change  {strength_change:6.3f}
    """
    print(txt)

    # dipole strength due to
    added_kick = dipole_due_to_quadrupole_gradient_change(
        dc.quad, strength_change - 1, offset
    )
    # Here only interested in x plane
    added_kick = added_kick.real
    added_kick /= theta_guess

    ax_nu.plot(nu_x, observed_displacement[:, 0] * 1e6, "bx")
    ax_nu.set_ylabel("x [um]")
    ax_s.plot(s, observed_displacement[:, 0] * 1e6, "bx")
    ax_s.set_ylabel("x [um]")
    ax_nu.plot(t_nu_x, co * 1e6 * added_kick, "c.-")
    ax_nu.plot(t_nu_x_bpm, co_bpm * 1e6 * added_kick, "c+")
    ax_s.plot(twiss.s, co * 1e6 * added_kick, "c.-")
    ax_s.plot(twiss_bpm.s, co_bpm * 1e6 * added_kick, "c+")


    ax_nu.plot(t_nu_x, forecast_displacement[:, 0] * offset.real * 1e4 /2, "g.-")
    ax_s.plot(twiss.s, forecast_displacement[:, 0] * offset.real * 1e4 / 2, "g.-")
    ax_nu.plot(t_nu_x, forecast_displacement2[:, 0] * offset.real * 1e4 / 2, "m.-")
    ax_s.plot(forecast_displacement2_s, forecast_displacement2[:, 0] * offset.real * 1e4 / 2, "m.-")
    # ax_nu.plot(nu_x, displaced_quad_effect.sel(dep="dx2", plane="x").values * 1, "r.")
    # ax_nu.plot(nu_x, displaced_quad_effect.sel(dep="dKdx", plane="x").values * 1 * 4, "m.")

    # ax_s.plot(s, displaced_quad_effect.sel(dep="dx", plane="x").values * 1, "g.")
    # ax_s.plot(s, displaced_quad_effect.sel(dep="dx2", plane="x").values * 1, "r.")
    # ax.plot(t_nu_x, twiss.twiss.sel(par="beta", plane="x") * 5, "-")
    fig.savefig("displaced_quad_effect_x.pdf")
    fig.savefig("displaced_quad_effect_x.png", dpi=600)
    return

    for x, d in zip(observed_displacement[:, 0], co_bpm):
        print(d / x)
        # predicted = displaced_quad_effect.sel(plane="x", bpm=name).values
        # print(name, index, predicted * offset / x)
    return

    for x, y, name, index in observed_displacement:
        predicted = displaced_quad_effect.sel(plane="x", bpm=name).values
        observed = x
        assert observed == pytest.approx(predicted * offset.real, abs=1e-12)


def test02_quadrupole_not_moved():
    """Check that test setup of first test is only local"""
    dc = create_accelerator()
    assert dc.quad.get_dx() == pytest.approx(0e0, abs=1e-12)


def bar():
    # get steerer reponse .... required to optimise beam
    response = xr.open_dataset("steerer_response_matrix.nc")
    corr_mat_x = response.effect.sel(plane="x", dep="K")
    print(corr_mat_x.shape)


def foo():
    ps = gtpsa.ss_vect_double(0e0)
    for cnt in range(1):
        find_fixed_point(acc)
        orbit = [ob.get_phase_space() for ob in observers]
        plt.plot(ob_info.s, [o.x for o in orbit], ".")
        orbit_bpm_x = [ob.get_phase_space().x for ob in bpm_observers]
        ridge = Ridge(alpha=0.5)
        ridge.fit(corr_mat_x.T, orbit_bpm_x)
        print(np.absolute(ridge.coef_).mean())
    pass


if __name__ == "__main__":
    test01_quadrupole_moved_steerer_compensation()
    plt.show()
