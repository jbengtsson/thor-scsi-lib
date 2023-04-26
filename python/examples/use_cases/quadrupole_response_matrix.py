# Find the lattice file
from pathlib import Path

import tqdm

from thor_scsi.factory import accelerator_from_config
import thor_scsi.lib as tslib
from thor_scsi.lib import Quadrupole, QuadrupoleTpsa, AcceleratorTpsa, ConfigType
from thor_scsi.utils.accelerator import instrument_with_standard_observers
import gtpsa.utils
import numpy as np
import pandas as pd
import copy
import os
import sys
import xarray as xr
from typing import Sequence
from dataclasses import dataclass

named_index_d = dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6, dx=7, dy=8)
named_index = gtpsa.IndexMapping(named_index_d)
# default_desc = gtpsa.desc(nv=6+3, mo=3)
# maximum order of variables
mo = 3
# maximum order of parameters ... at least one
po = 2
default_desc = gtpsa.desc(nv=6, mo=mo, np=3, po=po)
# default_desc = gtpsa.desc(nv=6, mo=3, np=144*3, po=1)

t_file = (
    Path(os.environ["HOME"])
    / "Devel"
    / "gitlab"
    / "dt4acc"
    / "lattices"
    / "b2_stduser_beamports_blm_tracy_corr.lat"
)


def convert_quad(a_quad: Quadrupole):
    config = a_quad.config()
    quad_tpsa = QuadrupoleTpsa(config)
    return quad_tpsa


def lattice_tpsa(acc):
    def convert_if_quad(elem):
        if isinstance(elem, Quadrupole):
            n_elem = convert_quad(elem)
        else:
            n_elem = elem
        return n_elem

    elements = [convert_if_quad(elem) for elem in acc]
    return AcceleratorTpsa(elements)


def make_quad_knobbable(quad, *, desc=default_desc, named_index=named_index):
    # Augment the quadrupole
    k_orig = quad.get_main_multipole_strength().to_object()
    k = gtpsa.ctpsa(desc, po, mapping=named_index)
    dx = gtpsa.tpsa(desc, po, mapping=named_index)
    dy = gtpsa.tpsa(desc, po, mapping=named_index)
    k.name = quad.name + "_K"
    dx.name = quad.name + "_dx"
    dx.name = quad.name + "_dy"
    using_knobs = True
    dxv = 1e-6
    dyv = 1e-6
    dxv = 0e0
    dyv = 0e0
    if using_knobs:
        k.set_knob(k_orig, "K")
        dx.set_knob(dxv, "dx")
        dy.set_knob(dyv, "dy")
    else:
        k.set_variable(k_orig, "K")
        dx.set_variable(dxv, "dx")
        dy.set_variable(dyv, "dy")
    # European Norm ... thus 2
    quad.get_multipoles().set_multipole(2, gtpsa.CTpsaOrComplex(k))
    quad.set_dx(gtpsa.TpsaOrDouble(dx))
    quad.set_dy(gtpsa.TpsaOrDouble(dy))
    return quad


def make_quad_unknobbable(quad):
    """Replace knobbed variables with standard values"""
    dx = float(quad.get_dx().to_object().get())
    dy = float(quad.get_dy().to_object().get())
    quad.set_dx(gtpsa.TpsaOrDouble(dx))
    quad.set_dy(gtpsa.TpsaOrDouble(dy))
    mul = complex(quad.get_multipoles().get_multipole(2).to_object().get())
    quad.get_multipoles().set_multipole(2, gtpsa.CTpsaOrComplex(mul))
    return quad



def response_of_one_quad(
    acc, quad, *, desc=default_desc, n_start=0, n_elements=2 ** 30, **kwargs
):
    quad = make_quad_knobbable(quad, desc=desc, **kwargs)
    ps = gtpsa.ss_vect_tpsa(desc, mo, 6, index_mapping=named_index)
    ps.set_identity()
    calc_config = ConfigType()
    acc.propagate(calc_config, ps, n_start, n_elements)
    make_quad_unknobbable(quad)
    return ps


print(f"Reading lattice file {t_file}")
acc = accelerator_from_config(t_file)
acc_tpsa = lattice_tpsa(acc)

# assuming that the fixed point corresponds to 0, 0,
# add observers before we propagate
# select bpms
n_start = 0
n_elements = 2000

compare_to_lattice_without_knobs = False
if compare_to_lattice_without_knobs:
    desc = default_desc
    ps = gtpsa.ss_vect_tpsa(desc, mo, 6, named_index)
    ps.set_identity()
    calc_config = ConfigType()
    # print(ps)
    acc.propagate(calc_config, ps, n_start, n_elements)
    # print(ps)
    ps2 = copy.copy(ps)

    print("tpsa")
    ps.set_identity()
    # print(ps)
    acc_tpsa.propagate(calc_config, ps, n_start, n_elements)
    # print(ps)

    print("Difference ")
    print(ps2 - ps)


for elem in acc_tpsa:
    if isinstance(elem, QuadrupoleTpsa):
        print(elem.name, elem.index)
        ps3 = response_of_one_quad(
            acc_tpsa, elem, n_start=n_start, n_elements=n_elements
        )
        break

print("#----------------------------------------------------------------------")
print("Replaced quad")
print("#----------------------------------------------------------------------")

if compare_to_lattice_without_knobs:
    np.set_printoptions(precision=3)
    print(ps3)
    print("Find difference")
    print("cst", ps3.cst() - ps.cst())
    print("ps\n", ps)
    print("ps3\n", ps3)
    jac = ps.jacobian()
    jac3 = ps3.jacobian()
    print("ps jac\n", jac)
    print("ps3 jac\n\n", jac3)
    n_jac = np.array(jac)
    n_jac3 = np.array(jac3)
    print("ps jac\n", n_jac)
    print("ps3 jac\n\n", n_jac3)
    nrows, ncols = jac.shape
    jac3 = jac3[:nrows, :ncols]
    djac = jac3 - jac
    # djac =  (djac*1e6).astype(int) / 1e6
    print("djac ps\n", djac[:6, :6])
    print("djac knobs?\n", djac[:, 6:])
    print("djac knobs?\n", djac[6:, :])
    hes = ps.hessian()
    hes = ps3.hessian()
    nslices, nrows, ncols = hes.shape
    # print(ps3.hessian()[:nslices, :nrows, :ncols] - hes)


bpm_elems = [elem for elem in acc if isinstance(elem, tslib.BPM)]
observers = instrument_with_standard_observers(bpm_elems, mapping=named_index)


def combine_responses(observers: Sequence[tslib.StandardObserver]):
    def extract_x_y(ob):
        ps = ob.get_truncated_power_series_a()
        data = [
            (ps.x.get(d), ps.y.get(d))
            for d in (dict(K=1), dict(K=2), dict(dx=1), dict(dy=1))
        ]

        ob_name = ob.get_observed_name()
        da = xr.DataArray(
            data=[data],
            dims=["bpm", "dep", "plane"],
            coords=[
                [ob_name],
                ["K", "K2", "dx", "dy"],
                ["x", "y"],
            ],
            name="bpm_observation",
        )

        idx = ob.get_observed_index()
        index = xr.DataArray(name="index", data=[idx], dims=["bpm"], coords=[[ob_name]])
        ds = xr.Dataset(dict(index=index, effect=da))
        return ds

    data = [extract_x_y(ob) for ob in observers]
    ds = xr.concat(data, dim="bpm")
    return ds


def abs(v: gtpsa.tpsa):
    if v.get() < 0e0:
        return v * -1
    return v


def acos2(sin, cos):
    """arcos(phi): 0 - 2 * pi.
    The sin part is used to determine the quadrant.
    """
    if abs(cos).get() > 1e0:
        raise ValueError(f"arg = {cos:22.15e} > 1e0")
    phi = gtpsa.acos(cos)
    if sin.get() < 0e0:
        phi = 2e0 * np.pi - phi
    return phi


def compute_phase_advance(
    ps: gtpsa.ss_vect_tpsa, canonical_coordinates_names: Sequence[str]
) -> gtpsa.tpsa:
    """compute phase advance

    Args:
        ps: phase space
        canonical_coordinates_names: name of the two coordinates

    Returns:
        phase advance as gtpsa.tpsa object
    """
    # names of the 2 coordinates
    c, d = canonical_coordinates_names
    # compute the trace
    tr = ps.loc[c].deriv(c) + ps.loc[d].deriv(d)
    tr_cst = tr.get()
    if not tr_cst < 2.0:
        raise ValueError(f"Trace not stable: trace = {tr}")
    if True:
        # a bit a lot of calculation just to find the quadrant
        s2 = ps.loc[c].deriv(d) * ps.loc[d].deriv(c)
        sgn = s2.get() > 0e0 or -1e0
        sin = gtpsa.sqrt(abs(s2))
        phase_advance = acos2(sgn * sin, tr / 2.0) / (2 * np.pi)
    else:
        # ignoring quadrant
        phase_advance = gtpsa.acos(tr / 2.0) / (2 * np.pi)
    return phase_advance


calc_config = ConfigType()
acc_tpsa.propagate(calc_config, ps3)
compute_phase_advance(ps3, ["x", "px"])
response_one_quad = combine_responses(observers)
print(response_one_quad.effect.sel(plane="x", dep="dx"))

all_responses = []


def process_one_quad(quad):
    nps = response_of_one_quad(acc_tpsa, quad)
    acc_tpsa.propagate(calc_config, nps)
    response_one_quad = combine_responses(observers)
    response_one_quad["name"] = quad.name
    # list comprehension as currently iloc not active?
    ps = [
        np.array(
            [
                nps.x.get(d),
                nps.px.get(d),
                nps.y.get(d),
                nps.py.get(d),
                nps.delta.get(d),
                nps.ct.get(d),
            ]
        )
        for d in [dict(K=1), dict(K=2), dict(dx=1), dict(dy=1)]
    ]
    phase_space = xr.DataArray(
        data=ps,
        name="phase_space",
        dims=["dep", "coor"],
        coords=[["K", "K2", "dx", "dy"], ["x", "px", "y", "py", "delta", "ct"]],
    )
    response_one_quad["phase_space"] = phase_space

    pa_x = compute_phase_advance(nps, ["x", "px"])
    pa_y = compute_phase_advance(nps, ["y", "py"])

    data = [
        np.array(
            [
                pa_x.get(d),
                pa_y.get(d),
            ]
        )
        for d in [dict(K=1), dict(K=2), dict(dx=1), dict(dy=1)]
    ]
    phase_advance = xr.DataArray(
        data=data,
        dims=["dep", "plane"],
        coords=[["K", "K2", "dx", "dy"], ["x", "y"]],
    )

    response_one_quad["phase_advance"] = phase_advance

    return response_one_quad


all_responses = {
    elem.name: process_one_quad(elem)
    for elem in tqdm.tqdm(acc_tpsa, total=len(acc_tpsa))
    if isinstance(elem, QuadrupoleTpsa)
}
tmp = [(key, item) for key, item in all_responses.items()]
db = xr.concat([t[1] for t in tmp], pd.Index([t[0] for t in tmp], name="quadrupole"))
db.to_netcdf("quadrupole_response_matrix2.nc")
