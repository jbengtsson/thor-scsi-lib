"""Use Case:
     Set the linear chromaticity.
"""


import os
import enum

import numpy as np

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, prm_class as pc
from thor_scsi.utils.output import mat2txt, vec2txt

import gtpsa

import thor_scsi.lib as ts

from thor_scsi.utils import knobs


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole  = 3


corresponding_types = {
    ts.Sextupole:         ts.SextupoleTpsa,
    ts.Quadrupole:        ts.QuadrupoleTpsa,
    ts.HorizontalSteerer: ts.HorizontalSteererTpsa,
    ts.VerticalSteerer:   ts.VerticalSteererTpsa,
}


def compute_nu_xi(lat_prop, M):
    planes = ["x", "y"]
    try:
        nu, xi = [np.zeros(2), np.zeros(2)]
        M_delta = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
        for k in range(2):
            tr_k = np.trace(M.jacobian()[2*k:2*(k+1)])
            if tr_k >= 2e0:
                print("\ncompute_nu_xi: unstable in plane {:s}".
                      format(planes[k]))
                raise ValueError
            M_delta.clear()
            if k == 0:
                M_delta += M.x.get(lo.ind_1(2*k))
                M_delta.set(
                    lo.ind_1(ind.delta), 0e0, M.x.get(lo.ind_2(2*k, ind.delta)))
                # m_22 + delta * m_26.
                M_delta += M.px.get(lo.ind_1(2*k+1))
                M_delta.set(
                    lo.ind_1(ind.delta), 1e0,
                    M.px.get(lo.ind_2(2*k+1, ind.delta)))
            elif k == 1:
                M_delta += M.y.get(lo.ind_1(2*k))
                M_delta.set(
                    lo.ind_1(ind.delta), 0e0,
                    M.y.get(lo.ind_2(2*k, ind.delta)))
                # m_22 + delta * m_26.
                M_delta += M.py.get(lo.ind_1(2*k+1))
                M_delta.set(
                    lo.ind_1(ind.delta), 1e0,
                    M.py.get(lo.ind_2(2*k+1, ind.delta)))
            # M_delta = m_11 + m_22.
            nu_tpsa = \
                lo.acos2_tpsa(M.jacobian()[2*k][2*k+1], M_delta/2e0)/(2e0*np.pi)
            nu[k], xi[k] = [nu_tpsa.get(), nu_tpsa.get(lo.ind_1(ind.delta))]
    except ValueError:
        print("\ncompute_nu_xi: unstable =", stable)
        return False, np.nan, np.nan
    else:
        return True, nu, xi


def prt_map(map, str, *, eps: float=1e-30):
    print(str)
    map.x.print("x", eps)
    map.px.print("p_x", eps)
    map.y.print("y", eps)
    map.py.print("p_y", eps)
    map.delta.print("delta", eps)
    map.ct.print("ct", eps)


def convert_magnet_to_knobbable(a_magnet: ts.Mpole) -> ts.MpoleTpsa:
    config = a_magnet.config()
    corresponding_type = corresponding_types[type(a_magnet)]
    return corresponding_type(config)


def mult_prm(lat_prop, elems, mult_family, mpole_n):
    for k in range(len(mult_family)):
        index = mult_family[k].index
        elem = convert_magnet_to_knobbable(elems[index])
        elems[index] = \
            knobs.make_magnet_knobbable(
                elem, po=1, desc=lat_prop._desc,
                named_index=lat_prop._named_index, multipole_number=mpole_n,
                offset=False
            )
        # While the RHS pointer can be recasted to:
        #   CellVoid
        #   ElemTypeKnobbed
        #   QuadrupoleType
        # the assignment of the LHS only gives:
        #   CellVoid
        # and the lack of:
        #   ElemtypeKnobbed
        # will lead to an assert on line 158 in:
        #   thor_scsi/std_machine/accelerator.cc
        #
    return elems


def lat_mult_prm(lat_prop, mult_prm_name, mpole_n):
    elems = [elem for elem in lat_prop._lattice]
    mult_family = lat_prop._lattice.elements_with_name(mult_prm_name)
    elems = mult_prm(lat_prop, elems, mult_family, mpole_n)
    return ts.AcceleratorTpsa(elems)


# Work-around for C++ virtual function -> Python mapping issue.
# See function above:
#   mult_prm
def compute_map(lat_prop, lat):
    M = gtpsa.ss_vect_tpsa(
        lat_prop._desc, lat_prop._no, lat_prop._nv,
        index_mapping=lat_prop._named_index)
    M.set_identity()
    for k in range(len(lat)):
        lat[k].propagate(lat_prop._model_state, M)
    return M


def zero_sext(lat_prop, b3_list):
    # Zero sextupoles.
    print("\nZeroing sextupoles.")
    for b3_name in b3_list:
        lat_prop.set_b_n_fam(b3_name, MultipoleIndex.sextupole, 0e0)


def param_dep(lat_prop, prm_name):
    lat_prop._named_index = gtpsa.IndexMapping(
        dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6))

    nv_prm = 1
    no_prm = no
    lat_prop._desc = gtpsa.desc(lat_prop._nv, lat_prop._no, nv_prm, no_prm)
    lat_ptc = lat_mult_prm(lat_prop, prm_name, 3)
    return lat_ptc


ind = ind.index_class()
# Number of phase-space coordinates.
nv = 6
# Variables max order.
no = 3
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_iv")
lat_name = "ake_2"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

b3_list = ["sfoh", "sfmh", "sd"]

zero_sext(lat_prop, b3_list)

lat_ptc = param_dep(lat_prop, b3_list[0])

M = compute_map(lat_prop, lat_ptc)

if False:
    prt_map(M, "\nM:")

stable, nu, xi = compute_nu_xi(lat_prop, M)

print("\nnu = [{:7.5f}, {:7.5f}]".format(nu[ind.X], nu[ind.Y]))
print("xi = [{:7.5f}, {:7.5f}]".format(xi[ind.X], xi[ind.Y]))
