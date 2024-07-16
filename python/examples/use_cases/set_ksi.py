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
    print("\nmult_prm ->")
    for k in range(len(mult_family)):
        index = mult_family[k].index
        elem = convert_magnet_to_knobbable(elems[index])
        elems[index] = \
            knobs.make_magnet_knobbable(
                elem, po=1, desc=lat_prop._desc,
                named_index=lat_prop._named_index, multipole_number=mpole_n,
                offset=False
            )
        print("\nmult_prm:", elems[index])
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
    print("<- mult_prm")
    return elems


def lat_mult_prm(lat_prop, mult_prm_name, mpole_n):
    print(lat_prop)
    elems = [elem for elem in lat_prop._lattice]
    mult_family = lat_prop._lattice.elements_with_name(mult_prm_name)
    elems = mult_prm(lat_prop, elems, mult_family, mpole_n)
    return ts.AcceleratorTpsa(elems)


# Work-around for C++ virtual function -> Python mapping issue.
# See function above:
#   mult_prm
def propagate(lat_ptc, lat):
    M = gtpsa.ss_vect_tpsa(
        lat_ptc._desc, lat_ptc._no, lat_ptc._nv,
        index_mapping=lat_ptc._named_index)
    M.set_identity()
    for k in range(len(lat)):
        lat[k].propagate(lat_ptc._model_state, M)
    return M


def zero_sext(lat_prop, b3_list):
    # Zero sextupoles.
    print("\nZeroing sextupoles.")
    for b3_name in b3_list:
        lat_prop.set_b_n_fam(b3_name, MultipoleIndex.sextupole, 0e0)


def param_dep(lat_prop):
    b3_list = ["sf_h", "sd1", "sf_e_h", "sd2"]
    zero_sext(lat_prop, b3_list)

    if False:
        lat_prop._named_index = gtpsa.IndexMapping(
            dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
        desc = gtpsa.desc(nv, no, nv_prm, no_prm)

        print("\n  C [m]     = {:9.7f}".format(C))
        print("  phi [deg] = {:5.3f}".format(compute_phi(lat)))
        print("  nu        = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))
        print("  xi        = [{:7.5f}, {:7.5f}]".format(xi[X_], xi[Y_]))

    lat_prop._named_index = gtpsa.IndexMapping(
        dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6))

    if False:
        nv_prm = 0
        no_prm = 0
        lat_prop._desc = gtpsa.desc(lat_prop._nv, lat_prop._no, nv_prm, no_prm)
        lat_ptc = lat_mult_prm("", lat_prop._lattice, 0)
    else:
        nv_prm = 1
        no_prm = no
        print(no, nv, nv_prm, no_prm)
        lat_prop._desc = gtpsa.desc(lat_prop._nv, lat_prop._no, nv_prm, no_prm)
        print(lat_prop._desc)
        if True:
            lat_ptc = lat_mult_prm(lat_prop, b3_list[0], 3)
        else:
            lat_ptc = lat_mult_prm(lat_prop, "uq3", 2)

    M = propagate(lat_prop, lat_ptc)

    print("\nM:", M)
    if not False:
        prt_map(M, "\nM:")

    if False:
        M2 = M
        # Busted.
        M2.getOrder(M2, 2)
        print("\n:\n", M2)

    if False:
        M_inv = ts.inv(M)

        M_M_inv = gtpsa.ss_vect_tpsa(
            desc, no, nv, index_mapping=lat_prop_named_index)
        M_M_inv.rcompose(M, M_inv)
        print("\nM:", M)
        print("\nM^-1:", M_inv)
        print("\nM*M^-1:", M_M_inv)
        print("\nM.x:\n", M.x)
        print("\nM_M_inv.x:\n", M_M_inv.x)
        # print("\nM*M^-1:", M_M_inv[0])
        assert False


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
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_4u")
lat_name = "max_4u_g_2"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

param_dep(lat_prop)
