"""Use Case:
     Set the linear chromaticity.
"""


import os
import enum

import numpy as np
from scipy import linalg as la

import gtpsa

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, knobs
from thor_scsi.utils.output import mat2txt, vec2txt


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole  = 3


def compute_map(lat_prop):
    M = lo.compute_map(
        lat_prop._lattice, lat_prop._model_state, desc=lat_prop._desc,
        tpsa_order=lat_prop._no)
    return M


def prt_map(map, str, *, eps: float=1e-30):
    print(str)
    map.x.print("x", eps)
    map.px.print("p_x", eps)
    map.y.print("y", eps)
    map.py.print("p_y", eps)
    map.delta.print("delta", eps)
    map.ct.print("ct", eps)


corresponding_types = {
    ts.Bending:           ts.BendingTpsa,
    ts.Sextupole:         ts.SextupoleTpsa,
    ts.Quadrupole:        ts.QuadrupoleTpsa,
    ts.HorizontalSteerer: ts.HorizontalSteererTpsa,
    ts.VerticalSteerer:   ts.VerticalSteererTpsa,
}


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
def compute_map_tps(lat_prop, lat):
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


def compute_nu_xi(lat_prop, M):
    # nu = acos( Tr{ M_x,y(delta; b_3) } / 2 ) / 2 pi
    planes = ["x", "y"]

    m_11 = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
            gtpsa.tpsa(lat_prop._desc, lat_prop._no)]
    m_22 = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
            gtpsa.tpsa(lat_prop._desc, lat_prop._no)]

    m_11[0].set(0e0, M.x.get([1, 0, 0, 0, 0, 0, 0]))
    m_11[0].set([0, 0, 0, 0, 1, 0, 0], 0e0, M.x.get([1, 0, 0, 0, 1, 0, 0]))

    m_22[0].set(0e0, M.px.get([0, 1, 0, 0, 0, 0, 0]))
    m_22[0].set([0, 0, 0, 0, 1, 0, 0], 0e0, M.px.get([0, 1, 0, 0, 1, 0, 0]))

    m_11[1].set(0e0, M.y.get([0, 0, 1, 0, 0, 0, 0]))
    m_11[1].set([0, 0, 0, 0, 1, 0, 0], 0e0, M.y.get([0, 0, 1, 0, 1, 0, 0]))

    m_22[1].set(0e0, M.py.get([0, 0, 0, 1, 0, 0, 0]))
    m_22[1].set([0, 0, 0, 0, 1, 0, 0], 0e0, M.py.get([0, 0, 0, 1, 1, 0, 0]))

    print()
    tr = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    xi = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
          gtpsa.tpsa(lat_prop._desc, lat_prop._no)]
    for k in range(2):
        tr = m_11[k] + m_22[k]
        if tr.get() < 2e0:
            xi[k] = lo.acos2_tpsa(M.jacobian()[2*k][2*k+1], tr/2e0)/(2e0*np.pi)
        else:
            xi[k] = np.nan
            print("\ncompute_nu_xi: unstable in plane {:s}".format(planes[k]))

    print("xi_x = {:12.5e} {:12.5e}".
          format(xi[ind.X].get(), xi[ind.X].get([0, 0, 0, 0, 1, 0, 0])))
    print("xi_y = {:12.5e} {:12.5e}".
          format(xi[ind.Y].get(), xi[ind.Y].get([0, 0, 0, 0, 1, 0, 0])))

    return xi


def compute_xi_prm(lat_prop, M):
    # nu = acos( Tr{ M_x,y(delta; b_3) } / 2 ) / 2 pi
    planes = ["x", "y"]

    m_11 = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
            gtpsa.tpsa(lat_prop._desc, lat_prop._no)]
    m_22 = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
            gtpsa.tpsa(lat_prop._desc, lat_prop._no)]

    m_11[0].set(0e0, M.x.get([1, 0, 0, 0, 0, 0, 0]))
    m_11[0].set([0, 0, 0, 0, 1, 0, 0], 0e0, M.x.get([1, 0, 0, 0, 1, 0, 0]))
    m_11[0].set([0, 0, 0, 0, 1, 0, 1], 0e0, M.x.get([1, 0, 0, 0, 1, 0, 1]))

    m_22[0].set(0e0, M.px.get([0, 1, 0, 0, 0, 0, 0]))
    m_22[0].set([0, 0, 0, 0, 1, 0, 0], 0e0, M.px.get([0, 1, 0, 0, 1, 0, 0]))
    m_22[0].set([0, 0, 0, 0, 1, 0, 1], 0e0, M.px.get([0, 1, 0, 0, 1, 0, 1]))

    m_11[1].set(0e0, M.y.get([0, 0, 1, 0, 0, 0, 0]))
    m_11[1].set([0, 0, 0, 0, 1, 0, 0], 0e0, M.y.get([0, 0, 1, 0, 1, 0, 0]))
    m_11[1].set([0, 0, 0, 0, 1, 0, 1], 0e0, M.y.get([0, 0, 1, 0, 1, 0, 1]))

    m_22[1].set(0e0, M.py.get([0, 0, 0, 1, 0, 0, 0]))
    m_22[1].set([0, 0, 0, 0, 1, 0, 0], 0e0, M.py.get([0, 0, 0, 1, 1, 0, 0]))
    m_22[1].set([0, 0, 0, 0, 1, 0, 1], 0e0, M.py.get([0, 0, 0, 1, 1, 0, 1]))

    print()
    tr = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    xi = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
          gtpsa.tpsa(lat_prop._desc, lat_prop._no)]
    for k in range(2):
        tr = m_11[k] + m_22[k]
        if tr.get() < 2e0:
            xi[k] = lo.acos2_tpsa(M.jacobian()[2*k][2*k+1], tr/2e0)/(2e0*np.pi)
        else:
            xi[k] = np.nan
            print("\ncompute_xi_prm: unstable in plane {:s}".format(planes[k]))

    print("xi_x = {:12.5e} {:12.5e} {:12.5e}".
          format(xi[ind.X].get(), xi[ind.X].get([0, 0, 0, 0, 1, 0, 0]),
                 xi[ind.X].get([0, 0, 0, 0, 1, 0, 1])))
    print("xi_y = {:12.5e} {:12.5e} {:12.5e}".
          format(xi[ind.Y].get(), xi[ind.Y].get([0, 0, 0, 0, 1, 0, 0]),
                 xi[ind.Y].get([0, 0, 0, 0, 1, 0, 1])))

    return xi


def compute_sext_resp_mat_num(lat_prop, b3_list):
    db_3xL = 1e0
    n = len(b3_list)
    A = np.zeros((n, n))
    M = compute_map(lat_prop)
    stable, _, xi = lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
    for k in range(n):
        b_3xL = \
            lat_prop.get_b_nxL_elem(b3_list[k], 0, MultipoleIndex.sextupole)
        lat_prop.set_b_nxL_fam(
            b3_list[k], MultipoleIndex.sextupole, b_3xL-db_3xL)
        M = compute_map(lat_prop)
        stable, _, xi_1 = lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
        lat_prop.set_b_nxL_fam(
            b3_list[k], MultipoleIndex.sextupole, b_3xL+db_3xL)
        M = compute_map(lat_prop)
        stable, _, xi_2 = lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
        a_ij = (xi_2-xi_1)/(2e0*db_3xL)
        A[ind.X, k] = a_ij[ind.X]
        A[ind.Y, k] = a_ij[ind.Y]
    return xi, A


def compute_sext_resp_mat(lat_prop, b3_list):
    n = len(b3_list)
    xi = np.zeros(2)
    A = np.zeros((n, n))
    for k in range(n):
        lat_ptc = param_dep(lat_prop, b3_list[k])
        M = compute_map_tps(lat_prop, lat_ptc)
        xi_tps = compute_xi_prm(lat_prop, M)
        for j in range(2):
            if k == 0:
                xi[j] = xi_tps[j].get([0, 0, 0, 0, 1, 0, 0])
            A[j, k] = xi_tps[j].get([0, 0, 0, 0, 1, 0, 1])
    return xi, A


def set_xi(lat_prop, xi_x, xi_y, b3_list):
    n = len(b3_list)
    xi, A = compute_sext_resp_mat(lat_prop, b3_list)
    A_inv = la.pinv(A)
    b_3 = A_inv @ (-xi)
    for k in range(len(b3_list)):
        lat_prop.set_b_n_fam(b3_list[k], MultipoleIndex.sextupole, b_3[k])

    print("\nxi:\n", xi)
    print("\nA:\n", A)
    print("\nA_inv:\n", A_inv)
    print("\nb_3:\n", b_3)

    lat_ptc = param_dep(lat_prop, b3_list[0])
    M = compute_map_tps(lat_prop, lat_ptc)
    nu, xi_corr = compute_nu_xi(lat_prop, M)


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

np.set_printoptions(formatter={"float": "{:13.5e}".format})

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

b3_list = ["sfoh", "sdqd"]

if False:
    zero_sext(lat_prop, b3_list)

set_xi(lat_prop, 0e0, 0e0, b3_list)
