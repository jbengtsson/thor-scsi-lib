"""Use Case:
     Module for implementing & optimising a unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

import os
import sys
from dataclasses import dataclass
import enum
from typing import Tuple

import numpy as np
from scipy import optimize as opt
from scipy import linalg as la

import gtpsa

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, knobs as kb, prm_class as pc
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()

L_min        = 0.10
L_max        = 0.50
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0

eps_x_des    = 139e-12
nu_uc        = [0.265, 0.0816]


class MpoleInd(enum.IntEnum):
    quad = 2
    sext  = 3


def prt_map(map, str, *, eps: float=1e-30):
    print(str)
    map.x.print("x", eps)
    map.px.print("p_x", eps)
    map.y.print("y", eps)
    map.py.print("p_y", eps)
    map.delta.print("delta", eps)
    map.ct.print("ct", eps)


def compute_map(lat_prop, no):
    M = lo.compute_map(
        lat_prop._lattice, lat_prop._model_state, desc=lat_prop._desc,
        tpsa_order=no)
    return M


def set_mult_prm(lat_prop, mult_prm_name, mpole_n):
    elems = [elem for elem in lat_prop._lattice]
    mult_family = lat_prop._lattice.elements_with_name(mult_prm_name)

    for k in range(len(mult_family)):
        index = mult_family[k].index
        elem = kb.convert_magnet_to_knobbable(elems[index])
        elems[index] = \
            kb.make_magnet_knobbable(
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


def set_mult_prm(lat_prop, mult_prm_name, mpole_n):
    elems = [elem for elem in lat_prop._lattice]
    mult_family = lat_prop._lattice.elements_with_name(mult_prm_name)

    for k in range(len(mult_family)):
        index = mult_family[k].index
        elem = kb.convert_magnet_to_knobbable(elems[index])
        elems[index] = \
            kb.make_magnet_knobbable(
                elem, po=1, desc=lat_prop._desc,
                named_index=lat_prop._named_index, multipole_number=mpole_n,
                offset=False
            )
    return elems


def unset_mult_prm(lat_ptc, mult_prm_name, mpole_n):
    elems = [elem for elem in lat_prop._lattice]
    mult_family = lat_ptc.elements_with_name(mult_prm_name)
    print(mult_family)
    assert False

    for k in range(len(mult_family)):
        index = mult_family[k].index
        elem = kb.convert_magnet_to_unknobbable(elems[index])
        elems[index] = \
            kb.make_magnet_unknobbable(
                elem, po=1, desc=lat_prop._desc,
                named_index=lat_prop._named_index, multipole_number=mpole_n,
                offset=False
            )
    return elems


# Work-around for C++ virtual function -> Python mapping issue.
# See function above:
#   set_mult_prm
def compute_map_prm(lat_prop, lat):
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
        lat_prop.set_b_n_fam(b3_name, MpoleInd.sext, 0e0)


def set_param_dep(lat_prop, prm_name):
    lat_prop._named_index = gtpsa.IndexMapping(
        dict(x=0, px=1, y=2, py=3, delta=4, ct=5, prm=6, K=7))
    nv_prm = 1
    no_prm = no
    lat_prop._desc = gtpsa.desc(lat_prop._nv, lat_prop._no, nv_prm, no_prm)
    lat_ptc = set_mult_prm(lat_prop, prm_name, 3)
    return lat_ptc


def unset_param_dep(lat_ptc, prm_name):
    lat_ptc = unset_mult_prm(lat_prop, prm_name, 3)

    nv_prm = 0
    lat_prop._named_index = gtpsa.IndexMapping(
        dict(x=0, px=1, y=2, py=3, delta=4, ct=5, prm=6))

    lat_prop._desc = gtpsa.desc(lat_prop._nv, lat_prop._no, nv_prm, no_prm)
    return lat_ptc


def compute_nu_tps_prm(lat_prop, M):
    # nu = acos( Tr{ M_x,y(delta; b_3) } / 2 ) / 2 pi
    planes = ["x", "y"]

    m_11 = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
            gtpsa.tpsa(lat_prop._desc, lat_prop._no)]
    m_22 = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
            gtpsa.tpsa(lat_prop._desc, lat_prop._no)]

    m_11[0].set(0e0, M.x.get([1, 0, 0, 0, 0, 0, 0, 0]))
    m_11[0].set(
        [0, 0, 0, 0, 1, 0, 0, 0], 0e0, M.x.get([1, 0, 0, 0, 1, 0, 0, 0]))
    m_11[0].set(
        [0, 0, 0, 0, 1, 0, 0, 1], 0e0, M.x.get([1, 0, 0, 0, 1, 0, 0, 1]))

    m_22[0].set(0e0, M.px.get([0, 1, 0, 0, 0, 0, 0, 0]))
    m_22[0].set(
        [0, 0, 0, 0, 1, 0, 0, 0], 0e0, M.px.get([0, 1, 0, 0, 1, 0, 0, 0]))
    m_22[0].set(
        [0, 0, 0, 0, 1, 0, 0, 1], 0e0, M.px.get([0, 1, 0, 0, 1, 0, 0, 1]))

    m_11[1].set(0e0, M.y.get([0, 0, 1, 0, 0, 0, 0, 0]))
    m_11[1].set(
        [0, 0, 0, 0, 1, 0, 0, 0], 0e0, M.y.get([0, 0, 1, 0, 1, 0, 0, 0]))
    m_11[1].set(
        [0, 0, 0, 0, 1, 0, 0, 1], 0e0, M.y.get([0, 0, 1, 0, 1, 0, 0, 1]))

    m_22[1].set(0e0, M.py.get([0, 0, 0, 1, 0, 0, 0, 0]))
    m_22[1].set(
        [0, 0, 0, 0, 1, 0, 0, 0], 0e0, M.py.get([0, 0, 0, 1, 1, 0, 0, 0]))
    m_22[1].set(
        [0, 0, 0, 0, 1, 0, 0, 1], 0e0, M.py.get([0, 0, 0, 1, 1, 0, 0, 1]))

    tr = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    xi = [gtpsa.tpsa(lat_prop._desc, lat_prop._no),
          gtpsa.tpsa(lat_prop._desc, lat_prop._no)]
    for k in range(2):
        tr = m_11[k] + m_22[k]
        if tr.get() < 2e0:
            xi[k] = lo.acos2_tpsa(M.jacobian()[2*k][2*k+1], tr/2e0)/(2e0*np.pi)
        else:
            xi[k] = np.nan
            print("\ncompute_nu_tps_prm: unstable in plane {:s}".
                  format(planes[k]))

    return xi


def compute_sext_resp_mat(lat_prop, b3_list):
    n = len(b3_list)
    xi = np.zeros(2)
    A = np.zeros((n, n))
    for k in range(n):
        lat_ptc = set_param_dep(lat_prop, b3_list[k])

        M = compute_map_prm(lat_prop, lat_ptc)
        xi_tps = compute_nu_tps_prm(lat_prop, M)
        for j in range(2):
            if k == 0:
                xi[j] = xi_tps[j].get([0, 0, 0, 0, 1, 0, 0, 0])
            A[j, k] = xi_tps[j].get([0, 0, 0, 0, 1, 0, 0, 1])

    return xi, A


def set_xi(lat_prop, xi_x, xi_y, b3_list):
    n = len(b3_list)

    xi, A = compute_sext_resp_mat(lat_prop, b3_list)
    A_inv = la.pinv(A)
    db_3 = A_inv @ (-xi)

    for k in range(len(b3_list)):
        b_3 = lat_prop.get_b_n_elem(b3_list[k], 0, MpoleInd.sext)
        lat_prop.set_b_n_fam(b3_list[k], MpoleInd.sext, b_3+db_3[k])

    M = gtpsa.ss_vect_tpsa(
        lat_prop._desc, lat_prop._no, lat_prop._nv,
        index_mapping=lat_prop._named_index)
    M.set_identity()
    lat_prop._lattice.propagate(lat_prop._model_state, M)
    nu_tps = compute_nu_tps_prm(lat_prop, M)

    print("\nnu_tps_x = {:12.5e} {:12.5e}".
          format(nu_tps[ind.X].get(),
                 nu_tps[ind.X].get([0, 0, 0, 0, 1, 0, 0, 0])))
    print("nu_tps_y = {:12.5e} {:12.5e}".
          format(nu_tps[ind.Y].get(),
                 nu_tps[ind.Y].get([0, 0, 0, 0, 1, 0, 0, 0])))

    lat_prop._named_index = gtpsa.IndexMapping(
        dict(x=0, px=1, y=2, py=3, delta=4, ct=5, prm=6))
    nv_prm = 0
    no_prm = no
    lat_prop._desc = gtpsa.desc(lat_prop._nv, lat_prop._no, nv_prm, no_prm)



def compute_linear_optics(lat_prop):
    # Compute the beam dynamics properties.
    alpha_rad_z_min = -1e-7

    try:
        # Compute Twiss parameters along lattice.
        if not lat_prop.comp_per_sol():
            print("\ncomp_per_sol: unstable")
            raise ValueError
        else:
            _, _, _, nu = lat_prop.get_Twiss(len(lat_prop._lattice)-1)

        # Compute radiation properties.
        stable, stable_rad = lat_prop.compute_radiation()
        if not stable:
            print("\ncompute_radiation: unstable")
            raise ValueError
        if lat_prop._alpha_rad[ind.Z] > alpha_rad_z_min:
            print("\ncompute_linear_optics:\n"+
                  " unstable in the longitudinal plane alpha_rad_z = {:9.3e}".
                  format(lat_prop._alpha_rad[ind.Z]))
            raise ValueError

        # Compute linear chromaticity.
        M = lo.compute_map(
            lat_prop._lattice, lat_prop._model_state,
            desc=lat_prop._desc, tpsa_order=2)
        stable, _, xi = \
            lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
        if not stable:
            print("\nf_sp - compute_nu_xi: unstable")
            raise ValueError
    except ValueError:
        return False, np.nan, np.nan
    else:
        return True, nu, xi


def compute_twoJ(A_max, beta_inj):
    twoJ = \
        np.array(
            [A_max[ind.X]**2/beta_inj[ind.X], A_max[ind.Y]**2/beta_inj[ind.Y]])
    return twoJ


def compute_Id_scl(lat_prop, twoJ):
    Id_scl = \
        gtpsa.ss_vect_tpsa(
            lat_prop._desc, lat_prop._no, index_mapping=lat_prop._named_index)
    Id_scl.set_identity()
    for k in range(4):
        Id_scl.iloc[k].set_variable(0e0, k+1, np.sqrt(twoJ[k//2]))
    Id_scl.delta.set_variable(0e0, 5, delta_max)
    return Id_scl


def compose_bs(h, map):
    Id = \
        gtpsa.ss_vect_tpsa(
        lat_prop._desc, lat_prop._no, index_mapping=lat_prop._named_index)
    t_map = \
        gtpsa.ss_vect_tpsa(
        lat_prop._desc, lat_prop._no, index_mapping=lat_prop._named_index)
    t_map.x = h
    t_map.compose(t_map, map)
    return t_map.x 


def compute_h(lat_prop, M):
    h    = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    h_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    h_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)

    M.M_to_h_DF(h)
    h.CtoR(h_re, h_im)
    return h_re, h_im


def compute_map_normal_form(lat_prop, M):
    A_0  = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
    A_1  = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
    R    = gtpsa.ss_vect_tpsa(lat_prop._desc, lat_prop._no)
    g    = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    g_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    g_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    K    = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    K_re = gtpsa.tpsa(lat_prop._desc, lat_prop._no)
    K_im = gtpsa.tpsa(lat_prop._desc, lat_prop._no)

    M.Map_Norm(A_0, A_1, R, g, K)
    K.CtoR(K_re, K_im)
    return A_0, A_1, R, g_re, g_im, K_re, K_im


h_dict = {
    "h_10002" : [1, 0, 0, 0, 2, 0, 0],
    "h_20001" : [2, 0, 0, 0, 1, 0, 0],
    "h_00201" : [0, 0, 2, 0, 1, 0, 0],

    "h_30000" : [3, 0, 0, 0, 0, 0, 0],
    "h_21000" : [2, 1, 0, 0, 0, 0, 0],
    "h_11100" : [1, 1, 1, 0, 0, 0, 0],
    "h_10200" : [1, 0, 2, 0, 0, 0, 0],
    "h_10020" : [1, 0, 0, 2, 0, 0, 0]
}

K_dict = {
    "K_22000" : [2, 2, 0, 0, 0, 0, 0],
    "K_11110" : [1, 1, 1, 1, 0, 0, 0],
    "K_00220" : [0, 0, 2, 2, 0, 0, 0],

    "K_11002" : [1, 1, 0, 0, 2, 0, 0],
    "K_00112" : [0, 0, 1, 1, 2, 0, 0]
}


def compute_rms(h, dict):
    var = 0e0
    for key in dict:
        var += h.get(dict[key])**2
    return np.sqrt(var)


def prt_nl(h_rms, K_rms, h_im, K_re):
    print("\n    h_im rms   = {:9.3e}".format(h_rms))
    print("    K_re rms   = {:9.3e}".format(K_rms))

    print()
    for key in h_dict:
        print("    {:s} = {:10.3e}".format(key, h_im.get(h_dict[key])))
        if key == "h_00201":
            print()
    print("\n    K_11001 = {:9.3e}".format(
        K_re.get([1, 1, 0, 0, 1, 0, 0]), K_re.get([0, 0, 1, 1, 1, 0, 0])))
    print("    K_00111 = {:9.3e}".format(
        K_re.get([1, 1, 0, 0, 1, 0, 0]), K_re.get([0, 0, 1, 1, 1, 0, 0])))
    print()
    for key in K_dict:
        print("    {:s} = {:10.3e}".format(key, K_re.get(K_dict[key])))
        if key == "K_00220":
            print()


def opt_uc \
        (lat_prop, prm_list, weight, b1_list, phi_lat, rb_1, phi_b, eps_x_des,
         nu_uc, Id_scl, b3_list):
    """Use Case: optimise super period.
    """
    chi_2_min = 1e30
    eta       = np.nan
    n_iter    = 0
    file_name = "opt_uc.txt"

    def prt_iter(prm, chi_2, nu, xi, h_im, K_re, h_rms, K_rms):
        nonlocal n_iter

        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        phi = lat_prop.compute_phi_lat()
        phi_b1 = compute_phi_bend(lat_prop, b1_list)
        phi_rb_1 = lat_prop.get_phi_elem(rb_1, 0)
        phi_bend = lat_prop.get_phi_elem(phi_b, 0)

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("    eps_x [pm.rad] = {:5.3f} [{:5.3f}]".
              format(1e12*lat_prop._eps[ind.X], 1e12*eps_x_des))
        print("    nu_uc          =  [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu[ind.X], nu[ind.Y], nu_uc[ind.X], nu_uc[ind.Y]))
        print("    xi             =  [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        print("\n    phi_sp         = {:8.5f}".format(phi))
        print("    C [m]          = {:8.5f}".format(lat_prop.compute_circ()))
        print("\n    phi_b1         = {:8.5f}".format(phi_b1))
        print("    phi_rb_1       = {:8.5f}".format(phi_rb_1))
        print("    phi_bend       = {:8.5f}".format(phi_bend))
        prt_nl(h_rms, K_rms, h_im, K_re)
        lat_prop.prt_rad()
        prm_list.prt_prm(prm)

    def compute_chi_2(nu, xi, h_rms, K_rms):
        prt = not False

        dchi_2 = weight[0]*(lat_prop._eps[ind.X]-eps_x_des)**2
        chi_2 = dchi_2
        if prt:
            print("\n  dchi2(eps_x)    = {:10.3e}".format(dchi_2))

        dchi_2 = weight[1]*(nu[ind.X]-nu_uc[ind.X])**2            
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_x)  = {:10.3e}".format(dchi_2))

        dchi_2 = weight[2]*(nu[ind.Y]-nu_uc[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_y)  = {:10.3e}".format(dchi_2))

        dchi_2 = weight[3]*(xi[ind.X]**2+xi[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(xi)       = {:10.3e}".format(dchi_2))

        dchi_2 = weight[4]*h_rms**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(h_rms)    = {:10.3e}".format(dchi_2))

        dchi_2 = weight[5]*K_rms**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(K_rms)    = {:10.3e}".format(dchi_2))

        return chi_2

    def f_sp(prm):
        nonlocal chi_2_min, n_iter

        n_iter += 1
        prm_list.set_prm(prm)
        phi_lat.set_phi_lat()

        stable, nu, xi = compute_linear_optics(lat_prop)

        if stable:
            set_xi(lat_prop, 0e0, 0e0, b3_list)

            M = compute_map(lat_prop, no)

            h_re, h_im = compute_h(lat_prop, M)
            A_0, A_1, R, g_re, g_im, K_re, K_im = \
                compute_map_normal_form(lat_prop, M)

            h_im = compose_bs(h_im, Id_scl)
            K_re = compose_bs(K_re, Id_scl)

            h_rms = compute_rms(h_im, h_dict)
            K_rms = compute_rms(K_re, K_dict)

            chi_2 = compute_chi_2(nu, xi, h_rms, K_rms)
            if chi_2 < chi_2_min:
                prt_iter(prm, chi_2, nu, xi, h_im, K_re, h_rms, K_rms)
                pc.prt_lat(lat_prop, "opt_uc.txt", prm_list, phi_lat=phi_lat)
                chi_2_min = min(chi_2, chi_2_min)
            return chi_2
        else:
            return 1e30


    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4
    g_tol    = 1e-5

    prm, bounds = prm_list.get_prm()

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol, xtol
    # CG     gtol
    minimum = opt.minimize(
        f_sp,
        prm,
        method="CG",
        # callback=prt_iter,
        bounds = bounds,
        # options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
        options={"gtol": g_tol, "maxiter": max_iter}
    )

    print("\n".join(minimum))


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 4
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps   = 1e-15
E_0       = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.prt_lat("max_4u_uc_lat.txt")

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

if not False:
    compute_linear_optics(lat_prop)

    lat_prop.prt_lat_param()
    lat_prop.prt_rad()
    lat_prop.prt_M()
    lat_prop.prt_M_rad()

# Weights.
weight = np.array([
    1e19, # eps_x.
    0e-2, # nu_uc_x.
    0e-2, # nu_uc_y.
    0e-7, # xi.
    1e8,  # h rms.
    1e12  # K rms.
])

b1_list = ["b1_0", "b1_1", "b1_2", "b1_3", "b1_4", "b1_5"]

b1_bend = pc.bend_class(lat_prop, b1_list, phi_max, b_2_max)

if True:
    prms = [
        ("qf1",      "b_2"),

        ("b1_0",     "b_2"),
        ("b1_1",     "b_2"),
        ("b1_2",     "b_2"),
        ("b1_3",     "b_2"),
        ("b1_4",     "b_2"),
        ("b1_5",     "b_2"),

        ("phi_bend", b1_bend),
        ("qf1",     "phi"),
    ]
else:
    prms = [
        ("qf1",      "b_2"),

        ("b1_0",     "b_2"),
        ("b1_1",     "b_2"),
        ("b1_2",     "b_2"),
        ("b1_3",     "b_2"),
        ("b1_4",     "b_2"),
        ("b1_5",     "b_2"),

        ("b1_0",     "phi"),
        ("b1_1",     "phi"),
        ("b1_2",     "phi"),
        ("b1_3",     "phi"),
        ("b1_4",     "phi"),
        # ("b1_5",     "phi"),

        ("qf1",     "phi"),
    ]

prm_list = pc.prm_class(lat_prop, prms, b_2_max)

rb_1    = "qf1"
phi_b   = "b1_5"
n_phi_b = 10
# To maintain the total bend angle.
phi_lat = pc.phi_lat_class(lat_prop, n_phi_b, phi_b)

twoJ = compute_twoJ(A_max, beta_inj)
Id_scl = compute_Id_scl(lat_prop, twoJ)

b3_list = ["sf_h", "sd1"]
# b3_list = ["sfoh", "sdqd"]

opt_uc \
    (lat_prop, prm_list, weight, b1_list, phi_lat, rb_1, phi_b, eps_x_des,
     nu_uc, Id_scl, b3_list)
