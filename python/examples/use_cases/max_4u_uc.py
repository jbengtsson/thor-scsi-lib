"""Use Case:
     Module for implementing & optimising a unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

from dataclasses import dataclass
import os
from typing import Tuple

import numpy as np
from scipy import optimize as opt

import gtpsa

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, prm_class as pc
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()

L_min        = 0.10
L_max        = 0.50
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0

eps_x_des    = 139e-12
nu_uc        = [0.265, 0.0816]


def compute_linear_optics(lat_prop):
    # Compute the beam dynamics properties.
    J_z_min = 1e-3

    try:
        # Compute Twiss parameters along lattice.
        if not lat_prop.comp_per_sol():
            print("\nf_sp - comp_per_sol: unstable")
            raise ValueError
        else:
            _, _, _, nu = lat_prop.get_Twiss(len(lat_prop._lattice)-1)

        # Compute radiation properties.
        stable, stable_rad = lat_prop.compute_radiation()
        if stable and not stable_rad:
            stable = lat_prop._J[ind.Z] > J_z_min
            print("\ncompute_linear_optics: J_z = {:9.3e}".format(
                lat_prop._J[ind.Z]))
        if not stable:
            print("\nf_sp - compute_radiation: unstable")
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


def prt_nl(h_im, K_re):
    print()
    for key in h_dict:
        print("    {:s} = {:10.3e}".format(key, h_im.get(h_dict[key])))
        if key == "h_00201":
            print()
    print()
    for key in K_dict:
        print("    {:s} = {:10.3e}".format(key, K_re.get(K_dict[key])))
        if key == "K_00220":
            print()


def compute_rms(h, dict):
    var = 0e0
    for key in dict:
        var += h.get(dict[key])**2
    return np.sqrt(var)


def opt_uc \
        (lat_prop, prm_list, weight, b1_list, phi_lat, rb_1, phi_b, eps_x_des,
         nu_uc, Id_scl):
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
        print("    h_im rms       =  {:9.3e}".format(h_rms))
        print("    K_re rms       =  {:9.3e}".format(K_rms))
        print("\n    phi_sp         = {:8.5f}".format(phi))
        print("    C [m]          = {:8.5f}".format(lat_prop.compute_circ()))
        print("\n    phi_b1         = {:8.5f}".format(phi_b1))
        print("    phi_rb_1       = {:8.5f}".format(phi_rb_1))
        print("    phi_bend       = {:8.5f}".format(phi_bend))
        prt_nl(h_im, K_re)
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
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_4u")
lat_name = "max_4u_g_1"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.prt_lat(lat_name+"_lat.txt")

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

opt_uc \
    (lat_prop, prm_list, weight, b1_list, phi_lat, rb_1, phi_b, eps_x_des,
     nu_uc, Id_scl)
