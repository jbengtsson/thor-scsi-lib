"""Use Case:
     Module for general linear optics matching to obtain a periodic solution.
"""


import os
import sys
import copy as _copy
from dataclasses import dataclass
from typing import ClassVar
from typing import Tuple

import numpy as np
from scipy import optimize as opt

import gtpsa

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, get_set_mpole as gs, \
    index_class as ind, linear_optics as lo, prm_class as pc, \
    courant_snyder as cs
from thor_scsi.utils.output import mat2txt, vec2txt


ind = ind.index_class()

phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0

eta_des   = [0.0, 0.0]   # [eta_x, eta'_x] at the centre of the straight.
alpha_des = [0.0, 0.0]   # Alpha_x,y at the centre of the straight.
beta_des  = [5.0, 3.0]   # Beta_x,y at the centre of the straight.
dnu_des   = [0.5, 0.25]  # Phase advance across the straight.


@dataclass
class gtpsa_prop:
    # GTPSA properties.
    # Number of variables - phase-space coordinates & 1 for parameter
    #dependence.
    nv: ClassVar[int] = 6 + 1
    # Max order.
    no: ClassVar[int] = 1
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]


def compute_periodic_cell(lat_prop, uc_list):
    M = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    M.set_identity()
    # The 3rd argument is the 1st element index & the 4th the number of elements
    # to propagate across.
    lat_prop._lattice.propagate(
        lat_prop._model_state, M, uc_list[0], uc_list[1]-uc_list[0]+1)
    stable, nu, A, A_inv, _ = lo.compute_M_diag(2, M.jacobian())
    eta, Twiss = lo.transform_matrix_extract_twiss(A)
    alpha = np.array((Twiss[ind.X][0], Twiss[ind.Y][0]))
    beta = np.array((Twiss[ind.X][1], Twiss[ind.Y][1]))
    Twiss = eta, alpha, beta
    print("\ncompute_periodic_cell:")
    lat_prop.prt_Twiss_param(Twiss)
    return Twiss, A


def match_straight(
        lat_prop, prm_list, uc_list, sp_list, weight, d1_bend, d2_bend,
        Twiss_des, rb):
    chi_2_min = 1e30
    n_iter    = 0
    A0        = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)

    def prt_iter(prm, chi_2, Twiss_des, Twiss_k):

        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        phi = lat_prop.compute_phi_lat()
        phi_d1 = compute_phi_bend(lat_prop, d1_bend._bend_list)
        phi_d2 = compute_phi_bend(lat_prop, d2_bend._bend_list)
        phi_rb = lat_prop.get_phi_elem(rb, 0)

        print("\n{:4d} chi_2 = {:9.3e}".format(n_iter, chi_2))
        print("    [eta_x, eta'_x] = "+
              "[{:10.3e}, {:10.3e}] ([{:10.3e}, {:10.3e}])".
              format(Twiss_k[0][ind.x], Twiss_k[0][ind.px],
                     Twiss_des[0][ind.x], Twiss_des[0][ind.px]))
        print("    alpha           = "
              +"[{:10.3e}, {:10.3e}] ([{:10.3e}, {:10.3e}])".
              format(Twiss_k[1][ind.X], Twiss_k[1][ind.Y],
                     Twiss_des[1][ind.X], Twiss_des[1][ind.Y]))
        print("    beta            = "+
              "[{:8.5f}, {:10.5f}]   ([{:8.5f}, {:8.5f}])".
              format(Twiss_k[2][ind.X], Twiss_k[2][ind.Y],
                     Twiss_des[2][ind.X], Twiss_des[2][ind.Y]))
        print("    dnu             = "+
              "[{:8.5f}, {:10.5f}]   ([{:8.5f}, {:8.5f}])".
              format(Twiss_k[3][ind.X], Twiss_k[3][ind.Y],
                     Twiss_des[3][ind.X], Twiss_des[3][ind.Y]))
        print("\n    phi_sp          = {:8.5f}".format(phi))
        print("    C [m]           = {:8.5f}".format(lat_prop.compute_circ()))
        print("\n    phi_d2          = {:8.5f}".format(phi_d2))
        print("    phi_d1          = {:8.5f}".format(phi_d1))
        print("    phi_rb          = {:8.5f}".format(phi_rb))
        prm_list.prt_prm(prm)

    def compute_chi_2_Twiss(Twiss_des, Twiss_k):
        prt = not False

        dchi_2 = weight[0]*(Twiss_k[0][ind.x]-Twiss_des[0][ind.x])**2
        if prt:
            print("\n  dchi2(eta_x)    = {:9.3e}".format(dchi_2))
        chi_2 = dchi_2

        dchi_2 = weight[0]*(Twiss_k[0][ind.px]-Twiss_des[0][ind.px])**2
        if prt:
            print("  dchi2(eta'_x)   = {:9.3e}".format(dchi_2))
        chi_2 += dchi_2

        dchi_2 = weight[1]*(Twiss_k[1][ind.X]-Twiss_des[1][ind.X])**2
        if prt:
            print("  dchi2(alpha_x)  = {:9.3e}".format(dchi_2))
        chi_2 += dchi_2

        dchi_2 = weight[1]*(Twiss_k[1][ind.Y]-Twiss_des[1][ind.Y])**2
        if prt:
            print("  dchi2(alpha_y)  = {:9.3e}".format(dchi_2))
        chi_2 += dchi_2

        dchi_2 = weight[2]*(Twiss_k[2][ind.X]-Twiss_des[2][ind.X])**2
        if prt:
            print("  dchi2(beta_x)   = {:9.3e}".format(dchi_2))
        chi_2 += dchi_2

        dchi_2 = weight[3]*(Twiss_k[2][ind.Y]-Twiss_des[2][ind.Y])**2
        if prt:
            print("  dchi2(beta_y)   = {:9.3e}".format(dchi_2))
        chi_2 += dchi_2

        dchi_2 = weight[4]*(Twiss_k[3][ind.X]-Twiss_des[3][ind.X])**2
        if prt:
            print("  dchi2(dnu_x)    = {:9.3e}".format(dchi_2))
        chi_2 += dchi_2

        dchi_2 = weight[5]*(Twiss_k[3][ind.Y]-Twiss_des[3][ind.Y])**2
        if prt:
            print("  dchi2(dnu_y)    = {:9.3e}".format(dchi_2))
        chi_2 += dchi_2

        return chi_2

    def compute_chi_2(A0, Twiss_des):

        A1 = _copy.copy(A0)
        lat_prop._lattice.propagate(
            lat_prop._model_state, A1, uc_list[1]+1, sp_list[1]-uc_list[1])
        Twiss_k = cs.compute_Twiss_A(A1.jacobian())

        A1 = _copy.copy(A0)
        lat_prop._lattice.propagate(
            lat_prop._model_state, A1, uc_list[1]+1, sp_list[0]-uc_list[1])
        _, _, _, dnu = cs.compute_Twiss_A(A1.jacobian())
        Twiss_k[3][:] = 2e0*(Twiss_k[3][:]-dnu[:])

        chi_2 = compute_chi_2_Twiss(Twiss_des, Twiss_k)

        return chi_2, Twiss_k

    def f_match(prm):
        nonlocal chi_2_min, n_iter, A0, Twiss_des

        n_iter += 1
        prm_list.set_prm(prm)
        if d1_bend != []:
            d1_bend.correct_bend_phi()

        chi_2, Twiss_k = compute_chi_2(A0, Twiss_des)
        if chi_2 < chi_2_min:
            prt_iter(prm, chi_2, Twiss_des, Twiss_k)
            pc.prt_lat(lat_prop, "match_lat_k.txt", prm_list, d1_bend=d1_bend)
            chi_2_min = min(chi_2, chi_2_min)

            # Problematic => system not time invariant.
            _, A = compute_periodic_cell(lat_prop, uc_list)
            A0.set_jacobian(A_7x7)
        else:
            print("\n{:3d} {:9.3e} ({:9.3e})".
                  format(n_iter, chi_2, chi_2-chi_2_min))
            prm_list.prt_prm(prm)

        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    g_tol    = 1e-5
    x_tol    = 1e-4

    Twiss_0, A = compute_periodic_cell(lat_prop, uc_list)

    print("\nmatch_straight:\n\nEntrance:")
    lat_prop.prt_Twiss_param(Twiss_0)

    A0.set_zero()
    A_7x7 = np.zeros((7, 7))
    A_7x7[:6, :6] = A
    A0.set_jacobian(A_7x7)

    prm, bounds = prm_list.get_prm()

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    if not True:
        minimum = opt.minimize(
            f_match,
            prm,
            method="Powell",
            # callback=prt_iter,
            bounds = bounds,
            options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
        )
    else:
        minimum = opt.minimize(
            f_match,
            prm,
            method="CG",
            # callback=prt_iter,
            jac=None,
            bounds = bounds,
            options={"gtol": g_tol, "maxiter": max_iter}
        )

    print("\n".join(minimum))


# TPSA max order.
gtpsa_prop.no = 2

cod_eps = 1e-15
E_0     = 3.0e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

lat_prop.prt_lat("max_4u_match_lat.txt")

uc_list = np.zeros(2, dtype=int)
uc_list[0] = lat_prop._lattice.find("r1_f1", 0).index
uc_list[1] = lat_prop._lattice.find("r1_f1", 1).index

sp_list = np.zeros(2, dtype=int)
sp_list[0] = lat_prop._lattice.find("lsborder", 1).index
sp_list[1] = lat_prop._lattice.find("cav", 0).index

print("\nunit cell entrance           {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_list[0]].name, uc_list[0]))
print("unit cell exit               {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_list[1]].name, uc_list[1]))
print("super period last sextupole  {:5s} loc = {:d}".
      format(lat_prop._lattice[sp_list[0]].name, sp_list[0]))
print("super period exit            {:5s} loc = {:d}".
      format(lat_prop._lattice[sp_list[1]].name, sp_list[1]))

Twiss_des = np.array([eta_des, alpha_des, beta_des, dnu_des])

# Weights.
weight = np.array([
    1e2,  # [eta_x, eta'_x] at the centre of the straight.
    1e0,  # alpha at the centre of the straight.
    1e-2, # beta_x at the centre of the straight.
    1e-3, # beta_y at the centre of the straight.
    0e-4, # dnu_x across the straight.
    0e-4  # dnu_y across the straight.
])

d2_list = ["d2_f1_sl_d0a", "d2_f1_sl_d0b", "d2_f1_sl_d0c", "d2_f1_sl_df1",
           "d2_f1_sl_df2", "d2_f1_sl_df3", "d2_f1_sl_df4", "d2_f1_sl_df5"]

d1_list = [
    "d1_f1_sl_ds6", "d1_f1_sl_ds5", "d1_f1_sl_ds4", "d1_f1_sl_ds3",
    "d1_f1_sl_ds2", "d1_f1_sl_ds1", "d1_f1_sl_ds0",
    "d1_f1_sl_dm1", "d1_f1_sl_dm2", "d1_f1_sl_dm3", "d1_f1_sl_dm4",
    "d1_f1_sl_dm5"
]

d1_bend = pc.bend_class(lat_prop, d1_list, phi_max, b_2_max)
d2_bend = pc.bend_class(lat_prop, d2_list, phi_max, b_2_max)

# Remark:
# Using the bend angles as parameters - maintaining the total - will change the
# initial conditions for the horizontal dipspersion; i.e., requires an iterative
# approach.

prms = [
    ("q1_f1", "b_2"),
    ("q2_f1", "b_2"),
    ("q3_f1", "b_2"),

    ("q2_f1", "phi"),
    ("b_2_bend", d1_bend),

    # ("d1_f1_sl_ds6", "phi"),
    # ("d1_f1_sl_ds5", "phi"),
    # ("d1_f1_sl_ds4", "phi"),
    # ("d1_f1_sl_ds3", "phi"),
    # ("d1_f1_sl_ds2", "phi"),
    # ("d1_f1_sl_ds1", "phi"),
    # ("d1_f1_sl_ds0", "phi"),
    # ("d1_f1_sl_dm1", "phi"),
    # ("d1_f1_sl_dm2", "phi"),
    # ("d1_f1_sl_dm3", "phi"),
    # ("d1_f1_sl_dm4", "phi"),
    # ("d1_f1_sl_dm5", "phi")
]

dprm_list = np.array([
    1e-3, 1e-3, 1e-3,
    1e-3, 1e-3, 1e-3, 1e-3,
    1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3,
    1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3
])

prm_list = pc.prm_class(lat_prop, prms, b_2_max)

rb = "r1_f1"

match_straight(
    lat_prop, prm_list, uc_list, sp_list, weight, d1_bend, d2_bend, Twiss_des,
    rb)
