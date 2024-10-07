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

L_min        = 0.10
L_max        = 0.50
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0

eta_des   = [0.0, 0.0]   # [eta_x, eta'_x] at the centre of the straight.
alpha_des = [0.0, 0.0]   # Alpha_x,y at the centre of the straight.
beta_des  = [9.2, 2.0]   # Beta_x,y at the centre of the straight.
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


def compute_periodic_cell(lat_prop, uc_0, uc_1):
    M = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    M.set_identity()
    # The 3rd argument is the 1st element index & the 4th the number of elements
    # to propagate across.
    lat_prop._lattice.propagate(
        lat_prop._model_state, M, uc_0, uc_1-uc_0+1)
    stable, nu, A, A_inv, _ = lo.compute_M_diag(2, M.jacobian())
    eta, Twiss = lo.transform_matrix_extract_twiss(A)
    alpha = np.array((Twiss[ind.X][0], Twiss[ind.Y][0]))
    beta = np.array((Twiss[ind.X][1], Twiss[ind.Y][1]))
    Twiss = eta, alpha, beta
    print("\ncompute_periodic_cell:")
    lat_prop.prt_Twiss_param(Twiss)
    return Twiss, A


def match_straight(
        lat_prop, prm_list, uc_0, uc_1, sp_1, sp_2, weight, phi_lat, Twiss_des):
    chi_2_min = 1e30
    n_iter    = 0
    A0        = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)

    def prt_iter(prm, chi_2, Twiss_des, Twiss_k):

        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        rb = "r1"

        phi = lat_prop.compute_phi_lat()
        phi_d2 = compute_phi_bend(lat_prop, d2_list)
        phi_d1 = compute_phi_bend(lat_prop, d1_list)
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
            lat_prop._model_state, A1, uc_1+1, sp_2-uc_1)
        Twiss_k = cs.compute_Twiss_A(A1.jacobian())
        print("\ncompute_chi_2:")

        A1 = _copy.copy(A0)
        lat_prop._lattice.propagate(
            lat_prop._model_state, A1, uc_1+1, sp_1-uc_1)
        _, _, _, dnu = cs.compute_Twiss_A(A1.jacobian())
        Twiss_k[3][:] = 2e0*(Twiss_k[3][:]-dnu[:])

        chi_2 = compute_chi_2_Twiss(Twiss_des, Twiss_k)

        return chi_2, Twiss_k

    def f_match(prm):
        nonlocal chi_2_min, n_iter, A0, Twiss_des

        n_iter += 1
        prm_list.set_prm(prm)
        phi_lat.set_phi_lat()

        chi_2, Twiss_k = compute_chi_2(A0, Twiss_des)
        if chi_2 < chi_2_min:
            prt_iter(prm, chi_2, Twiss_des, Twiss_k)
            pc.prt_lat(lat_prop, "match_lat_k.txt", prm_list, phi_lat=phi_lat)
            chi_2_min = min(chi_2, chi_2_min)

            # Problematic => system not time invariant.
            _, A = compute_periodic_cell(lat_prop, uc_0, uc_1)
            A0.set_jacobian(A_7x7)
        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4
    g_tol    = 1e-5

    Twiss_0, A = compute_periodic_cell(lat_prop, uc_0, uc_1)

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

    # Powell ftol, xtol
    # CG     gtol
    minimum = opt.minimize(
        f_match,
        prm,
        method="CG",
        # callback=prt_iter,
        # bounds = bounds,
        # options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
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

if True:
    uc_0 = lat_prop._lattice.find("d2_0", 7).index
    uc_1 = lat_prop._lattice.find("d2_0", 8).index
else:
    uc_0 = lat_prop._lattice.find("q4_h", 7).index
    uc_1 = lat_prop._lattice.find("q4_h", 8).index
sp_1 = lat_prop._lattice.find("s1", 1).index
sp_2 = lat_prop._lattice.find("cav", 0).index
print("\nunit cell entrance           {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_0].name, uc_0))
print("unit cell exit               {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_1].name, uc_1))
print("super period last sextupole  {:5s} loc = {:d}".
      format(lat_prop._lattice[sp_1].name, sp_1))
print("super period exit            {:5s} loc = {:d}".
      format(lat_prop._lattice[sp_2].name, sp_2))

Twiss_des = np.array([eta_des, alpha_des, beta_des, dnu_des])

# Weights.
weight = np.array([
    1e1,  # eta_x at the centre of the straight.
    1e0,  # alpha at the centre of the straight.
    0e0,  # beta_x at the centre of the straight.
    0e0,  # beta_y at the centre of the straight.
    1e-4, # dnu_x across the straight.
    1e-4  # dnu_y across the straight.
])

d2_list = ["d2_0", "d2_1", "d2_2", "d2_3", "d2_4", "d2_5"]

d1_list = [
    "d1_u6", "d1_u5", "d1_u4", "d1_u3", "d1_u2", "d1_u1", "d1_0",
    "d1_d1", "d1_d2", "d1_d3", "d1_d4", "d1_d5"
]

d2_bend = pc.bend_class(lat_prop, d2_list, phi_max, b_2_max)
d1_bend = pc.bend_class(lat_prop, d1_list, phi_max, b_2_max)

# Remark:
# Using the bend angles as parameters - maintaining the total - will change the
# initial conditions for the horizontal dipspersion; i.e., requires an iterative
# approach.

prm_list = [
    ("q1",       "b_2"),
    ("q2",       "b_2"),
    ("q3",       "b_2"),
    ("b_2_bend", d1_bend),

    ("phi_bend", d2_bend)
]

prm_list = pc.prm_class(lat_prop, prm_list, b_2_max)

# To maintain the total bend angle.
phi_lat = pc.phi_lat_class(lat_prop, 2, d1_bend)

match_straight(
    lat_prop, prm_list, uc_0, uc_1, sp_1, sp_2, weight, phi_lat, Twiss_des)
