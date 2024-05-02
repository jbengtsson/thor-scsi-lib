"""Use Case:
     Parametric scans/evaluations for a unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

import copy as _copy
from dataclasses import dataclass
import os
from typing import Tuple

import numpy as np
from scipy import optimize as opt

import gtpsa

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, get_set_mpole as gs, \
    index_class as ind, linear_optics as lo, prm_class as pc, \
    courant_snyder as cs
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()

L_min        = 0.10
L_max        = 0.50
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0


def match_straight(lat_prop, prm_list, Twiss_0, Twiss_1, weight):

    chi_2_min = 1e30
    n_iter    = 0
    A0        = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)

    def prt_iter(prm, chi_2, Twiss_k):
        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        b_list = [
            "b2u_6", "b2u_5", "b2u_4", "b2u_3", "b2u_2", "b2u_1", "b2_0",
            "b2d_1", "b2d_2", "b2d_3", "b2d_4", "b2d_5"]
        rb = "qf1"

        phi = lat_prop.compute_phi_lat()
        phi_b = compute_phi_bend(lat_prop, b_list)
        phi_rb = lat_prop.get_phi_elem(rb, 0)

        print("\n{:4d} chi_2 = {:9.3e}".format(n_iter, chi_2))
        lat_prop.prt_Twiss_param(Twiss_k[:3])
        print("\n  phi_sp =  {:8.5f}".format(phi))
        print("  C [m]  =  {:8.5f}".format(lat_prop.compute_circ()))
        print("\n  phi_b  =  {:8.5f}".format(phi_b))
        print("  phi_rb =  {:8.5f}".format(phi_rb))
        prm_list.prt_prm(prm)

    def compute_chi_2_Twiss(Twiss_k):
        chi_2 = 0e0
        for j in range(3):
            for k in range(2):
                chi_2 += weight[j] * ((Twiss_1[j][k] - Twiss_k[j][k]) ** 2)
        return chi_2

    def compute_chi_2():
        nonlocal A0

        A1 = _copy.copy(A0)
        lat_prop._lattice.propagate(
            lat_prop._model_state, A1, 0, len(lat_prop._lattice))
        Twiss_k = cs.compute_Twiss_A(A1.jacobian())
        chi_2 = compute_chi_2_Twiss(Twiss_k)
        return chi_2, Twiss_k

    def f_match(prm):
        nonlocal chi_2_min, n_iter

        n_iter += 1
        prm_list.set_prm(lat_prop, prm)

        chi_2, Twiss_1 = compute_chi_2()
        if chi_2 < chi_2_min:
            prt_iter(prm, chi_2, Twiss_1)
            pc.prt_lat(lat_prop, "match_lat_k.txt", prm_list)
            chi_2_min = min(chi_2, chi_2_min)
        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    print("\nmatch_straight:\n\nEntrance:")
    lat_prop.prt_Twiss_param(Twiss_0)
    print("\nDesired:")
    lat_prop.prt_Twiss_param(Twiss_1)

    A0.set_zero()
    # Use *-operator to unpack the list of arguments.
    A_7x7 = np.zeros((7, 7))
    A_7x7[:6, :6] = cs.compute_A(*Twiss_0[:3])
    A0.set_jacobian(A_7x7)

    prm, bounds = prm_list.get_prm(lat_prop)

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
        bounds = bounds,
        options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
    )

    print("\n".join(minimum))


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 2
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_4U")
lat_name = "max_4u_jb_2"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.get_types()

print("\nTotal bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))
print("Circumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))

# Entrance Twiss parameters.
eta   = np.array([0.03074, 0.0])
alpha = np.array([0.0, 0.0])
beta  = np.array([0.92126, 9.72555])
Twiss_0 = eta, alpha, beta

# Desired exit Twiss parameters.
eta   = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta  = np.array([9.2, 2.0])
Twiss_1 = eta, alpha, beta

# Weights.
weight = np.array([
    1e8,  # eta.
    1e4,  # alpha.
    0*1e0 # beta.
])

bend_list = [
    "b1_0", "b1_1", "b1_2", "b1_3", "b1_4", "b1_5",
    "b2u_6", "b2u_5", "b2u_4", "b2u_3", "b2u_2", "b2u_1", "b2d_1", "b2d_2",
    "b2d_3", "b2d_4", "b2d_5"
]

opt_phi = pc.opt_phi_class(lat_prop, "b2_0", bend_list, phi_max)

prm_list = [
    ("qd",  "b_2"),
    ("qf2", "b_2"),

    ("opt_phi",  opt_phi)
]

prm_list = pc.prm_class(lat_prop, prm_list, b_2_max)

match_straight(lat_prop, prm_list, Twiss_0, Twiss_1, weight)


if False:
    dip_list = [bend]
    dip_list.extend(bend_list)
    prt_lat(lat_prop, "opt_sp.txt", dip_list)
