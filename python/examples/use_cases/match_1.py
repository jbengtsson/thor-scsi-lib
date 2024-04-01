"""Use Case:
     Matching the centre of a straight section with a triplet.
   Input:
     Lattice with one super period.
   Constraints:
    [alpha_x,y, beta_x,y, eta_x, eta'_x] at the centre of the straight,
    minimize linear chromaticity.
   Parameters:
     Triplet gradients (3) & spacing (3).
"""
import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

from dataclasses import dataclass
import os
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import gtpsa

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    get_set_mpole as gs, index_class as ind, courant_snyder as cs

from thor_scsi.utils.output import vec2txt


# check that a is a vector, as an example
# _, = a.shape
# compare to
# np.atleast_2d


@dataclass
class MatchingState:
    n_iter: float = 0
    #: :math:`\chiˆ2`
    chi_2: float = np.nan


def match_straight(
        lat_prop, prm_list, bounds, Twiss_0, Twiss_1, Twiss_1_design, weight):

    chi_2_min = 1e30
    n_iter    = 0
    A0        = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)

    def get_prm(lat_prop, prm_list):
        # Dictionary of parameter types and corresponding get functions.
        how_to_get_prm = \
            {"b_2": lat_prop.get_b_2_elem, "L": lat_prop.get_L_elem}
        prms = []
        for k in range(len(prm_list)):
            prms.append(how_to_get_prm[prm_list[k][1]](prm_list[k][0], 0))
        return prms

    def set_prm(lat_prop, prm_list, prms):
        # Dictionary of parameter types and corresponding set functions.
        how_to_set_prm = {"b_2": lat_prop.set_b_2_fam, "L": lat_prop.set_L_fam}
        for k in range(len(prm_list)):
            how_to_set_prm[prm_list[k][1]](prm_list[k][0], prms[k])

    def prt_iter(prms, chi_2, Twiss_k):
        lat_prop.prt_Twiss_param(Twiss_k)
        print(
            f"\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n\n  prms   =" + vec2txt(prms))

    def compute_chi_2(prms):
        nonlocal A0

        A1 = A0
        lat_prop._lattice.propagate(
            lat_prop._model_state, A1, 0, len(lat_prop._lattice))
        chi_2 = 0e0
        Twiss_k = cs.compute_Twiss_A(A1.jacobian())
        for k in range(3):
            chi_2 += weight[k] * np.sum((Twiss_k[k] - Twiss_1_design[k]) ** 2)
        return chi_2, Twiss_k

    def f_match(prms):
        nonlocal chi_2_min
        nonlocal n_iter

        n_iter += 1
        set_prm(lat_prop, prm_list, prms)
        chi_2, Twiss_k = compute_chi_2(prms)
        if True or chi_2 < chi_2_min:
            prt_iter(prms, chi_2, Twiss_k)
            chi_2_min = min(chi_2, chi_2_min)
        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    print("\nmatch_straight:\n\nEntrance:")
    lat_prop.prt_Twiss_param(Twiss_0)
    print("\nExit:")
    lat_prop.prt_Twiss_param(Twiss_1)
    print("\nDesired:")
    lat_prop.prt_Twiss_param(Twiss_1_design)

    A0.set_zero()
    # Use *-operator to unpack the list of arguments.
    A_7x7 = np.zeros((7, 7))
    A_7x7[:6, :6] = cs.compute_A(*Twiss_0[:3])
    A0.set_jacobian(A_7x7)

    prms = get_prm(lat_prop, prm_list)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol
    # CG     gtol
    minimum = optimize.minimize(
        f_match,
        prms,
        method="CG",
        bounds=bounds,
        options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter},
    )


def chk_lat(lat_name, lat_prop, Twiss_0):
    n_dof = 2
    A0 = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    A0.set_zero()
    # Use *-operator to unpack the list of arguments.
    A_7x7 = np.zeros((7, 7))
    A_7x7[:6, :6] = cs.compute_A(*Twiss_0[:3])
    A0.set_jacobian(A_7x7)

    A1 = A0
    lat_prop._lattice.propagate(
        lat_prop._model_state, A1, 0, len(lat_prop._lattice))
    lat_prop._data = \
        lo.compute_Twiss_along_lattice(
            n_dof, lat_prop._lattice, lat_prop._model_state, A=A1,
            desc=lat_prop._desc, mapping=lat_prop._named_index)

    Twiss_k = cs.compute_Twiss_A(A1.jacobian())
    print("\nchk_lat - Exit:")
    lat_prop.prt_Twiss()
    lat_prop.plt_Twiss(lat_name+".png", not False)


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 1
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

ind = ind.index_class()

if not True:
    home_dir = os.path.join(
        os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "BESSY-III",
        "ipac_2024")
    lat_name = "2024Mar"
else:
    home_dir = \
        os.path.join(
            os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_4U")
    lat_name = "max_4u_match"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

types = lat_prop.get_types()

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".
      format(lat_prop.compute_phi()))

# Entrance.
eta   = np.array([0.02273, 0.0])
alpha = np.array([0.0, 0.0])
beta  = np.array([0.92011, 9.72555])
nu    = np.array([np.nan, np.nan])
Twiss_0 = eta, alpha, beta, nu

if False:
    chk_lat(lat_name, lat_prop, Twiss_0)

# Exit.
eta   = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta  = np.array([9.2, 2.0])
nu    = np.array([0e0, 0e0])
Twiss_1 = eta, alpha, beta, nu

# Desired.
eta   = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta  = np.array([9.2, 2.0])
nu    = np.array([0e0, 0e0])
Twiss_1_design = eta, alpha, beta, nu

# Weights: eta, alpha, beta.
weight = np.array([
    1e0, # eta.
    1e3, # alpha.
    1e0  # beta.
])

# Triplet parameter family names & type.
prm_list = [
    ("qd",  "b_2"),
    ("qf2", "b_2"),
    # ("ds5", "b_2"),
    # ("ds4", "b_2"),
    # ("ds3", "b_2"),
    # ("ds2", "b_2"),
    # ("ds1", "b_2"),
    # ("ds0", "b_2"),
    # ("dm1", "b_2"),
    # ("dm2", "b_2"),
    # ("dm3", "b_2"),
    # ("dm4", "b_2"),
    # ("dm5", "b_2"),
    ("d5",  "L")
]

# Max parameter range.
b_2_max = 10.0
L_min   = 0.1

bounds = [
    (-b_2_max, 0.0),     # qd b_2.
    ( 0.0,     b_2_max), # qf2 b_2.
    # (-b_2_max, 0.0),   # ds5 b_2.
    # (-b_2_max, 0.0),   # ds4 b_2.
    # (-b_2_max, 0.0),   # ds3 b_2.
    # (-b_2_max, 0.0),   # ds2 b_2.
    # (-b_2_max, 0.0),   # ds1 b_2.
    # (-b_2_max, 0.0),   # ds0 b_2.
    # (-b_2_max, 0.0),   # dm1 b_2.
    # (-b_2_max, 0.0),   # dm2 b_2.
    # (-b_2_max, 0.0),   # dm3 b_2.
    # (-b_2_max, 0.0),   # dm4 b_2.
    # (-b_2_max, 0.0),   # dm5 b_2.
    ( L_min,   0.25)     # d5 L.
]

match_straight(
    lat_prop, prm_list, bounds, Twiss_0, Twiss_1, Twiss_1_design, weight)

if not False:
    plt.show()
