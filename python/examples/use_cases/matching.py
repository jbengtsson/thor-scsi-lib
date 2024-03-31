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

import copy as _copy
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


# Global variables.
n_iter    = 0
chi_2_min = 1e30


def get_prm(lat_prop, prm_list):
    # Dictionary of parameter types and corresponding get functions.
    how_to_get_prm = {"b_2": lat_prop.get_b_2_elem, "L": lat_prop.get_L_elem}

    prms = []
    for k in range(len(prm_list)):
        prms.append(how_to_get_prm[prm_list[k][1]](prm_list[k][0], 0))
    return prms


def set_prm(lat_prop, prm_list, prms):
    # Dictionary of parameter types and corresponding set functions.
    how_to_set_prm = {"b_2": lat_prop.set_b_2_fam, "L": lat_prop.set_L_fam}

    for k in range(len(prm_list)):
        how_to_set_prm[prm_list[k][1]](prm_list[k][0], prms[k])


@dataclass
class MatchingState:
    n_iter: float = 0
    #: :math:`\chiˆ2`
    chi_2: float = np.nan


def match_straight(
        lat_prop, loc, prm_list, bounds, Twiss_0, Twiss_1, Twiss_1_design,
        weight):
    """ """

    def compute_chi_2(
            A0: np.ndarray,
            chi_2_min: float,
            prms, copy: bool = True
    ) -> Tuple[float, float]:
        """Computes weighted sum square
        """
        # Local copy.
        if copy:
            A0 = _copy.copy(A0)

        lat_prop._lattice.propagate(
            lat_prop._model_state, A0, loc, len(lat_prop._lattice) - loc)

        chi_2 = 0e0

        # returns Twiss parameters [alpha, beta, eta] and dnu
        Twiss_k = cs.compute_Twiss_A(A0.jacobian())
        # todo: sum to twiis
        for k in range(3):
            chi_2 += weight[k] * np.sum((Twiss_k[k] - Twiss_1_design[k]) ** 2)

        # py ?
        if weight[3] != 0e0:
            M = \
                lo.compute_map(
                    lat_prop._lattice, lat_prop._model_state,
                    desc=lat_prop._desc, tpsa_order=2)
            stable, _, xi = lo.compute_nu_xi(lat_prop._desc, 2, M)
            if stable:
                dchi_2 = weight[3] * np.sum(xi ** 2)
            else:
                dchi_2 = 1e30
            chi_2 += dchi_2

        if chi_2 < chi_2_min:
            chi_2_min = chi_2
            print(f"\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n\n  prms   ="
                  + vec2txt(prms))
            lat_prop.prt_Twiss_param(Twiss_k)
            if weight[3] != 0e0:
                print(f"\n  xi     = [{xi[ind.X]:5.3f}, {xi[ind.Y]:5.3f}]")
                print(f"  dchi_2 = {dchi_2:9.3e}")
        return chi_2, chi_2_min

    def f_match(prms):
        global n_iter
        global chi_2_min

        n_iter += 1
        set_prm(lat_prop, prm_list, prms)
        chi_2, chi_2_min = compute_chi_2(A0, chi_2_min, prms)

        return chi_2

    def prt_result(f_match, prms0, minimum):
        global n_iter
        global chi_2_min

        # Compute new Twiss parameters along lattice.
        stable = lat_prop.comp_per_sol()
        lat_prop.plt_Twiss("after.png", True)

        print("\nInitial parameter values:")
        n_iter = 0
        chi_2_min = 1e30
        f_match(prms0)
        print("\nFinal parameter values:")
        n_iter = 0
        chi_2_min = 1e30
        f_match(minimum["x"])
        print("\n Minimum:\n", minimum)

    global n_iter
    global chi_2_min

    max_iter = 1000
    f_tol = 1e-4
    x_tol = 1e-4

    print("\nmatch_straight:\n\nloc = ", loc)
    print("\nTwiss parameters at centre of super period:")
    lat_prop.prt_Twiss_param(Twiss_0)
    print("\nTwiss parameters at centre of straight:")
    lat_prop.prt_Twiss_param(Twiss_1)
    print("\nDesired Twiss parameters at the centre of straight:")
    lat_prop.prt_Twiss_param(Twiss_1_design)

    A0 = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    A0.set_zero()
    # Use *-operator to unpack the arguments from a list.
    A_7x7 = np.zeros((7, 7))
    A_7x7[:6, :6] = cs.compute_A(*Twiss_0[:3])
    A0.set_jacobian(A_7x7)

    # Initialise parameters.
    prms1 = prms0 = get_prm(lat_prop, prm_list)
    print("\nprms = ", prms1)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.
    minimum = optimize.minimize(
        f_match,
        prms1,
        method="Powell",
        bounds=bounds,
        options={"xtol": x_tol, "ftol": f_tol, "maxiter": max_iter},
    )

    prt_result(f_match, prms0, minimum)


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
    lat_name = "max_4u_sup_per"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
if not stable:
    assert False

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".
      format(lat_prop.compute_phi()))
lat_prop.prt_M()
lat_prop.prt_lat_param()

types = lat_prop.get_types()
lat_prop.prt_Twiss()
lat_prop.plt_Twiss("before.png", False)

# Get dispersion & Twiss parameters at centre of the last unit cell.
loc = lat_prop._lattice.find("sf_h", 10).index
print(f"\nCentre of super period: {loc:1d} {lat_prop._lattice[loc].name:s}")
Twiss_0 = lat_prop.get_Twiss(loc)
lat_prop.prt_Twiss_param(Twiss_0)

# Get dispersion & Twiss parameters at centre of straight.
# Last element is the RF cavity.
loc = len(lat_prop._lattice) - 2
Twiss_1 = lat_prop.get_Twiss(loc)
print(f"\nCentre of straight: {loc:1d} {lat_prop._lattice[loc].name:s}")
lat_prop.prt_Twiss_param(Twiss_1)

# Desired dispersion & Twiss parameters at centre of straight.
eta   = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta  = np.array([2.5, 2.5])
nu    = np.array([0e0, 0e0])
Twiss_1_design = eta, alpha, beta, nu

# Weights: eta, alpha, beta, xi.
weight = np.array([
    1e0,  # eta.
    1e3,  # alpha.
    1e-5, # beta.
    1e-10 # xi.
])

# Triplet parameter family names & type.
prm_list = [
    ("qd",  "b_2"),
    ("qf2", "b_2"),
    ("d5", "L"),
    ("d6", "L"),
]

# Max parameter range.
b_2_max = 10.0
L_min   = 0.1

bounds = [
    (-b_2_max, 0.0),     # qd b_2.
    ( 0.0,     b_2_max), # qf2 b_2.
    ( L_min,   0.25),    # d5 L.
    ( L_min,   0.25)     # d6 L.
]

# Zero sextopoles.
# Compute linear chromaticity.

lat_prop.set_b_n_fam("sf_h", 3, 0e0)
lat_prop.set_b_n_fam("sd1",  3, 0e0)
lat_prop.set_b_n_fam("sd2",  3, 0e0)

M = \
    lo.compute_map(
        lat_prop._lattice, lat_prop._model_state, desc=lat_prop._desc,
        tpsa_order=no)
stable, nu, xi = lo.compute_nu_xi(lat_prop._desc, no, M)
print("\n  nu = [{:7.5f}, {:7.5f}]".format(nu[ind.X], nu[ind.Y]))
print("  xi = [{:7.5}, {:7.5}]".format(xi[ind.X], xi[ind.Y]))

match_straight(
    lat_prop, loc, prm_list, bounds, Twiss_0, Twiss_1, Twiss_1_design, weight)

if not False:
    plt.show()
print("\nPlots saved as: before.png & after.png")
