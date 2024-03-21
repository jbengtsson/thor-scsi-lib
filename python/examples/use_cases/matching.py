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
logging.basicConfig(level="INFO")
logger = logging.getLogger("thor_scsi")

import copy as _copy
from dataclasses import dataclass
import os
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    get_set_mpole as gs, index_class as ind


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole = 3


# check that a is a vector, as an example
# _, = a.shape
# compare to
# np.atleast_2d


# Global variables.
n_iter = 0
chi_2_min = 1e30


def get_prm(lat, prm_list):
    # Dictionary of parameter types and corresponding get functions.
    how_to_get_prm = {"b_2": get_b_2_elem, "L": get_L_elem}

    prms = []
    for k in range(len(prm_list)):
        prms.append(how_to_get_prm[prm_list[k][1]](lat, prm_list[k][0], 0))
    return prms


def set_prm(lat, prm_list, prms):
    # Dictionary of parameter types and corresponding set functions.
    how_to_set_prm = {"b_2": set_b_2_fam, "L": set_L_fam}

    for k in range(len(prm_list)):
        how_to_set_prm[prm_list[k][1]](lat, prm_list[k][0], prms[k])


@dataclass
class MatchingState:
    n_iter: float = 0
    #: :math:`\chiˆ2`
    chi_2: float = np.nan


def match_straight(
        lat_prop, loc, prm_list, bounds, Twiss0, Twiss1, Twiss1_design, weight):
    """ """

    def compute_chi_2(
            A0: np.ndarray, chi_2_min: float, prms, copy: bool = True
    ) -> Tuple[float, float]:
        """Computes weighted sum square
        """
        # Local copy.
        if copy:
            A0 = _copy.copy(A0)

        lat.propagate(model_state, A0, loc, len(lat) - loc)

        chi_2 = 0e0

        # returns Twiss parameters [alpha, beta, eta] and dnu
        Twiss_k = cs.compute_Twiss_A(A0.jacobian())[:3]
        # todo: sum to twiis
        for k in range(3):
            chi_2 += weight[k] * np.sum((Twiss_k[k] - Twiss1_design[k]) ** 2)

        # py ?
        if weight[3] != 0e0:
            M = \
                lo.compute_map(lat, model_state, desc=desc,
                               tpsa_order=tpsa_order)
            stable, _, xi = lo.compute_nu_xi(desc, tpsa_order, M)
            if stable:
                dchi_2 = weight[3] * np.sum(xi ** 2)
            else:
                dchi_2 = 1e30
            chi_2 += dchi_2

        if chi_2 < chi_2_min:
            chi_2_min = chi_2
            print(f"\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n\n  prms   ="
                  + vec2txt(prms))
            prt_Twiss("\n", Twiss_k)
            if weight[3] != 0e0:
                print(f"\n  xi     = [{xi[X_]:5.3f}, {xi[Y_]:5.3f}]")
                print(f"  dchi_2 = {dchi_2:9.3e}")
        return chi_2, chi_2_min

    def f_match(prms):
        global n_iter
        global chi_2_min

        n_iter += 1
        set_prm(lat, prm_list, prms)
        chi_2, chi_2_min = compute_chi_2(A0, chi_2_min, prms)

        return chi_2

    def prt_result(f_match, prms0, minimum):
        global n_iter
        global chi_2_min

        # Compute new Twiss parameters along lattice.
        M, A, data = compute_periodic_solution(lat, model_state)
        plt_Twiss(data, "after.png", "Linear Optics - After")

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
    prt_Twiss("\nTwiss parameters at centre of super period:\n", Twiss0)
    prt_Twiss("\nTwiss parameters at centre of straight:\n", Twiss1)
    prt_Twiss("\nDesired Twiss parameters at the centre of straight:\n",
              Twiss1_design)

    A0 = gtpsa.ss_vect_tpsa(desc, 1)
    A0.set_zero()
    A0.set_jacobian(cs.compute_A(*Twiss0))

    # Initialise parameters.
    prms1 = prms0 = get_prm(lat, prm_list)
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

t_dir = \
    os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB",
                 "BESSY-III", "ipac_2024")
file_name = os.path.join(t_dir, "2024Mar.lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)
get_set = gs.get_set_mpole_class()

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".
      format(get_set.compute_phi(lat_prop._lattice)))
lat_prop.prt_M()
if stable:
    lat_prop.prt_Twiss_param()
else:
    assert False

types = lat_prop.get_types()
lat_prop.prt_Twiss(types)

lat_prop.plt_Twiss("before.png", False)

# Get dispersion & Twiss parameters at centre of the last unit cell.
loc = lat_prop._lattice.find("sf1", 1).index
print(f"\nCentre of super period: {loc:1d} {lat_prop._lattice[loc].name:s}")
Twiss_0 = lat_prop.get_Twiss(loc)

# Compute linear chromaticity.
M = \
    lo.compute_map(lat_prop._lattice, lat_prop._model_state,
                   desc=lat_prop._desc, tpsa_order=no)
stable, nu, xi = lo.compute_nu_xi(lat_prop._desc, no, M)
print("\n  nu = [{:7.5f}, {:7.5f}]".format(nu[ind.X], nu[ind.Y]))
print("  xi = [{:7.5f}, {:7.5f}]".format(xi[ind.X], xi[ind.Y]))

# Twiss parameters at centre of straight.
Twiss_1 = lat_prop.get_Twiss(len(lat_prop._lattice)-1)

# Desired Twiss parameters at centre of straight.
eta   = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta  = np.array([2.5, 2.5])
Twiss_1_design = eta, alpha, beta

# Weights: eta, alpha, beta, xi.
weight = np.array([
    1e0,
    1e3,
    1e0,                     # xi.
    1e-5
])

# Triplet parameter family names & type.
prm_list = [
    ("uq1", "b_2"),
    ("uq2", "b_2"),
    ("uq3", "b_2"),
    ("ul1", "L"),
    ("ul2", "L"),
    ("ul3", "L"),
]

# Max parameter range.
b_2_max = 10.0
L_min   = 0.1

bounds = [
    (-b_2_max, 0.0),          # uq1 b_2.
    ( 0.0,     b_2_max),      # uq2 b_2.
    (-b_2_max, 0.0),          # uq3 b_2.
    ( L_min,   0.25),         # ul1 L.
    ( L_min,   0.3),          # ul2 L.
    ( L_min,   0.25)          # ul3 L.
]

# Zero sextopoles.
# Compute linear chromaticity.
get_set.set_b_n_fam(lat_prop._lattice, "sf_h", MultipoleIndex.sextupole, 0e0)
get_set.set_b_n_fam(lat_prop._lattice, "sd_h", MultipoleIndex.sextupole, 0e0)
M = \
    lo.compute_map(
        lat_prop._lattice, lat_prop._model_state, desc=lat_prop._desc,
        tpsa_order=no)
stable, nu, xi = lo.compute_nu_xi(lat_prop._desc, no, M)
print("\n  nu = [{:7.5f}, {:7.5f}]".format(nu[ind.X], nu[ind.Y]))
print("  xi = [{:7.5}, {:7.5}]".format(xi[ind.X], xi[ind.Y]))

match_straight(
    lat_prop._lattice, loc, prm_list, bounds, Twiss_0, Twiss_1, Twiss_1_design,
    weight)

if not False:
    plt.show()
print("\nPlots saved as: before.png & after.png")
