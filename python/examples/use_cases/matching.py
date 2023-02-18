"""Use Case:
     Super period straight section matching with a triplet.
   Constraints:
    [alpha_x,y, beta_x,y, eta_x, eta'_x] at the centre of the straigth.
   Parameters:
     triplet gradients (3) & spacing (3).
"""
import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="INFO")
logger = logging.getLogger("thor_scsi")


from scipy import optimize

import copy as _copy
from dataclasses import dataclass
import math
import os
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as tslib

from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs

# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

tpsa_order = 2

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(6, tpsa_order)

# Configuration space coordinates.
X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta,
]


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole = 3


def plt_Twiss(ds, file_name, title):
    # Turn interactive mode off.
    plt.ioff()

    fig, (gr_1, gr_2) = plt.subplots(2)

    gr_1.set_title(title)
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel(r"$\beta_{x,y}$ [m]")
    gr_1.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), "b-", label=r"$\beta_x$")
    gr_1.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), "r-", label=r"$\beta_y$")
    gr_1.legend()

    gr_2.set_xlabel("s [m]")
    gr_2.set_ylabel(r"$\eta_x$ [m]")
    gr_2.plot(ds.s, ds.dispersion.sel(phase_coordinate="x"), label=r"$\eta_x$")
    fig.tight_layout()
    plt.savefig(file_name)


def compute_periodic_solution(lat, model_state):
    """
    Todo:
        model_state: rename to calculation_configuration or calc_config
    """
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    M, A = lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc)
    print("\nM:\n", mat2txt(M.jacobian()))
    res= cs.compute_Twiss_A(A)
    Twiss = res[:3]
    prt_Twiss("\nTwiss:\n", Twiss)
    print("\nend of twiss\n")
    ds = lo.compute_Twiss_along_lattice(n_dof, lat, model_state, A=A, desc=desc)

    return M, A, ds

# check that a is a vector, as an example
# _, = a.shape
# compare to
# np.atleast_2d

def prt_Twiss(str, Twiss):
    """

    todo:
        rename str e.g. header? prefix?

    """
    # eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    # that way I also check that Twiss has exactly three parameters
    eta, alpha, beta = Twiss
    print(str, end="")
    print(f"  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
    print(f"  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
    print(f"  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")


def get_b_n_elem(lat, fam_name, kid_num, n):
    mp = lat.find(fam_name, kid_num)
    return mp.get_multipoles().get_multipole(n).real


def set_b_n_elem(lat, fam_name, kid_num, n, b_n):
    mp = lat.find(fam_name, kid_num)
    mp.get_multipoles().set_multipole(n, b_n)


def set_b_n_fam(lat, fam_name, n, b_n):
    for mp in lat.elements_with_name(fam_name):
        mp.get_multipoles().set_multipole(n, b_n)


def get_L_elem(lat, fam_name, n_kid):
    elem = lat.find(fam_name, n_kid)
    return elem.get_length()


def set_L_fam(lat, fam_name, L):
    for elem in lat.elements_with_name(fam_name):
        elem.set_length(L)


def get_b_2_elem(lat, fam_name, kid_num):
    return get_b_n_elem(lat, fam_name, kid_num, MultipoleIndex.quadrupole)


def set_b_2_fam(lat, fam_name, b_2):
    set_b_n_fam(lat, fam_name, MultipoleIndex.quadrupole, b_2)


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


def compute_chi_2(lat: tslib.Accelerator, model_state: tslib.ConfigType, A0: np.ndarray, chi_2_min: float, prms, copy: bool = True) -> Tuple[float, float]:
    """Computes weighted sum square

    Args:
        A0: start map
    Returns:
         :math:`\chiˆ2`

    Todo:
    check selection of dimensions
        model state - calc_config?
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
        M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
        stable, _, xi = lo.compute_nu_xi(desc, tpsa_order, M)
        if stable:
            dchi_2 = weight[3] * np.sum(xi**2)
        else:
            dchi_2 = 1e30
        chi_2 += dchi_2

    if chi_2 < chi_2_min:
        chi_2_min = chi_2
        logger.info(f"\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n\n  prms   =" + vec2txt(prms))
        prt_Twiss("\n", Twiss_k)
        if weight[3] != 0e0:
            logger.info(f"\n  xi     = [{xi[X_]:5.3f}, {xi[Y_]:5.3f}]")
            logger.info(f"  dchi_2 = {dchi_2:9.3e}")
    return chi_2, chi_2_min


@dataclass
class MatchingState:
    n_iter: float = 0
    #: :math:`\chiˆ2`
    chi_2: float = np.nan


def match_straight(lat, loc, prm_list, bounds, Twiss0, Twiss1, Twiss1_design, weight):
    """ """

    def f_match(prms):
        global n_iter
        global chi_2_min

        n_iter += 1
        set_prm(lat, prm_list, prms)
        chi_2, chi_2_min = compute_chi_2(lat, model_state, A0, chi_2_min, prms)

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
    prt_Twiss("\nDesired Twiss parameters at the centre of straight:\n", Twiss1_design)

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


def get_Twiss(loc):
    eta = np.array(
        [
            data.dispersion.sel(phase_coordinate="x").values[loc],
            data.dispersion.sel(phase_coordinate="px").values[loc],
            data.dispersion.sel(phase_coordinate="y").values[loc],
            data.dispersion.sel(phase_coordinate="py").values[loc],
        ]
    )
    alpha = np.array(
        [
            data.twiss.sel(plane="x", par="alpha").values[loc],
            data.twiss.sel(plane="y", par="alpha").values[loc],
        ]
    )
    beta = np.array(
        [
            data.twiss.sel(plane="x", par="beta").values[loc],
            data.twiss.sel(plane="y", par="beta").values[loc],
        ]
    )
    return eta, alpha, beta


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
# t_file = os.path.join(t_dir, 'b3_tst.lat')
t_file = os.path.join(t_dir, "b3_sf_40Grad_JB.lat")

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

# Compute Twiss parameters along lattice.
M, A, data = compute_periodic_solution(lat, model_state)
plt_Twiss(data, "before.png", "Linear Optics - Before")

# Get dispersion & Twiss parameters at centre of super period.
loc = lat.find("sf_m", 3).index
print(f"\nCentre of super period: {loc:1d} {lat[loc].name:s}")
Twiss0 = get_Twiss(loc)

# Compute linear chromaticity.
M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
stable, nu, xi = lo.compute_nu_xi(desc, tpsa_order, M)
print("\n  nu = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))
print("  xi = [{:7.5}, {:7.5}]".format(xi[X_], xi[Y_]))

# Twiss parameters at centre of straight.
Twiss1 = get_Twiss(len(lat) - 1)

# Desired Twiss parameters at centre of straight.
eta = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta = np.array([2.5, 2.5])
Twiss1_design = eta, alpha, beta

weight = np.array([1e0, 1e3, 1e0, 1e-5])  # Eta.  # Alpha.  # Beta.  # xi.

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
L_min = 0.1

bounds = [
    (-b_2_max, 0.0),
    (0.0, b_2_max),
    (-b_2_max, 0.0),
    (L_min, 0.25),
    (L_min, 0.3),
    (L_min, 0.25),
]

# Zero sextopoles.
# Compute linear chromaticity.
set_b_n_fam(lat, "sf_h", MultipoleIndex.sextupole, 0e0)
set_b_n_fam(lat, "sd_h", MultipoleIndex.sextupole, 0e0)
M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
stable, nu, xi = lo.compute_nu_xi(desc, tpsa_order, M)
print("\n  nu = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))
print("  xi = [{:7.5}, {:7.5}]".format(xi[X_], xi[Y_]))

match_straight(lat, loc, prm_list, bounds, Twiss0, Twiss1, Twiss1_design, weight)

if not False:
    plt.show()
print("\nPlots saved as: before.png & after.png")
