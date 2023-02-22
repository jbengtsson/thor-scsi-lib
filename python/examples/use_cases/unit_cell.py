"""Use Case:
     Optimisation of unit cell.
   Input:
     Lattice with one unit cell.
   Constraints:
     horizontal emittance,
     cell tune,
     minimize linear chromaticity,
     linear momentum compaction.
   Parameters:
     dipole with bend angle,
     dipole with reverse bend angle & focusing gradient,
     and quadrupole with defocusing gradient.
"""
import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
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

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
    radiate as rad

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
    sextupole  = 3


def plt_Twiss(ds, file_name, title):
    # Turn interactive mode off.
    plt.ioff()

    fig, (gr_1, gr_2) = plt.subplots(2)

    gr_1.set_title(title)
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel(r"$\beta_{x,y}$ [m]")
    gr_1.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), "b-",
              label=r"$\beta_x$")
    gr_1.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), "r-",
              label=r"$\beta_y$")
    gr_1.legend()

    gr_2.set_xlabel("s [m]")
    gr_2.set_ylabel(r"$\eta_x$ [m]")
    gr_2.plot(ds.s, ds.dispersion.sel(phase_coordinate="x"), label=r"$\eta_x$")
    fig.tight_layout()
    plt.savefig(file_name)


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


def compute_periodic_solution(lat, model_state, prt):
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
    if prt:
        print("\nM:\n", mat2txt(M.jacobian()))
    res= cs.compute_Twiss_A(A)
    Twiss = res[:3]
    if prt:
        prt_Twiss("\nTwiss:\n", Twiss)
    ds = lo.compute_Twiss_along_lattice(n_dof, lat, model_state, A=A, desc=desc)

    return M, A, ds

# check that a is a vector, as an example
# _, = a.shape
# compare to
# np.atleast_2d

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


def get_phi_elem(lat, fam_name, n_kid):
    elem = lat.find(fam_name, n_kid)
    return elem.get_length() * elem.get_curvature() * 180e0 / np.pi


def set_phi_elem(lat, fam_name, kid_num, phi):
    b = lat.find(fam_name, kid_num)
    L = b.get_length()
    h = phi * np.pi / (L * 180e0)
    b.set_curvature(h)


def set_phi_fam(lat, fam_name, phi):
    b = lat.find(fam_name, 0)
    L = b.get_length()
    h = phi * np.pi / (L * 180e0)
    for b in lat.elements_with_name(fam_name):
        b.set_curvature(h)


def compute_phi(lat):
    """Compute total bend angle.
    """
    prt = False
    phi = 0e0
    for k in range(len(lat)):
        if type(lat[k]) == tslib.Bending:
            dphi = get_phi_elem(lat, lat[k].name, 0)
            phi += dphi
            if prt:
                print("{:8s} {:5.3f} {:6.3f}".
                      format(lat[k].name, lat[k].get_length(), dphi))
    return phi


def get_prm(lat, prm_list):
    # Dictionary of parameter types and corresponding get functions.
    how_to_get_prm = {
        "dphi": get_phi_elem,
        "b_2":  get_b_2_elem,
        "L":    get_L_elem
    }

    prms = []
    for k in range(len(prm_list)):
        prms.append(how_to_get_prm[prm_list[k][1]](lat, prm_list[k][0], 0))
    return np.array(prms)


def set_prm(lat, prm_list, prms):
    # Dictionary of parameter types and corresponding set functions.
    how_to_set_prm = {
        "dphi": set_phi_fam,
        "b_2":  set_b_2_fam,
        "L":    set_L_fam
    }

    for k in range(len(prm_list)):
        how_to_set_prm[prm_list[k][1]](lat, prm_list[k][0], prms[k])


@dataclass
class MatchingState:
    n_iter: float = 0
    #: :math:`\chiˆ2`
    chi_2: float = np.nan


def set_phi_cell(name, phi):
    print("set_phi_cell: {name:8s} {phi:5.3f}")


# Global additional variables needed for optimisation function.
n_iter    = 0
chi_2_min = 1e30
n_iter_min  = 0
prms_min  = []

def opt_unit_cell(lat, prm_list, C, bounds, phi_des, eps_x_des, nu_des, weight,
                  rbend, phi):
    """Use Case: optimise unit cell.
    """
    def compute_chi_2(
            lat: tslib.Accelerator, model_state: tslib.ConfigType,
            eps_x_des: float, nu_des, prms
    ) -> Tuple[float, float]:
        """Computes weighted sum square
        """

        global n_iter
        global chi_2_min
        global n_iter_min
        global prms_min

        M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
        stable = lo.check_if_stable_3D(M.jacobian())

        if stable[0] and stable[1] and stable[2]:
            alpha_c = M[ct_].get(lo.ind_1(delta_))/C
            _, nu, xi = lo.compute_nu_xi(desc, tpsa_order, M)
            U_0, J, tau, eps = \
                rad.compute_radiation(lat, model_state, 2.5e9, 1e-15, desc=desc)

            n = len(weight)
            dchi_2 = np.zeros(n, float)
            dchi_2[0] = weight[0] * np.sum((eps[X_] - eps_x_des) ** 2)
            dchi_2[1] = weight[1] * np.sum((nu - nu_des) ** 2)
            dchi_2[2] = weight[2] * np.sum(xi ** 2)
            dchi_2[3] = weight[3] * alpha_c ** 2
            chi_2 = 0e0
            for k in range(n):
                chi_2 += dchi_2[k]
        else:
            chi_2 = 1e30

        if chi_2 < chi_2_min:
            chi_2_min = chi_2
            n_iter_min = n_iter
            prms_min = prms
            print(f"\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n\n  prms   =", prms)
            if stable[0] and stable[1] and stable[2]:
                print(f"  dchi_2 =", dchi_2)
            print(f"\n  phi     = {compute_phi(lat):5.3f}")
            print(f"""  bb_h    = {get_phi_elem(lat, "bb_h", 0):5.3f}""")
            print(f"  dphi    = {get_phi_elem(lat, rbend, 0):7.5f}")
            if stable[0] and stable[1] and stable[2]:
                print(f"  eps_x   = {1e12*eps[X_]:5.3f}")
                print(f"  nu      = [{nu[X_]:7.5f}, {nu[Y_]:7.5f}]")
                print(f"  xi      = [{xi[X_]:5.3f}, {xi[Y_]:5.3f}]")
                print(f"  alpha_c = {alpha_c:9.3e}")
            else:
                print("\n Unstable.")
        return chi_2

    def f_unit_cell(prms):
        global n_iter
        global chi_2_min
        global prms_min

        n_iter += 1
        set_prm(lat, prm_list, prms)
        dphi = phi/2e0 - prms[0]
        set_phi_fam(lat, rbend, dphi)
        chi_2 = compute_chi_2(lat, model_state, eps_x_des, nu_des, prms)
        return chi_2

    def prt_result(f_unit_cell, prm_list, prms0, minimum):
        global n_iter
        global chi_2_min
        global n_iter_min
        global prms_min

        M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
        stable = lo.check_if_stable_3D(M.jacobian())
        if not stable[0] or not stable[1] or not stable[2]:
            print("\nCell unstable:")
            print("  prms_min =", prms_min)
            print("  Setting prs_min for plot.")
            set_prm(lat, prm_list, prms_min)
        # Compute new Twiss parameters along lattice.
        M, A, data = compute_periodic_solution(lat, model_state, False)
        plt_Twiss(data, "after.png", "Linear Optics - After")

        print("\nInitial parameter values:")
        n = n_iter_min
        n_iter = 0
        chi_2_min = 1e30
        f_unit_cell(prms0)
        print("\nFinal parameter values:")
        n_iter = n - 1
        chi_2_min = 1e30
        f_unit_cell(minimum["x"])
        print("\n Minimum:\n", minimum)

    max_iter = 1000
    f_tol = 1e-4
    x_tol = 1e-4

    print("\nopt_unit_cell:\n")
    # Initialise parameters.
    prms1 = prms0 = get_prm(lat, prm_list)
    print("\nprms = ", prms1)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.
    minimum = optimize.minimize(
        f_unit_cell,
        prms1,
        method="Nelder-Mead",
        bounds=bounds,
        options={"xtol": x_tol, "ftol": f_tol, "maxiter": max_iter},
    )

    prt_result(f_unit_cell, prm_list, prms0, minimum)


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

if False:
    print(compute_phi(lat))
    set_phi_elem(lat, "br", 0, -0.20)
    print(compute_phi(lat))
    raise Exception("")

# Compute Twiss parameters along lattice.
M, A, data = compute_periodic_solution(lat, model_state, True)
plt_Twiss(data, "before.png", "Linear Optics - Before")
if False:
    plt.show()

# Zero sextopoles.
print("\nZeroing sextupoles.")
set_b_n_fam(lat, "sf_h", MultipoleIndex.sextupole, 0e0)
set_b_n_fam(lat, "sd_h", MultipoleIndex.sextupole, 0e0)
# Compute linear chromaticity.
M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
stable, nu, xi = lo.compute_nu_xi(desc, tpsa_order, M)
# Compute radiation properties.
U_0, J, tau, eps = \
    rad.compute_radiation(lat, model_state, 2.5e9, 1e-15, desc=desc)

C = rad.compute_circ(lat)

print("\n  C [m]          = {:5.3f}".format(C))
print("  phi [deg]      = {:5.3f}".format(compute_phi(lat)))
print("  nu             = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))
print("  eps_x [pm.rad] = {:3.1f}".format(1e12*eps[X_]))
print("  xi             = [{:7.5f}, {:7.5f}]".format(xi[X_], xi[Y_]))

# Design parameters.
phi   = 4.0                   # Total bend angle.
eps_x = 100e-12               # Horizontal emittance.
nu    = np.array([0.4, 0.1])  # Cell tune.

# Weights.
weight = np.array([
    1e15,                     # eps_x.
    1e0,                      # nu.
    1e-6,                     # xi.
    1e0                       # alpha_c
])

# Parameter family names & type.
prm_list = [
    ("bb_h", "dphi"),
    ("br",   "b_2"),
    ("qd", "b_2")
]

# Max parameter range.
b_2_max = 12.0

bounds = [
    (1.5,      2.5),          # dphi.
    (0.0,      b_2_max),      # br b_2.
    (-b_2_max, 0.0),          # qd b_2.
]

opt_unit_cell(lat, prm_list, C, bounds, phi, eps_x, nu, weight, "br", phi)

if not False:
    plt.show()
print("\nPlots saved as: before.png & after.png")
