"""Use Case:
     Optimisation of super period.
   Input:
     Lattice with one super period.
   Constraints:
     horizontal emittance,
     cell tune,
     minimize linear chromaticity,
     linear momentum compaction.
   Parameters:
     dipole bend angle,
     focusing gradient gradient reverse bend angle,
     defocusing quadrupole,
     triplet gradients (3) & spacing (3).
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


def plot_Twiss(ds, file_name, title):
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


def print_Twiss_params(str, Twiss):
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


def print_Twiss(lat, data):
    """
    Print Twiss parameters along the lattice.
    """
    L = 0e0
    nu = np.zeros(2, dtype=float)
    print("")
    for k in range(len(data.index)):
        L += lat[k].get_length()
        nu[X_] += data.twiss.sel(plane="x", par="dnu").values[k]
        nu[Y_] += data.twiss.sel(plane="y", par="dnu").values[k]
        print("{:3d} {:8s} {:7.3f} {:8.5f} {:8.5f} {:7.5f} {:8.5f} {:8.5f} {:8.5f} {:7.5f} {:8.5f}".
              format(k, lat[k].name, L,
                     data.twiss.sel(plane="x", par="alpha").values[k],
                     data.twiss.sel(plane="x", par="beta").values[k],
                     nu[X_], data.dispersion.sel(phase_coordinate="x")[k],
                     data.twiss.sel(plane="y", par="alpha").values[k],
                     data.twiss.sel(plane="y", par="beta").values[k],
                     nu[Y_], data.dispersion.sel(phase_coordinate="y")[k]))


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

    stable, M, A = lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc)
    if stable:
        if prt:
            print("\nM:\n", mat2txt(M.jacobian()))
        res= cs.compute_Twiss_A(A)
        Twiss = res[:3]
        if prt:
             print_Twiss_params("\nTwiss:\n", Twiss)
        ds = \
            lo.compute_Twiss_along_lattice(
                n_dof, lat, model_state, A=A, desc=desc)
    else:
        ds = np.nan
        print("\ncompute_periodic_solution: unstable")

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


def get_prm(lat, param_list):
    # Dictionary of parameter types and corresponding get functions.
    how_to_get_prm = {
        "dphi":     get_phi_elem,
        "b_2":      get_b_2_elem,
        "L":        get_L_elem
    }

    prms = []
    for k in range(len(param_list)):
        print(param_list[k][0], param_list[k][1])
        prms.append(how_to_get_prm[param_list[k][1]](lat, param_list[k][0], 0))
    return np.array(prms)


def set_prm(lat, param_list, prms):
    # Dictionary of parameter types and corresponding set functions.
    how_to_set_prm = {
        "dphi":     set_phi_fam,
        "b_2":      set_b_2_fam,
        "L":        set_L_fam
    }

    for k in range(len(param_list)):
        how_to_set_prm[param_list[k][1]](lat, param_list[k][0], prms[k])


@dataclass
class MatchingState:
    n_iter: float = 0
    #: :math:`\chiˆ2`
    chi_2: float = np.nan


def set_phi_cell(name, phi):
    print("set_phi_cell: {name:8s} {phi:5.3f}")


# Global additional variables needed for optimisation function.
n_iter     = 0
chi_2_min  = 1e30
n_iter_min = 0
prms_min   = []
rbend      = ""

def opt_super_period(lat, param_list, C, bounds, phi_des, eps_x_des, nu_des,
                     beta_straight_des, weights, phi):
    """Use Case: optimise unit cell.
    """

    def compute_nu_cell(data):
        nu_cell = np.zeros(2, dtype=float)
        loc1 = lat.find("sf", 1).index
        loc2 = lat.find("sf", 3).index
        for k in range(loc1, loc2):
            nu_cell[X_] += data.twiss.sel(plane="x", par="dnu").values[k]
            nu_cell[Y_] += data.twiss.sel(plane="y", par="dnu").values[k]
        return nu_cell

    def compute_chi_2(
            lat: tslib.Accelerator, model_state: tslib.ConfigType,
            eps_x_des: float, nu_des, prms
    ) -> Tuple[float, float]:
        """Computes weightsed sum square
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

            M, A, data = compute_periodic_solution(lat, model_state, False)

            beta_straight = np.zeros(2, dtype=float)
            beta_straight[X_] = data.twiss.sel(plane="x", par="beta").values[0]
            beta_straight[Y_] = data.twiss.sel(plane="y", par="beta").values[0]
            eta_straight_x    = \
                data.dispersion.sel(phase_coordinate="x").values[0]
            nu_cell = compute_nu_cell(data)

            stable, U_0, J, tau, eps = \
                rad.compute_radiation(lat, model_state, 2.5e9, 1e-15, desc=desc)

            if stable:
                n = len(weights)
                dchi_2 = np.zeros(n, dtype=float)
                dchi_2[0] = \
                    weights["eps_x"] * np.sum((eps[X_] - eps_x_des) ** 2)
                dchi_2[1] = weights["nu_cell"] * np.sum((nu_cell - nu_des) ** 2)
                dchi_2[2] = \
                    weights["beta_straight"] \
                    * np.sum((beta_straight - beta_straight_des) ** 2)
                dchi_2[3] = weights["eta_straight_x"] * eta_straight_x ** 2
                dchi_2[4] = weights["xi"] * np.sum(xi ** 2)
                dchi_2[5] = weights["alpha_c"] / alpha_c ** 2
                chi_2 = 0e0
                for k in range(n):
                    chi_2 += dchi_2[k]
            else:
                chi_2 = 1e30
        else:
            chi_2 = 1e30

        if chi_2 < chi_2_min:
            chi_2_min = chi_2
            n_iter_min = n_iter
            prms_min = prms
            print(f"\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n\n  prms   =\n", prms)
            if stable:
                print(f"\n  dchi_2 =\n", dchi_2)
            print(f"\n  phi            = {compute_phi(lat):8.5f}")
            [b1, b2, b3] = \
                [param_list[0][0], param_list[1][0], param_list[2][0]]
            print(f"\n  {b1:5s}          = {get_phi_elem(lat, b1, 0):8.5f}")
            print(f"  {b2:5s}          = {get_phi_elem(lat, b2, 0):8.5f}")
            print(f"  {b3:5s}          = {get_phi_elem(lat, b3, 0):8.5f}")
            [b1, b2, b3, b4] = \
                [param_list[3][0], param_list[4][0], param_list[5][0],
                 param_list[6][0]]
            print(f"\n  {b1:5s}          = {get_phi_elem(lat, b1, 0):8.5f}")
            print(f"  {b2:5s}          = {get_phi_elem(lat, b2, 0):8.5f}")
            print(f"  {b3:5s}          = {get_phi_elem(lat, b3, 0):8.5f}")
            print(f"  {b4:5s}          = {get_phi_elem(lat, b4, 0):8.5f}")
            [b1, b2] = [rbend, param_list[8][0]]
            print(f"\n  {b1:5s}          = {get_phi_elem(lat, b1, 0):8.5f}")
            print(f"  {b2:5s}          = {get_phi_elem(lat, b2, 0):8.5f}")
            if stable:
                print(f"\n  eps_x          = {1e12*eps[X_]:5.3f}")
                print("  nu_cell        = [{:7.5f}, {:7.5f}]".
                      format(nu_cell[X_], nu_cell[Y_]))
                print("  beta_straight  = [{:7.5f}, {:7.5f}]".
                      format(beta_straight[X_], beta_straight[Y_]))
                print(f"  eta_straight_x = {eta_straight_x:10.3e}")
                print(f"  xi             = [{xi[X_]:5.3f}, {xi[Y_]:5.3f}]")
                print(f"  alpha_c        = {alpha_c:9.3e}")
            else:
                print("\n Unstable.")
        return chi_2

    def f_unit_cell(prms):
        global n_iter
        global chi_2_min
        global prms_min
        global rbend

        n_iter += 1
        set_prm(lat, param_list, prms)
        dphi = \
            (phi
             - 8*(prms[0] + prms[1] + prms[2])
             - 2*(prms[3] + prms[4] + prms[5] + prms[6])
             - 2*prms[8])/8e0
        set_phi_fam(lat, rbend, dphi)
        chi_2 = compute_chi_2(lat, model_state, eps_x_des, nu_des, prms)
        return chi_2

    def prt_result(f_unit_cell, param_list, prms0, minimum):
        global n_iter
        global chi_2_min
        global n_iter_min
        global prms_min

        M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
        stable = lo.check_if_stable(3, M.jacobian())
        if not stable:
            print("\nCell unstable:")
            print("  prms_min =", prms_min)
            print("  Setting prs_min for plot.")
            set_prm(lat, param_list, prms_min)
        # Compute new Twiss parameters along lattice.
        M, A, data = compute_periodic_solution(lat, model_state, False)
        plot_Twiss(data, "after.png", "Linear Optics - After")

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

        print_Twiss(lat, data)

    global rbend

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    print("\nopt_super_period:\n")
    # Initialise parameters.
    prms1 = prms0 = get_prm(lat, param_list)
    rbend = param_list[7][0]
    print("\nrbend = ", rbend)
    print("\nprms = ", prms1)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol
    # CG     gtol
    minimum = optimize.minimize(
        f_unit_cell,
        prms1,
        method="CG",
        # bounds=bounds,
        # options={"xtol": x_tol, "maxiter": max_iter},
        options={"gtol": f_tol, "maxiter": max_iter},
    )

    prt_result(f_unit_cell, param_list, prms0, minimum)


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


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB")
# t_file = os.path.join(t_dir, 'b3_tst.lat')
# t_file = os.path.join(t_dir, "b3_sf_40Grad_JB.lat")
# t_file = os.path.join(t_dir, "b3_sfsf_4Quads_unitcell.lat")
# t_file = os.path.join(t_dir, "b3_sfsf4Q_tracy_jb.lat")
t_file = os.path.join(t_dir, "b3_sfsf4Q_tracy_jb_2.lat")

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

if False:
    # Set RF cavity phase for negative alpha_c.
    cav = lat.find("cav", 0)
    cav.set_phase(180.0)
    print("\nRF Cavity:")
    print(f"  f [MHz]   = {1e-6*cav.get_frequency():3.1f}")
    print(f"  V [MV]    = {1e-6*cav.get_voltage():3.1f}")
    print(f"  h         = {cav.get_harmonic_number():1d}")
    print(f"  phi [deg] = {cav.get_phase():3.1f}")

if False:
    print(compute_phi(lat))
    set_phi_elem(lat, "br", 0, -0.20)
    print(compute_phi(lat))
    raise Exception("")

# Compute Twiss parameters along lattice.
M, A, data = compute_periodic_solution(lat, model_state, True)
plot_Twiss(data, "before.png", "Linear Optics - Before")
if False:
    print_Twiss(lat, data)
if False:
    plt.show()

# Zero sextopoles.
print("\nZeroing sextupoles.")
set_b_n_fam(lat, "sf", MultipoleIndex.sextupole, 0e0)
set_b_n_fam(lat, "sd", MultipoleIndex.sextupole, 0e0)
# Compute linear chromaticity.
M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
stable, nu, xi = lo.compute_nu_xi(desc, tpsa_order, M)
# Compute radiation properties.
stable, U_0, J, tau, eps = \
    rad.compute_radiation(lat, model_state, 2.5e9, 1e-15, desc=desc)

C = rad.compute_circ(lat)

print("\n  C [m]          = {:5.3f}".format(C))
print("  phi [deg]      = {:5.3f}".format(compute_phi(lat)))
print("  nu             = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))
print("  eps_x [pm.rad] = {:3.1f}".format(1e12*eps[X_]))
print("  xi             = [{:7.5f}, {:7.5f}]".format(xi[X_], xi[Y_]))

# Design parameters.
phi   = 22.5                  # Total bend angle.
eps_x = 0*100e-12               # Horizontal emittance.
nu    = np.array([0.4, 0.1])  # Cell tune.
beta  = np.array([3.0, 3.0])  # Beta functions at the centre of the straight.

# Weights.
weights = {
    "eps_x":          1e2*1e15,
    "nu_cell":        1e0,
    "beta_straight":  1e-5,
    "eta_straight_x": 1e3,
    "xi":             1e-6,
    "alpha_c":        0*1e-14
}

# Parameter family names & type.
param_list = [
    ("bb_h1", "dphi"),        # Dipole.
    ("bb_h2", "dphi"),        # Dipole.
    ("bb_h3", "dphi"),        # Dipole.

    ("mbb1",  "dphi"),        # Dipole.
    ("mbb2",  "dphi"),        # Dipole.
    ("mbb3",  "dphi"),        # Dipole.
    ("mbb4",  "dphi"),        # Dipole.

#    ("br",   "not used"),     # Focusing reverse bend;
    ("br",    "b_2"),          # used to set total bend angle.

    ("mbr",   "dphi"),        # Focusing reverse bend.
    ("mbr",   "b_2"),


    ("qd",    "b_2"),          # Defocusing quadrupole.
    ("mqd",   "b_2"),          # Defocusing quadrupole.

    ("uq1",   "b_2"),          # Triplet for matching section.
    ("uq2",   "b_2"),
    ("uq3",   "b_2"),
    ("uq4",   "b_2"),

    ("ul1",   "L"),
    ("ul2",   "L"),
    ("ul3",   "L"),
    ("ul4",   "L")
]

# Max parameter range.
b_2_max = 12.0
L_min   = 0.05

bounds = [
    (0.5,      7.0),          # bb phi.
    (0.5,      7.0),          # bb phi.
    (0.5,      7.0),          # bb phi.

    (0.5,      7.0),          # mbb phi.
    (0.5,      7.0),          # mbb phi.
    (0.5,      7.0),          # mbb phi.
    (0.5,      7.0),          # mbb phi.

    # (-0.5,     1.0),          # br phi.
    ( 0.0,     b_2_max),      # br b_2.

    (-1.0,     0.0),          # mbr phi.
    (0.0,      b_2_max),      # mbr b_2.

    (-b_2_max, 0.0),          # qd b_2.
    (-b_2_max, 0.0),          # mqd b_2.

    ( 0.0,     b_2_max),      # uq1 b_2.
    (-b_2_max, 0.0),          # uq2 b_2.
    ( 0.0,     b_2_max),      # uq3 b_2.
    (-b_2_max, 0.0),          # uq4 b_2.

    ( L_min,   0.5),          # ul1 L.
    ( L_min,   0.5),          # ul2 L.
    ( L_min,   0.5),          # ul3 L.
    ( L_min,   0.5)           # ul4 L.
]

opt_super_period(lat, param_list, C, bounds, phi, eps_x, nu, beta, weights, phi)

if not False:
    plt.show()
print("\nPlots saved as: before.png & after.png")
