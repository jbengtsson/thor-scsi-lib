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

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    get_set_mpole as gs, index_class as ind, courant_snyder as cs

from thor_scsi.utils.output import mat2txt, vec2txt


# check that a is a vector, as an example
# _, = a.shape
# compare to
# np.atleast_2d


def compute_phi_bend(lat_prop, bend_list):
    phi_b = 0e0
    for k in range(len(bend_list)):
        phi_b += lat_prop.get_phi_elem(bend_list[k], 0)
    return phi_b


def prt_lat(lat_prop, file_name, prm_list, bend):
    outf = open(file_name, 'w')

    def prt_drift(name):
        L = lat_prop.get_L_elem(name, 0)
        print(("{:4s}: Drift, L = {:7.5f};").format(name, L), file=outf)

    def prt_bend(name):
        L = lat_prop.get_L_elem(name, 0)
        phi = lat_prop.get_phi_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:4s}: Bending, L = {:7.5f}, T = {:8.5f}, K = {:8.5f}"
               ", T1 = 0.0, T2 = 0.0,\n      N = n_bend;")
              .format(name, L, phi, b_2), file=outf)

    def prt_quad(name):
        L = lat_prop.get_L_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:4s}: Quadrupole, L = {:7.5f}, K = {:8.5f}, N = n_quad;")
              .format(name, L, b_2), file=outf)

    for k in range(len(prm_list)):
        elem = lat_prop._lattice.find(prm_list[k][0], 0)
        if  type(elem) == ts.Drift:
            prt_drift(prm_list[k][0])
        elif  type(elem) == ts.Bending:
            prt_bend(prm_list[k][0])
        elif  type(elem) == ts.Quadrupole:
            prt_quad(prm_list[k][0])

    prt_bend(bend)


@dataclass
class MatchingState:
    n_iter: float = 0
    #: :math:`\chiˆ2`
    chi_2: float = np.nan


def match_straight(
        lat_prop, prm_list, phi_b, bend, bend_list, bounds, Twiss_0, Twiss_1,
        weight):

    chi_2_min = 1e30
    n_iter    = 0
    A0        = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)

    def get_prm():
        # Dictionary of parameter types and corresponding get functions.
        how_to_get_prm = \
            {"L":   lat_prop.get_L_elem,
             "phi": lat_prop.get_phi_elem,
             "b_2": lat_prop.get_b_2_elem}
        prms = []
        for k in range(len(prm_list)):
            prms.append(how_to_get_prm[prm_list[k][1]](prm_list[k][0], 0))
        return prms

    def set_prm(prms):
        # Dictionary of parameter types and corresponding set functions.
        how_to_set_prm = {
            "L":   lat_prop.set_L_fam,
            "phi": lat_prop.set_phi_fam,
            "b_2": lat_prop.set_b_2_fam}
        for k in range(len(prm_list)):
            how_to_set_prm[prm_list[k][1]](prm_list[k][0], prms[k])
        lat_prop.set_phi_fam(bend, phi_b-compute_phi_bend(lat_prop, bend_list))
        
    def prt_prms(prms):
        n_prt = 5
        print("\n  prms =")
        for k in range(len(prms)):
            print(" {:12.5e}".format(prms[k]), end="")
            if k % n_prt == 4:
                print()
        if k % n_prt != 4:
            print()

    def prt_iter(prms, chi_2, Twiss_k):
        phi_dip = lat_prop.get_phi_elem(bend, 0)
        phi_bend = compute_phi_bend(lat_prop, bend_list) + phi_dip
        print("\n{:4d} chi_2 = {:9.3e}".format(n_iter, chi_2))
        print("  phi_b   = {:7.5f}".format(phi_bend))
        print("  phi_{:s} = {:7.5f}\n".format(bend, phi_dip))
        lat_prop.prt_Twiss_param(Twiss_k[:3])
        prt_prms(prms)

    def compute_chi_2_Twiss(Twiss_k):
        dchi_2 = 0e0
        for j in range(3):
            for k in range(2):
                dchi_2 += weight[j] * ((Twiss_1[j][k] - Twiss_k[j][k]) ** 2)
        return dchi_2

    def compute_chi_2(prms):
        nonlocal A0

        A1 = _copy.copy(A0)
        lat_prop._lattice.propagate(
            lat_prop._model_state, A1, 0, len(lat_prop._lattice))
        Twiss_k = cs.compute_Twiss_A(A1.jacobian())
        chi_2 = compute_chi_2_Twiss(Twiss_k)
        return chi_2, Twiss_k

    def f_match(prms):
        nonlocal chi_2_min
        nonlocal n_iter

        n_iter += 1
        set_prm(prms)
        chi_2, Twiss_1 = compute_chi_2(prms)
        if chi_2 < chi_2_min:
            prt_iter(prms, chi_2, Twiss_1)
            prt_lat(lat_prop, "match_lat_k.txt", prm_list, bend)
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

    prms = get_prm()

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol
    # CG     gtol
    minimum = optimize.minimize(
        f_match,
        prms,
        method="Powell",
        bounds=bounds,
        options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
    )

    print("\n".join(minimum))

    prms = minimum["x"]
    set_prm(prms)

    lat_prop._data = \
                  lo.compute_Twiss_along_lattice(
                      lat_prop._n_dof, lat_prop._lattice, lat_prop._model_state,
                      A=A0, desc=lat_prop._desc, mapping=lat_prop._named_index)

    prt_lat(lat_prop, "match_lat.txt", prm_list, bend)
    lat_prop.prt_Twiss("match_twiss.txt")

    lat_prop.plt_Twiss("match_twiss.png", True)


def chk_lat(lat_name, lat_prop, Twiss_0):
    n_dof = 2
    A = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    A.set_zero()
    # Use *-operator to unpack the list of arguments.
    A_7x7 = np.zeros((7, 7))
    A_7x7[:6, :6] = cs.compute_A(*Twiss_0[:3])
    A.set_jacobian(A_7x7)

    lat_prop._data = \
        lo.compute_Twiss_along_lattice(
            n_dof, lat_prop._lattice, lat_prop._model_state, A=A,
            desc=lat_prop._desc, mapping=lat_prop._named_index)

    lat_prop.prt_Twiss("chk_lat_twiss.txt")

    Twiss_1 = cs.compute_Twiss_A(A.jacobian())[:3]

    print("\nchk_lat - Entrance:")
    lat_prop.prt_Twiss_param(Twiss_0)
    print("\nchk_lat - Exit:")
    lat_prop.prt_Twiss_param(Twiss_1)

    lat_prop.plt_Twiss("chk_lat.png", not False)


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
      format(lat_prop.compute_phi_lat()))

# Entrance.
eta   = np.array([0.02273, 0.0])
alpha = np.array([0.0, 0.0])
beta  = np.array([0.92011, 9.72555])
Twiss_0 = eta, alpha, beta

if False:
    chk_lat(lat_name, lat_prop, Twiss_0)
    assert False

# Desired at exit.
eta   = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta  = np.array([9.2, 2.0])
Twiss_1 = eta, alpha, beta

# Weights: eta, alpha, beta.
weight = np.array([
    1e8, # eta.
    1e4, # alpha.
    0*1e0  # beta.
])

# End dipole.

phi_b = 1.5

bend = "ds0"

bend_list = [
    "ds6", "ds5", "ds4", "ds3", "ds2", "ds1", "dm1", "dm2", "dm3", "dm4",
    "dm5"]

dip_list = [bend]
dip_list.extend(bend_list)

print("\nphi_b = {:7.5f}".format(compute_phi_bend(lat_prop, dip_list)))

# Parameters: family name & type.
prm_list = [
    ("qf1e", "b_2"),
    ("qd",   "b_2"),
    ("qf2",  "b_2"),
    # ("d5",   "L"),

    ("ds6",  "b_2"),
    ("ds5",  "b_2"),
    ("ds4",  "b_2"),
    ("ds3",  "b_2"),
    ("ds2",  "b_2"),
    ("ds1",  "b_2"),
    ("ds0",  "b_2"),
    ("dm1",  "b_2"),
    ("dm2",  "b_2"),
    ("dm3",  "b_2"),
    ("dm4",  "b_2"),
    ("dm5",  "b_2"),

    ("ds6",  "phi"),
    ("ds5",  "phi"),
    ("ds4",  "phi"),
    ("ds3",  "phi"),
    ("ds2",  "phi"),
    ("ds1",  "phi"),
    ("dm1",  "phi"),
    ("dm2",  "phi"),
    ("dm3",  "phi"),
    ("dm4",  "phi"),
    ("dm5",  "phi")
]

# Max parameter range.
L_min        = 0.1
phi_max      = 0.7
quad_b_2_max = 10.0
bend_b_2_max = 1.0

bounds = [
    ( 0.0,          quad_b_2_max),
    (-quad_b_2_max, 0.0),
    ( 0.0,          quad_b_2_max),
    # ( 0.15,         0.30),

    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),
    (-bend_b_2_max, bend_b_2_max),

    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max),
    (0.0, phi_max)
]

match_straight(
    lat_prop, prm_list, phi_b, bend, bend_list, bounds, Twiss_0, Twiss_1,
    weight)

if not False:
    plt.show()
