"""Use Case:
     Parametric scans/evaluations for a unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

from dataclasses import dataclass
import os
from typing import Tuple

import numpy as np
from scipy import optimize as opt

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, prm_class as pc
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()

L_min        = 0.10
L_max        = 0.50
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0


def opt_sp(Lat_prop, prm_list, weight):
    """Use Case: optimise super period.
    """

    loc       = lat_prop._lattice.find("sf_h", 2).index
    # Dispersion at the unit cell end.
    eta_uc_x  = 0.06255
    # Beta functions at the unit cell end.
    beta_uc   = [3.65614, 3.68868]
    nu_sp     = [2.110, 0.814]

    chi_2_min = 1e30
    eta       = np.nan
    n_iter    = 0
    file_name = "opt_sp.txt"

    def prt_iter(prm, chi_2, eta, beta, nu, xi):
        nonlocal n_iter

        phi = lat_prop.compute_phi_lat()
        # Twiss functions at end of super period.
        _, _, beta_id, _ = lat_prop.get_Twiss(-1)

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:5.3f}".format(1e12*lat_prop._eps[ind.X]))
        print("  U_0 [keV]      = {:5.3f}".format(1e-3*lat_prop._U_0))
        print("  eta_x [m]      =  {:10.3e} ({:9.3e})".
              format(eta[ind.x], eta_uc_x))
        print("  beta_match [m] =  [{:5.3f}, {:5.3f}] ([{:5.3f}, {:5.3f}])".
              format(beta[ind.X], beta[ind.Y],
                     beta_uc[ind.X], beta_uc[ind.Y]))
        print("  beta_id        =  [{:5.3f}, {:5.3f}]".
              format(beta_id[ind.X], beta_id[ind.Y]))
        print("  nu_sp          =  [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu[ind.X], nu[ind.Y], nu_sp[ind.X], nu_sp[ind.Y]))
        print("  xi             =  [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        print("  phi_sp         =  {:8.5f}".format(phi))
        prm_list.prt_prm(prm)

    def compute_chi_2(eta, beta, nu, xi):
        nonlocal loc, eta_uc_x

        prt = not False

        dchi_2 = weight[0]*lat_prop._eps[ind.X]**2
        chi_2 = dchi_2
        if prt:
            print("\n  dchi2(eps_x)    = {:10.3e}".format(dchi_2))

        dchi_2 = weight[1]*lat_prop._U_0**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(U_0)      = {:10.3e}".format(dchi_2))

        dchi_2 = weight[2]*(eta[ind.x]-eta_uc_x)**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(eta_uc_x) = {:10.3e}".format(dchi_2))

        dchi_2 = \
            weight[3]*(
                (beta[ind.X]-beta_uc[ind.X])**2
                +(beta[ind.Y]-beta_uc[ind.Y])**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(beta)     = {:10.3e}".format(dchi_2))

        dchi_2 = \
            weight[4]*(
                (nu[ind.X]-nu_sp[ind.X])**2
                +(nu[ind.Y]-nu_sp[ind.Y])**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_sp)    = {:10.3e}".format(dchi_2))

        dchi_2 = weight[5]*(xi[ind.X]**2+xi[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(xi)       = {:10.3e}".format(dchi_2))

        return chi_2

    def f_sp(prm):
        nonlocal chi_2_min, n_iter

        n_iter += 1
        prm_list.set_prm(lat_prop, prm)

        # Compute the beam dynamics properties.
        try:
            if not lat_prop.comp_per_sol():
                print("\nf_sp - comp_per_sol: unstable")
                raise ValueError

            if not lat_prop.compute_radiation():
                print("\nf_sp - compute_radiation: unstable")
                raise ValueError

            eta, _, beta, _ = lat_prop.get_Twiss(loc)
            _, _, _, nu = lat_prop.get_Twiss(-1)

            M = lo.compute_map(
                lat_prop._lattice, lat_prop._model_state,
                desc=lat_prop._desc, tpsa_order=2)
            stable, _, xi = \
                lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
            if not stable:
                raise ValueError
        except ValueError:
            print("\nf_sp - compute_nu_xi: unstable")
            return 1e30
        else:
            chi_2 = compute_chi_2(eta, beta, nu, xi)
            if chi_2 < chi_2_min:
                prt_iter(prm, chi_2, eta, beta, nu, xi)
                pc.prt_lat(lat_prop, "opt_sp.txt", prm_list)
                chi_2_min = min(chi_2, chi_2_min)
            return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    prm, bounds = prm_list.get_prm(lat_prop)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol, xtol
    # CG     gtol
    minimum = opt.minimize(
        f_sp,
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
lat_name = "max_4u_sp_1"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.get_types()

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncompute_chi_2 - comp_per_sol: unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop.prt_lat_param()
    lat_prop.prt_M()

try:
    # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation - unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop.prt_rad()
    lat_prop.prt_M_rad()

if False:
    lat_prop.prt_lat_param()
    lat_prop.prt_Twiss("max_4u_uc.txt")

# Weights.
weight = np.array([
    1e13,  # eps_x.
    1e-13, # U_0.
    1e2,   # eta_uc_x.
    1e-2,  # beta_uc.
    0e-1,  # nu_sp.
    1e-7   # xi.
])

dip_list = ["b1_1", "b1_2", "b1_3", "b1_4", "b1_5"]
b1 = \
    pc.bend_prm_class(
        lat_prop, "b1_0", dip_list, not False, phi_max, b_2_bend_max)

dip_list  = [
    "b2u_6", "b2u_5", "b2u_4", "b2u_3", "b2u_2", "b2u_1", "b2d_1", "b2d_2",
    "b2d_3", "b2d_4", "b2d_5"
]
b2 = pc.bend_prm_class(lat_prop, "b2_0", dip_list, True, phi_max, b_2_bend_max)

rb_list  = ["qf1"]
dip_list = ["b1_0"]
rb1 = pc.rev_bend_prm_class(lat_prop, rb_list, dip_list, phi_max)

rb_list  = ["qf1_e"]
dip_list = ["b2_0"]
rb2 = pc.rev_bend_prm_class(lat_prop, rb_list, dip_list, phi_max)

prm_list = [
    ("qf1",    "b_2"),

    ("qf1_e",    "b_2"),
    ("qd",       "b_2"),
    ("qf2",      "b_2"),
    ("b1_0",     "L_b"),
    ("b2_0",     "L_b"),
    ("bend",     b1),
    ("bend",     b2),
    ("rev_bend", rb1),
    ("rev_bend", rb2)
]

prm_list = pc.prm_class(lat_prop, prm_list, b_2_max)

if not False:
    opt_sp(lat_prop, prm_list, weight)

if False:
    dip_list = [bend]
    dip_list.extend(bend_list)
    prt_lat(lat_prop, "opt_sp.txt", dip_list)
