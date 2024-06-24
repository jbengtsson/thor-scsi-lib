"""Use Case:
     Module for implementing & optimising a unit cell.
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

eps_x_des    = 125e-12
nu_uc        = [0.29, 0.18]


def opt_uc(lat_prop, prm_list, weight, b1_list, phi_lat, eps_x_des, nu_uc):
    """Use Case: optimise super period.
    """
    chi_2_min = 1e30
    eta       = np.nan
    n_iter    = 0
    file_name = "opt_uc.txt"

    def prt_iter(prm, chi_2, nu, xi):
        nonlocal n_iter

        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        rb_1     = "qf1"
        phi_bend = "b1_5"

        phi = lat_prop.compute_phi_lat()
        phi_b1 = compute_phi_bend(lat_prop, b1_list)
        phi_rb_1 = lat_prop.get_phi_elem(rb_1, 0)
        phi_bend = lat_prop.get_phi_elem(phi_bend, 0)

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("    eps_x [pm.rad] = {:5.3f} [{:5.3f}]".
              format(1e12*lat_prop._eps[ind.X], 1e12*eps_x_des))
        print("    nu_uc          =  [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu[ind.X], nu[ind.Y], nu_uc[ind.X], nu_uc[ind.Y]))
        print("    xi             =  [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        print("\n    phi_sp         = {:8.5f}".format(phi))
        print("    C [m]          = {:8.5f}".format(lat_prop.compute_circ()))
        print("\n    phi_b1         = {:8.5f}".format(phi_b1))
        print("    phi_rb_1       = {:8.5f}".format(phi_rb_1))
        print("    phi_bend       = {:8.5f}".format(phi_bend))
        prm_list.prt_prm(prm)

    def compute_chi_2(nu, xi):
        prt = not False

        dchi_2 = weight[0]*(lat_prop._eps[ind.X]-eps_x_des)**2
        chi_2 = dchi_2
        if prt:
            print("\n  dchi2(eps_x)    = {:10.3e}".format(dchi_2))

        dchi_2 = weight[1]*(nu[ind.X]-nu_uc[ind.X])**2            
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_x)  = {:10.3e}".format(dchi_2))

        dchi_2 = weight[2]*(nu[ind.Y]-nu_uc[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_y)  = {:10.3e}".format(dchi_2))

        dchi_2 = weight[3]*(xi[ind.X]**2+xi[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(xi)       = {:10.3e}".format(dchi_2))

        return chi_2

    def f_sp(prm):
        nonlocal chi_2_min, n_iter

        n_iter += 1
        prm_list.set_prm(prm)
        phi_lat.set_phi_lat()

        # Compute the beam dynamics properties.
        try:
            if not lat_prop.comp_per_sol():
                print("\nf_sp - comp_per_sol: unstable")
                raise ValueError

            if not lat_prop.compute_radiation():
                print("\nf_sp - compute_radiation: unstable")
                raise ValueError

            # Compute linear chromaticity.
            M = lo.compute_map(
                lat_prop._lattice, lat_prop._model_state,
                desc=lat_prop._desc, tpsa_order=2)
            stable, nu, xi = \
                lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
            if not stable:
                print("\nf_sp - compute_nu_xi: unstable")
                raise ValueError
        except ValueError:
            return 1e30
        else:
            # _, _, _, nu = lat_prop.get_Twiss(-1)

            chi_2 = compute_chi_2(nu, xi)
            if chi_2 < chi_2_min:
                prt_iter(prm, chi_2, nu, xi)
                pc.prt_lat(lat_prop, "opt_uc.txt", prm_list, phi_lat=phi_lat)
                chi_2_min = min(chi_2, chi_2_min)
            return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4
    g_tol    = 1e-5

    prm, bounds = prm_list.get_prm()

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
        # options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
        options={"gtol": g_tol, "maxiter": max_iter}
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
lat_name = "max_4u_f_4"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol - unstable")
        raise ValueError

    # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation - unstable")
        raise ValueError
except ValueError:
    exit

# Weights.
weight = np.array([
    1e15, # eps_x.
    1e-1, # nu_uc_x.
    1e-1, # nu_uc_y.
    1e-7  # xi.
])

b1_list = [
    "b1_0", "b1_1", "b1_2", "b1_3", "b1_4", "b1_5"]

b1_bend = pc.bend_class(lat_prop, b1_list, phi_max, b_2_max)

prms = [
    ("qf1",      "b_2"),

    ("b1_0",     "b_2"),
    ("b1_1",     "b_2"),
    ("b1_2",     "b_2"),
    ("b1_3",     "b_2"),
    ("b1_4",     "b_2"),
    ("b1_5",     "b_2"),

    # ("phi_bend", b1_bend),

    ("b1_0",     "phi"),
    ("b1_1",     "phi"),
    ("b1_2",     "phi"),
    ("b1_3",     "phi"),
    ("b1_4",     "phi"),
    # ("b1_5",     "phi"),

    ("qf1",     "phi"),
]

prm_list = pc.prm_class(lat_prop, prms, b_2_max)

# To maintain the total bend angle.
phi_lat = pc.phi_lat_class(lat_prop, 2, "b1_5")

opt_uc(lat_prop, prm_list, weight, b1_list, phi_lat, eps_x_des, nu_uc)
