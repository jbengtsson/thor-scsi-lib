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


from thor_scsi.utils import lattice_properties as lp, index_class as ind


ind = ind.index_class()


def set_phi_rb(lat_prop, phi_rb):
    # Optimum reverse bend angle is:
    #   phi_rb = -0.37
    phi_b  = 3.0
    b_name = ["b0", "b1", "b2", "b3", "b4", "b5"]
    b_scl = np.array(
        [1.094181/phi_b, 0.151199/phi_b, 0.151101/phi_b, 0.101861/phi_b,
         0.001569/phi_b, 0.000089/phi_b]
    )

    dphi = 3.0 - 2.0*phi_rb

    lat_prop.set_phi_fam("qf", phi_rb, True)
    for k in range(len(b_scl)):
        lat_prop.set_phi_fam(b_name[k], b_scl[k]*dphi, True)


# Global additional variables needed for optimisation function.
phi_b      = np.nan
bend_list  = ["b1", "b2", "b3", "b4", "b5"]
chi_2_min  = 1e30
n_iter     = 0


def opt_var_bend_radius(lat_prop):
    """Use Case: optimise unit cell.
    """

    def set_phi_var_bend_rad_rb(phi):
        global bend_list

        phi_tot = 0e0
        n = len(bend_list)
        for k in range(n-1):
            lat_prop.set_phi_fam(bend_list[k], phi[k], False)
            phi_tot += phi[k]
        lat_prop.set_phi_fam("b0", phi_b-phi_tot, False)

    def prt_iter(phi, chi_2):
        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:14.10f}".format(1e12*lat_prop._eps[ind.X]))
        print("  phi_tot       = {:10.3e}".format(lat_prop.compute_phi()))
        print("  phi           = ", end="")
        print("[{:7.5f}, ".format(lat_prop.get_phi_elem("b0", 0)), end="")
        for k in range(1, len(phi)):
            if k < len(phi)-1:
                print("{:7.5f}, ".format(phi[k]), end="")
            else:
                print("{:7.5f}]".format(phi[k]))

    def compute_chi_2(phi):
        global n_iter

        n_iter += 1
        set_phi_var_bend_rad_rb(phi)
        stable = lat_prop.comp_per_sol()
        if stable:
            stable = lat_prop.compute_radiation()
            if stable:
                chi_2 = lat_prop._eps[ind.X]**2
            else:
                print("\ncompute_chi_2 - unstable")
                chi_2 = 1e30
        else:
            print("\ncompute_chi_2 - unstable")
            chi_2 = 1e30
        return chi_2

    def f_super_per(phi):
        global chi_2_min
        global n_iter

        chi_2 = compute_chi_2(phi)
        if chi_2 < chi_2_min:
            prt_iter(phi, chi_2)
            chi_2_min = min(chi_2, chi_2_min)
        return chi_2

    dhi_max  = 1e-4
    max_iter = 1000
    f_tol    = 1e-30
    x_tol    = 1e-30

    print("\nopt_var_bend_radius:")

    phi_b = lat_prop.get_phi_elem("b0", 0)
    n_phi = len(bend_list)
    phi = []
    for k in range(n_phi):
        dphi = lat_prop.get_phi_elem(bend_list[k], 0)
        phi.append(dphi)
        phi_b += dphi
    phi = np.array(phi)
    print("  phi_b = {:7.5f}".format(phi_b))
    
    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol
    # CG     gtol
    minimum = opt.minimize(
        f_super_per,
        phi,
        method="Powell",
        # callback=prt_iter,
        # bounds=bounds,
        # bounds = opt.Bounds(np.full(n_phi, -dhi_max), np.full(n_phi, dhi_max))
        # options={"xtol": x_tol, "maxiter": max_iter},
        # options={"gtol": f_tol, "maxiter": max_iter},
    )


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
file_name = os.path.join(home_dir, "max_4u.lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi()))
lat_prop.prt_M()
if not stable:
    assert False
Twiss = lat_prop.get_Twiss(len(lat_prop._lattice)-1)
lat_prop.prt_Twiss_param(Twiss)

# Compute radiation properties.
stable = lat_prop.compute_radiation()
lat_prop.prt_rad()
lat_prop.prt_M_rad()

lat_prop.get_types()
lat_prop.prt_Twiss()

if False:
    lat_prop.plt_Twiss("unit_cell.png", False)

if not False:
    opt_var_bend_radius(lat_prop)
