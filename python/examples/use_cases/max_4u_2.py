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


def prt_bend_lat(lat_prop, file_name, bend_list):
    chi_2_min = 1e30
    n_iter    = 0
    outf      = open(file_name, 'w')

    def prt_bend(lat_prop, name):
        L = lat_prop.get_L_elem(name, 0)
        phi = lat_prop.get_phi_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:2s}: Bending, L = {:7.5f}, T = {:7.5f}, K = {:8.5f}, T1 = 0.0"
               ", T2 = 0.0,\n    N = nbend;").format(name, L, phi, b_2),
              file=outf)

    prt_bend(lat_prop, "b0")
    for k in range(len(bend_list)):
        prt_bend(lat_prop, bend_list[k])


def opt_var_bend_radius(lat_prop, bend_list):
    """Use Case: optimise unit cell.
    """

    file_name = "opt_var_bend_radius.txt"

    def prt_bend_list(lat_prop, bend_list):
        phi = lat_prop.get_phi_elem("b0", 0)
        print("[{:7.5f}, ".format(phi), end="")
        n = len(bend_list)
        for k in range(n):
            phi = lat_prop.get_phi_elem(bend_list[k], 0)
            if k < n-1:
                print("{:7.5f}, ".format(phi), end="")
            else:
                print("{:7.5f}]".format(phi))

    def set_phi_var_bend_rad_rb(phi):
        phi_tot = 0e0
        for k in range(len(bend_list)):
            lat_prop.set_phi_fam(bend_list[k], phi[k], False)
            phi_tot += phi[k]
        lat_prop.set_phi_fam("b0", phi_b-phi_tot, False)

    def prt_iter(phi, chi_2):
        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:14.10f}".format(1e12*lat_prop._eps[ind.X]))
        print("  nu             = [{:7.5f}, {:7.5f}]".
              format(lat_prop._nu[ind.X], lat_prop._nu[ind.Y]))
        print("  phi_tot        = {:7.5f}".format(lat_prop.compute_phi()))
        print("  phi            = ", end="")
        prt_bend_list(lat_prop, bend_list)

    def compute_chi_2(phi):
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
    phi = []
    for k in range(len(bend_list)):
        dphi = lat_prop.get_phi_elem(bend_list[k], 0)
        phi.append(dphi)
        phi_b += dphi
    phi = np.array(phi)
    print("\nphi_b = {:7.5f}".format(phi_b))

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

    prt_bend_lat(lat_prop, file_name, bend_list)


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
file_name = os.path.join(home_dir, "max_4u_1.lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi()))
lat_prop.prt_M()
if not stable:
    assert False
lat_prop.prt_lat_param()

# Compute radiation properties.
stable = lat_prop.compute_radiation()
lat_prop.prt_rad()
lat_prop.prt_M_rad()

lat_prop.get_types()
lat_prop.prt_Twiss()

if False:
    lat_prop.plt_Twiss("unit_cell.png", False)

if not False:
    bend_list  = ["b1", "b2", "b3", "b4", "b5"]

    opt_var_bend_radius(lat_prop, bend_list)
