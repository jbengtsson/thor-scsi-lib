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
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()


def compute_phi_bend(bend_list):
    phi_b = 0e0
    for k in range(len(bend_list)):
        phi_b += lat_prop.get_phi_elem(bend_list[k], 0)
    return phi_b


def prt_lat(lat_prop, file_name, bend_list):
    outf = open(file_name, 'w')

    def prt_bend(lat_prop, name):
        L = lat_prop.get_L_elem(name, 0)
        phi = lat_prop.get_phi_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:s}: Bending, L = {:7.5f}, T = {:8.5f}, K = {:8.5f}"
               ", T1 = 0.0, T2 = 0.0,\n     N = n_bend;")
              .format(name, L, phi, b_2), file=outf)

    for k in range(len(bend_list)):
        prt_bend(lat_prop, bend_list[k])
    print(file=outf)


def opt_bend(lat_prop, bend, bend_list):
    """Use Case: optimise matching cell.
    """

    chi_2_min = 1e30
    n_iter    = 0
    file_name = "opt_bend.txt"
    dip_list  = [bend]
    dip_list.extend(bend_list)
    phi_bend  = compute_phi_bend(dip_list)

    def set_phi(prm):
        phi = 0e0
        for k in range(len(bend_list)):
            lat_prop.set_phi_fam(bend_list[k], prm[k])
            phi += prm[k]
        lat_prop.set_phi_fam(bend, phi_bend-phi)

    def prt_prms(prm):
        n_prt = 5
        print("\n  prm =")
        for k in range(len(prm)):
            print(" {:12.5e}".format(prm[k]), end="")
            if k % n_prt == 4:
                print()
        if k % n_prt != 4:
            print()

    def prt_iter(prm, chi_2):
        nonlocal n_iter

        phi_0 = lat_prop.get_phi_elem(bend, 0)
        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:5.3f}".format(1e12*lat_prop._eps[ind.X]))
        print("  phi_bend       = {:7.5f}".format(lat_prop.compute_phi_lat()))
        prt_prms(prm)

    def compute_chi_2(prm):
        nonlocal n_iter

        n_iter += 1
        set_phi(prm)
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

    def f_super_per(prm):
        nonlocal chi_2_min

        chi_2 = compute_chi_2(prm)
        if chi_2 < chi_2_min:
            prt_iter(prm, chi_2)
            prt_lat(lat_prop, "opt_bend.txt", dip_list)
            chi_2_min = min(chi_2, chi_2_min)
        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    phi_b = lat_prop.get_phi_elem(bend, 0)
    prm = []
    for k in range(len(bend_list)):
        dphi = lat_prop.get_phi_elem(bend_list[k], 0)
        prm.append(dphi)
        phi_b += dphi
    prm = np.array(prm)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol, xtol
    # CG     gtol
    minimum = opt.minimize(
        f_super_per,
        prm,
        method="Powell",
        # callback=prt_iter,
        # bounds = opt.Bounds(np.full(n_phi, -dhi_max), np.full(n_phi, dhi_max))
        options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
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
# lat_name = "max_iv"
# lat_name = "max_4u_uc"
lat_name = "max_4u_match"
# lat_name = "max_4u_sup_per"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.get_types()

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))
lat_prop.prt_M()
if not stable:
    assert False

# Compute radiation properties.
stable = lat_prop.compute_radiation()
lat_prop.prt_rad()
lat_prop.prt_M_rad()

if not False:
    lat_prop.prt_lat_param()
    lat_prop.prt_Twiss("max_4u_uc.txt")

bend_list  = [
    "ds6", "ds5", "ds4", "ds3", "ds2", "ds1", "dm1", "dm2", "dm3", "dm4", "dm5"
]
bend = "ds0"

if False:
    lat_prop.plt_Twiss(lat_name+".png", not False)

if not False:
    opt_bend(lat_prop, bend, bend_list)

if False:
    dip_list = [bend]
    dip_list.extend(bend_list)
    prt_lat(lat_prop, "opt_bend.txt", dip_list)
