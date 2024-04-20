"""Use Case:
     Optimise varying bend radius dipole.
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


def compute_phi_b(bend_list):
    phi_b = 0e0
    for k in range(len(bend_list)):
        phi_b += lat_prop.get_phi_elem(bend_list[k], 0)
    return phi_b


def prt_bend_lat(lat_prop, file_name, bend_list, rbend_name):
    outf = open(file_name, 'w')

    def prt_b(lat_prop, name):
        L = lat_prop.get_L_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:s}: Bending, L = {:7.5f}, T = {:s}_scl*phi, K = {:8.5f}"
               ", T1 = 0.0, T2 = 0.0,\n    N = n_bend;")
              .format(name, L, name, b_2), file=outf)

    def prt_rb(lat_prop, name):
        L = lat_prop.get_L_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:s}: Bending, L = {:7.5f}, T = phi_rb, K = {:8.5f}"
               ", N = n_bend;").format(name, L, b_2), file=outf)

    print("phi_b  = {:7.5f};\n".
          format(compute_phi_b(bend_list)), file=outf)
    for k in range(len(bend_list)):
        phi = lat_prop.get_phi_elem(bend_list[k], 0)
        print("{:s}_scl = {:7.5f}/phi_b;".format(bend_list[k], phi), file=outf)
    print("\nphi_rb = 0.0;", file=outf)
    print("phi    = phi_b - phi_rb;\n", file=outf)
    for k in range(len(bend_list)):
        prt_b(lat_prop, bend_list[k])
    print(file=outf)
    prt_rb(lat_prop, rbend_name)


def opt_var_bend_radius(lat_prop, bend_list, rbend_name, phi_uc, rbend):
    """Use Case: optimise unit cell.
    """

    chi_2_min = 1e30
    n_iter    = 0
    file_name = "opt_var_bend_radius.txt"

    def set_phi_var_bend_rad_rb(phi):
        dphi = 0e0
        for k in range(len(bend_list)):
            lat_prop.set_phi_fam(bend_list[k], phi[k])
            dphi += phi[k]
        if rbend:
            lat_prop.set_phi_fam(rbend_name, phi[-1])
            dphi += phi[-1]
        lat_prop.set_phi_fam("b0", phi_uc/2e0-dphi)

    def prt_iter(phi, chi_2):
        nonlocal n_iter

        phi_0 = lat_prop.get_phi_elem("b0", 0)
        if rbend:
            phi_rb = lat_prop.get_phi_elem(rbend_name, 0)
        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:14.10f}".format(1e12*lat_prop._eps[ind.X]))
        print("  nu             = [{:7.5f}, {:7.5f}]".
              format(lat_prop._nu[ind.X], lat_prop._nu[ind.Y]))
        print("  phi_uc         = {:7.5f}".format(lat_prop.compute_phi()))
        print("  phi            = {:12.5e}".format(phi_0)+vec2txt(phi))

    def compute_chi_2(phi):
        nonlocal n_iter

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
        nonlocal chi_2_min

        chi_2 = compute_chi_2(phi)
        if chi_2 < chi_2_min:
            prt_iter(phi, chi_2)
            chi_2_min = min(chi_2, chi_2_min)
        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    print("\nopt_var_bend_radius:")

    phi_b = lat_prop.get_phi_elem("b1_0", 0)
    phi = []
    for k in range(len(bend_list)):
        dphi = lat_prop.get_phi_elem(bend_list[k], 0)
        phi.append(dphi)
        phi_b += dphi
    if rbend:
        dphi = lat_prop.get_phi_elem(rbend_name, 0)
        phi.append(dphi)
    print("\nphi_b  = {:7.5f}".format(phi_b))
    print("phi_uc = {:7.5f}".format(phi_uc))
    phi = np.array(phi)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol, xtol
    # CG     gtol
    minimum = opt.minimize(
        f_super_per,
        phi,
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
# lat_name = "max_4u_match"
lat_name = "max_4u_sp_1"
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
    print("\ncomp_per_sol - unstable")
    assert False
lat_prop.prt_lat_param()

# Compute radiation properties.
stable = lat_prop.compute_radiation()
lat_prop.prt_rad()
lat_prop.prt_M_rad()

lat_prop.prt_Twiss("max_4u_uc.txt")

bend_list  = ["b1_1", "b1_2", "b1_3", "b1_4", "b1_5"]
rbend_name = "qf1"
phi_uc     = 3.0

if not False:
    lat_prop.plt_Twiss(lat_name+".png", not False)

if False:
    opt_var_bend_radius(lat_prop, bend_list, rbend_name, phi_uc, not False)

if False:
    dip_list = ["b1_0"]
    dip_list.extend(bend_list)
    prt_bend_lat \
        (lat_prop, "opt_var_bend_radius.txt", dip_list, rbend_name)
