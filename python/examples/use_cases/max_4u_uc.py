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


def opt_sp(lat_prop, prm_list, weight):
    """Use Case: optimise super period.
    """

    nu_uc     = [2.0/6.0+0.005, 1.0/12.0]

    chi_2_min = 1e30
    eta       = np.nan
    n_iter    = 0
    file_name = "opt_sp.txt"

    def prt_iter(prm, chi_2, nu, xi):
        nonlocal n_iter

        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        b_list = [
            "b1_0", "b1_1", "b1_2", "b1_3", "b1_4", "b1_5"]
        rb = "qf1"

        phi = lat_prop.compute_phi_lat()
        phi_b = compute_phi_bend(lat_prop, b_list)
        phi_rb = lat_prop.get_phi_elem(rb, 0)

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:5.3f}".format(1e12*lat_prop._eps[ind.X]))
        print("  nu_uc          =  [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu[ind.X], nu[ind.Y], nu_uc[ind.X], nu_uc[ind.Y]))
        print("  xi             =  [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        print("\n  phi_sp         =  {:8.5f}".format(phi))
        print("  C [m]          =  {:8.5f}".format(lat_prop.compute_circ()))
        print("\n  phi_b         =  {:8.5f}".format(phi_b))
        print("  phi_rb        =  {:8.5f}".format(phi_rb))
        prm_list.prt_prm(prm)

    def compute_chi_2(nu, xi):
        prt = not False

        dchi_2 = weight[0]*lat_prop._eps[ind.X]**2
        chi_2 = dchi_2
        if prt:
            print("\n  dchi2(eps_x)    = {:10.3e}".format(dchi_2))

        dchi_2 = \
            weight[1]*(
                (nu[ind.X]-nu_uc[ind.X])**2
                +(nu[ind.Y]-nu_uc[ind.Y])**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc)    = {:10.3e}".format(dchi_2))

        dchi_2 = weight[2]*(xi[ind.X]**2+xi[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(xi)       = {:10.3e}".format(dchi_2))

        return chi_2

    def f_sp(prm):
        nonlocal chi_2_min, n_iter

        n_iter += 1
        prm_list.set_prm(prm)

        # Compute the beam dynamics properties.
        try:
            if not lat_prop.comp_per_sol():
                print("\nf_sp - comp_per_sol: unstable")
                raise ValueError

            if not lat_prop.compute_radiation():
                print("\nf_sp - compute_radiation: unstable")
                raise ValueError

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
            chi_2 = compute_chi_2(nu, xi)
            if chi_2 < chi_2_min:
                prt_iter(prm, chi_2, nu, xi)
                pc.prt_lat(lat_prop, "opt_sp.txt", prm_list)
                chi_2_min = min(chi_2, chi_2_min)
            return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

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
lat_name = "max_4u_sp_jb_2"
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
    1e14,  # eps_x.
    1e2,   # eta_uc_x.
    1e-7   # xi.
])

bend_list = [
    "b1_1", "b1_2", "b1_3", "b1_4", "b1_5",
    "qf1"]

# opt_phi = pc.opt_phi_class(lat_prop, "b1_0", bend_list, phi_max)

prm_list = [
    ("qf1",      "b_2"),

    ("b1_0",     "b_2"),
    ("b1_1",     "b_2"),
    ("b1_2",     "b_2"),
    ("b1_3",     "b_2"),
    ("b1_4",     "b_2"),
    ("b1_5",     "b_2"),

    # ("opt_phi",  opt_phi)
]

prm_list = pc.prm_class(lat_prop, prm_list, b_2_max)

opt_sp(lat_prop, prm_list, weight)

if False:
    dip_list = [bend]
    dip_list.extend(bend_list)
    prt_lat(lat_prop, "opt_sp.txt", dip_list)
