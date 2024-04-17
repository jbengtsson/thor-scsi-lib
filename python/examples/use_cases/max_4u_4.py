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


class bend_prm_class:
    # Private
    def __init__(self, lat_prop, dip_0, dip_list, b_2_prm):
        self._dip_0    = dip_0
        self._phi_0    = lat_prop.get_phi_elem(dip_0, 0)
        self._dip_list = dip_list
        self._phi      = self.compute_phi_bend()
        self._b_2_prm  = b_2_prm

        print("\n{:s}:".format(dip_0))
        print("  phi   = {:8.5f}".format(self._phi))
        print("  phi_0 = {:8.5f}".format(self._phi_0))

    # Public.
    def init(self):
        prm = []
        for k in range(len(self._dip_list)):
            dphi = lat_prop.get_phi_elem(self._dip_list[k], 0)
            prm.append(dphi)
        prm = np.array(prm)
        return prm

    def compute_phi_bend(self):
        phi = 0e0
        for k in range(len(self._dip_list)):
            phi += lat_prop.get_phi_elem(self._dip_list[k], 0)
        phi += lat_prop.get_phi_elem(self._dip_0, 0)
        return phi

    def set_phi(self, prm):
        phi = 0e0
        for k in range(len(self._dip_list)):
            lat_prop.set_phi_fam(self._dip_list[k], prm[k])
            phi += prm[k]
        lat_prop.set_phi_fam(self._dip_0, self._phi-phi)


def prt_lat(lat_prop, file_name, bend):
    outf = open(file_name, 'w')

    def prt_bend(lat_prop, name):
        L = lat_prop.get_L_elem(name, 0)
        phi = lat_prop.get_phi_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:s}: Bending, L = {:7.5f}, T = {:8.5f}, K = {:8.5f}"
               ", T1 = 0.0, T2 = 0.0,\n     N = n_bend;")
              .format(name, L, phi, b_2), file=outf)

    prt_bend(lat_prop, bend._dip_0)
    for k in range(len(bend._dip_list)):
        prt_bend(lat_prop, bend._dip_list[k])


def opt_bend(lat_prop, prm_list, weight):
    """Use Case: optimise super period.
    """

    loc       = lat_prop._lattice.find("sf_h", 0).index
    eta       = np.nan
    chi_2_min = 1e30
    n_iter    = 0
    file_name = "opt_sp.txt"

    def prt_prms(prm):
        n_prt = 5
        print("\n  prm =")
        for k in range(len(prm)):
            print("  " if k == 0 else "", "{:12.5e}".format(prm[k]),
                  "\n  " if (k != 0) and (k % n_prt == n_prt-1) else "", end="")
        if k % n_prt != n_prt-1:
            print()

    def prt_iter(prm, chi_2):
        nonlocal loc, n_iter

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:5.3f}".format(1e12*lat_prop._eps[ind.X]))
        print("  eta'_x [m]     = {:10.3e}".format(eta[ind.px]))
        print("  phi_lat        = {:7.5f}".format(lat_prop.compute_phi_lat()))
        prt_prms(prm)

    def compute_chi_2(prm):
        nonlocal loc, eta, n_iter

        prt = not False

        n_iter += 1
        prm_list[0][1].set_phi(prm)
        stable = lat_prop.comp_per_sol()
        if stable:
            stable = lat_prop.compute_radiation()
        else:
            print("\ncompute_chi_2 - unstable")
        dchi_2 = weight[0]*lat_prop._eps[ind.X]**2 if stable else 1e30
        if prt:
            print("\n  dchi2(eps_x)  = {:10.3e}".format(dchi_2))
        chi_2 = dchi_2
        eta = lat_prop.get_Twiss(loc)[0]
        dchi_2 = weight[1]*eta[ind.px]**2
        if prt:
            print("  dchi2(eta'_x) = {:10.3e}".format(dchi_2))
        chi_2 += dchi_2
        return chi_2

    def f_super_per(prm):
        nonlocal chi_2_min

        chi_2 = compute_chi_2(prm)
        if chi_2 < chi_2_min:
            prt_iter(prm, chi_2)
            prt_lat(lat_prop, "opt_bend.txt", prm_list[0][1])
            chi_2_min = min(chi_2, chi_2_min)
        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    prm = prm_list[0][1].init()

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
lat_name = "max_4u_match"
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

if False:
    lat_prop.prt_lat_param()
    lat_prop.prt_Twiss("max_4u_uc.txt")

# Weights.
weight = np.array([
    1e0,  # eps_x.
    1e-14 # eta'_x.
])

dip_list = ["b1", "b2", "b3", "b4", "b5"]
b1 = bend_prm_class(lat_prop, "b0", dip_list, True)

dip_list  = [
    "ds6", "ds5", "ds4", "ds3", "ds2", "ds1", "dm1", "dm2", "dm3", "dm4", "dm5"
]
b2 = bend_prm_class(lat_prop, "ds0", dip_list, True)

prm_list = [
    ("bend", b2)
]

if not False:
    opt_bend(lat_prop, prm_list, weight)

if False:
    dip_list = [bend]
    dip_list.extend(bend_list)
    prt_lat(lat_prop, "opt_bend.txt", dip_list)
