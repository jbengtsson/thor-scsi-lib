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
    linear_optics as lo
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()

quad = 2

L_min        = 0.10
L_max        = 0.60
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0

class bend_prm_class:
    # Private
    def __init__(self, lat_prop, dip_0, dip_list, b_2_prm):
        self._dip_0        = dip_0
        self._phi_0        = lat_prop.get_phi_elem(self._dip_0, 0)
        self._dip_list     = dip_list
        self._phi          = self.compute_phi_bend()
        self._b_2_prm      = b_2_prm

        print("\n{:s}:".format(dip_0))
        print("  phi   = {:8.5f}".format(self._phi))
        print("  phi_0 = {:8.5f}".format(self._phi_0))

    # Public.
    def init(self):
        prm = []
        bounds = []
        for k in range(len(self._dip_list)):
            dphi = lat_prop.get_phi_elem(self._dip_list[k], 0)
            prm.append(dphi)
            bounds.append((-phi_max, phi_max))
            if self._b_2_prm:
                b_2 = lat_prop.get_b_n_elem(self._dip_list[k], 0, quad)
                prm.append(b_2)
                bounds.append((-b_2_max, b_2_max))
        if self._b_2_prm:
            b_2 = lat_prop.get_b_n_elem(self._dip_0, 0, quad)
            prm.append(b_2)
            bounds.append((-b_2_max, b_2_max))
        return prm, bounds

    def compute_phi_bend(self):
        phi = 0e0
        for k in range(len(self._dip_list)):
            phi += lat_prop.get_phi_elem(self._dip_list[k], 0)
        phi += lat_prop.get_phi_elem(self._dip_0, 0)
        return phi

    def set_bend(self, prm, prm_ind):
        phi = 0e0
        for k in range(len(self._dip_list)):
            lat_prop.set_phi_fam(self._dip_list[k], prm[prm_ind])
            phi += prm[prm_ind]
            prm_ind += 1
            if self._b_2_prm:
                lat_prop.set_b_n_fam(self._dip_list[k], quad, prm[prm_ind])
                prm_ind += 1
        self._phi_0 = self._phi - phi
        if self._b_2_prm:
            lat_prop.set_b_n_fam(self._dip_0, quad, prm[prm_ind])
            prm_ind += 1
        lat_prop.set_phi_fam(self._dip_0, self._phi_0)
        return prm_ind


def prt_lat(lat_prop, file_name, bend):
    outf = open(file_name, 'w')

    def prt_drift(name):
        L = lat_prop.get_L_elem(name, 0)
        print(("{:4s}: Drift, L = {:7.5f};").format(name, L), file=outf)

    def prt_bend(lat_prop, name):
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

    prt_bend(lat_prop, bend._dip_0)
    for k in range(len(bend._dip_list)):
        prt_bend(lat_prop, bend._dip_list[k])


def opt_bend(lat_prop, prm_list, weight):
    """Use Case: optimise super period.
    """

    loc       = lat_prop._lattice.find("sf_h", 0).index
    eta_x     = 0.06255
    eta       = np.nan
    chi_2_min = 1e30
    n_iter    = 0
    file_name = "opt_sp.txt"

    # Dictionary of parameter types and corresponding get functions.
    get_prm_func_dict = {
        "L":   lat_prop.get_L_elem,
        "L_b": lat_prop.get_L_elem,
        "phi": lat_prop.get_phi_elem,
        "b_2": lat_prop.get_b_2_elem
    }

    # Dictionary of parameter types and corresponding set functions.
    set_prm_func_dict = {
        "L":   lat_prop.set_L_fam,
        "L_b": lat_prop.set_L_bend_fam,
        "phi": lat_prop.set_phi_fam,
        "b_2": lat_prop.set_b_2_fam
    }

    # Dictionary of parameter bounds.
    prm_bounds_dict = {
        "L":   ( L_min,   L_max),
        "L_b": ( L_min,   L_max),
        "phi": (-phi_max, phi_max),
        "b_2": (-b_2_max, b_2_max)
    }

    def get_prm():
        prm = []
        bounds = []
        for k in range(len(prm_list)):
            if prm_list[k][0] == "bend":
                p, b = prm_list[k][1].init()
                prm.extend(p)
                bounds.extend(b)
            else:
                prm.append(get_prm_func_dict[prm_list[k][1]](prm_list[k][0], 0))
                bounds.append((-self._b_2_max, self._b_2_max))
        return np.array(prm), bounds

    def set_prm(prm):
        prm_ind = 0
        for k in range(len(prm_list)):
            if prm_list[k][0] == "bend":
                prm_ind = prm_list[k][1].set_bend(prm, prm_ind)
            else:
                set_prm_func_dict[prm_list[k][1]](prm_list[k][0], prm[k])
                prm_ind += 1
        
    def prt_prm(prm):
        n_prt = 5
        print("\n  prm =")
        for k in range(len(prm)):
            print("  " if k == 0 else "", "{:12.5e}".format(prm[k]),
                  "\n  " if (k != 0) and (k % n_prt == n_prt-1) else "", end="")
        if k % n_prt != n_prt-1:
            print()

    def prt_iter(prm, chi_2, eta, xi):
        nonlocal n_iter

        phi_0 = lat_prop.get_phi_elem(prm_list[0][1]._dip_0, 0)

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:5.3f}".format(1e12*lat_prop._eps[ind.X]))
        print("  eta_x [m]      =  {:10.3e} ({:9.3e})".
              format(eta[ind.x], eta_x))
        print("  eta'_x [m]     =  {:10.3e}".format(eta[ind.px]))
        print("  {:s}_phi        =  {:6.3f}".
              format(prm_list[0][1]._dip_0, phi_0))
        print("  xi             =  [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        prt_prm(prm)

    def compute_chi_2():
        nonlocal loc, eta_x

        prt = not False

        try:
            if not lat_prop.comp_per_sol():
                print("\ncompute_chi_2 - comp_per_sol: unstable")
                raise ValueError
            if not lat_prop.compute_radiation():
                print("\ncompute_chi_2 - compute_radiation: unstable")
                raise ValueError
        except ValueError:
            chi_2 = 1e30
            eta = np.nan
            xi = np.nan
        else:
            dchi_2 = weight[0]*lat_prop._eps[ind.X]**2
            chi_2 = dchi_2
            if prt:
                print("\n  dchi2(eps_x)  = {:10.3e}".format(dchi_2))

            eta = lat_prop.get_Twiss(loc)[0]
            dchi_2 = weight[1]*(eta[ind.x]-eta_x)**2
            chi_2 += dchi_2
            if prt:
                print("  dchi2(eta_x)  = {:10.3e}".format(dchi_2))

            dchi_2 = weight[2]*eta[ind.px]**2
            chi_2 += dchi_2
            if prt:
                print("  dchi2(eta'_x) = {:10.3e}".format(dchi_2))

            M = lo.compute_map(
                lat_prop._lattice, lat_prop._model_state,
                desc=lat_prop._desc, tpsa_order=2)
            stable, nu, xi = \
                lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
            if not stable:
                xi[ind.X] = 1e30
                xi[ind.Y] = 1e30
            dchi_2 = weight[3]*(xi[ind.X]**2+xi[ind.Y]**2)
            chi_2 += dchi_2
            if prt:
                print("  dchi2(xi)     = {:10.3e}".format(dchi_2))

        return chi_2, eta, xi

    def f_super_per(prm):
        nonlocal chi_2_min, n_iter

        n_iter += 1
        set_prm(prm)
        chi_2, eta, xi = compute_chi_2()
        if chi_2 < chi_2_min:
            prt_iter(prm, chi_2, eta, xi)
            prt_lat(lat_prop, "opt_bend.txt", prm_list[0][1])
            chi_2_min = min(chi_2, chi_2_min)
        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    prm, bounds = get_prm()

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
        bounds = bounds,
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
lat_name = "max_4u_sp"
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

# Compute radiation properties.
stable = lat_prop.compute_radiation()
if not stable:
    print("\ncompute_radiation - unstable")
    assert False
lat_prop.prt_rad()
lat_prop.prt_M_rad()

if False:
    lat_prop.prt_lat_param()
    lat_prop.prt_Twiss("max_4u_uc.txt")

# Weights.
weight = np.array([
    1e14, # eps_x.
    1e-1, # eta_x.
    1e0,  # eta'_x.
    1e-7  # xi.
])

dip_list = ["b1", "b2", "b3", "b4", "b5"]
b1 = bend_prm_class(lat_prop, "b0", dip_list, False)

dip_list  = [
    "ds6", "ds5", "ds4", "ds3", "ds2", "ds1", "dm1", "dm2", "dm3", "dm4", "dm5"
]
b2 = bend_prm_class(lat_prop, "ds0", dip_list, True)

prm_list = [
    # ("qf1e", "b_2"),
    # ("qd",   "b_2"),
    # ("qf2",  "b_2"),
    # ("ds0",  "L_b"),
    ("bend", b2)
]

if not False:
    opt_bend(lat_prop, prm_list, weight)

if False:
    dip_list = [bend]
    dip_list.extend(bend_list)
    prt_lat(lat_prop, "opt_bend.txt", dip_list)
