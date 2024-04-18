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
L_max        = 0.50
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


def prt_lat(lat_prop, file_name, prm_list):
    outf = open(file_name, 'w')

    def prt_drift(name):
        L = lat_prop.get_L_elem(name, 0)
        print(("{:5s}: Drift, L = {:7.5f};").format(name, L), file=outf)

    def prt_dip(name):
        L = lat_prop.get_L_elem(name, 0)
        phi = lat_prop.get_phi_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:5s}: Bending, L = {:7.5f}, T = {:8.5f}, K = {:8.5f}"
               ", T1 = 0.0, T2 = 0.0,\n       N = n_bend;")
              .format(name, L, phi, b_2), file=outf)

    def prt_bend(bend):
        prt_dip(bend._dip_0)
        for k in range(len(bend._dip_list)):
            prt_dip(bend._dip_list[k])

    def prt_quad(name):
        L = lat_prop.get_L_elem(name, 0)
        b_2 = lat_prop.get_b_n_elem(name, 0, 2)
        print(("{:5s}: Quadrupole, L = {:7.5f}, K = {:8.5f}, N = n_quad;")
              .format(name, L, b_2), file=outf)

    # Dictionary of parameter types and corresponding print functions.
    get_prm_func_dict = {
        "L":   prt_drift,
        "L_b": prt_dip,
        "phi": prt_dip,
        "b_2": prt_dip
    }

    for k in range(len(prm_list)):
        if prm_list[k][0] == "bend":
             prt_bend(prm_list[k][1])
        else:
            get_prm_func_dict[prm_list[k][1]](prm_list[k][0])


def opt_bend(lat_prop, prm_list, weight):
    """Use Case: optimise super period.
    """

    loc       = lat_prop._lattice.find("sf_h", 2).index
    # Dispersion at the unit cell end.
    eta_uc_x  = 0.06255
    # Beta functions at the unit cell end.
    beta_uc   = [3.65614, 3.68868]

    chi_2_min = 1e30
    eta       = np.nan
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
                bounds.append((-b_2_max, b_2_max))
        return np.array(prm), bounds

    def set_prm(prm):
        prm_ind = 0
        for k in range(len(prm_list)):
            if prm_list[k][0] == "bend":
                prm_ind = prm_list[k][1].set_bend(prm, prm_ind)
            else:
                set_prm_func_dict[prm_list[k][1]](prm_list[k][0], prm[prm_ind])
                prm_ind += 1
        
    def prt_prm(prm):
        n_prt = 5
        print("\n  prm =")
        for k in range(len(prm)):
            print("  " if k == 0 else "", "{:12.5e}".format(prm[k]),
                  "\n  " if (k != 0) and (k % n_prt == n_prt-1) else "", end="")
        if k % n_prt != n_prt-1:
            print()

    def prt_iter(prm, chi_2, eta, beta, xi):
        nonlocal n_iter

        phi = lat_prop.compute_phi_lat()
        # Twiss functions at end of super period.
        _, _, beta_id, _ = lat_prop.get_Twiss(-1)

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("  eps_x [pm.rad] = {:5.3f}".format(1e12*lat_prop._eps[ind.X]))
        print("  eta_x [m]      =  {:10.3e} ({:9.3e})".
              format(eta[ind.x], eta_uc_x))
        print("  beta           =  [{:5.3f}, {:5.3f}] ([{:5.3f}, {:5.3f}])".
              format(beta[ind.X], beta[ind.Y],
                     beta_uc[ind.X], beta_uc[ind.Y]))
        print("  beta_id        =  [{:5.3f}, {:5.3f}]".
              format(beta_id[ind.X], beta_id[ind.Y]))
        print("  xi             =  [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        print("  phi_sp         =  {:8.5f}".format(phi))
        prt_prm(prm)

    def compute_chi_2():
        nonlocal loc, eta_uc_x

        prt = not False

        try:
            if not lat_prop.comp_per_sol():
                print("\ncompute_chi_2 - comp_per_sol: unstable")
                raise ValueError

            if not lat_prop.compute_radiation():
                print("\ncompute_chi_2 - compute_radiation: unstable")
                raise ValueError

            M = lo.compute_map(
                lat_prop._lattice, lat_prop._model_state,
                desc=lat_prop._desc, tpsa_order=2)
            stable, nu, xi = \
                lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
            if not stable:
                print("\ncompute_chi_2 - compute_nu_xi: unstable")
                raise ValueError
        except ValueError:
            chi_2 = 1e30
            eta = np.nan
            beta = np.nan
            xi = np.nan
        else:
            dchi_2 = weight[0]*lat_prop._eps[ind.X]**2
            chi_2 = dchi_2
            if prt:
                print("\n  dchi2(eps_x)     = {:10.3e}".format(dchi_2))

            eta, alpha, beta, nu = lat_prop.get_Twiss(loc)

            dchi_2 = weight[1]*(eta[ind.x]-eta_uc_x)**2
            chi_2 += dchi_2
            if prt:
                print("  dchi2(eta_uc_x)  = {:10.3e}".format(dchi_2))

            dchi_2 = \
                weight[2]*(
                    (beta[ind.X]-beta_uc[ind.X])**2
                    +(beta[ind.Y]-beta_uc[ind.Y])**2)
            chi_2 += dchi_2
            if prt:
                print("  dchi2(beta)      = {:10.3e}".format(dchi_2))

            dchi_2 = weight[3]*(xi[ind.X]**2+xi[ind.Y]**2)
            chi_2 += dchi_2
            if prt:
                print("  dchi2(xi)        = {:10.3e}".format(dchi_2))

        return chi_2, eta, beta, xi

    def f_super_per(prm):
        nonlocal chi_2_min, n_iter

        n_iter += 1
        set_prm(prm)
        chi_2, eta, beta, xi = compute_chi_2()
        if chi_2 < chi_2_min:
            prt_iter(prm, chi_2, eta, beta, xi)
            prt_lat(lat_prop, "opt_bend.txt", prm_list)
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
lat_name = "max_4u_sp"
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
finally:
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
    1e-2 , # beta_uc.
    1e-6   # xi.
])

dip_list = ["b1_1", "b1_2", "b1_3", "b1_4", "b1_5"]
b1 = bend_prm_class(lat_prop, "b1_0", dip_list, False)

dip_list  = [
    "b2u_6", "b2u_5", "b2u_4", "b2u_3", "b2u_2", "b2u_1", "b2d_1", "b2d_2",
    "b2d_3", "b2d_4", "b2d_5"
]
b2 = bend_prm_class(lat_prop, "b2_0", dip_list, True)

prm_list = [
    ("qf1_e", "b_2"),
    ("qd",    "b_2"),
    ("qf2",   "b_2"),
    ("b1_0",  "L_b"),
    ("b2_0",  "L_b"),
    ("bend",  b1),
    ("bend",  b2)
]

if not False:
    opt_bend(lat_prop, prm_list, weight)

if False:
    dip_list = [bend]
    dip_list.extend(bend_list)
    prt_lat(lat_prop, "opt_bend.txt", dip_list)
