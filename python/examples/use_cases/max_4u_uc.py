"""Use Case:
     Optimising a periodic unit cell.
"""


import os
import sys
from dataclasses import dataclass
from typing import ClassVar
import enum
from typing import Tuple

import numpy as np
from scipy import optimize as opt
from scipy import linalg as la

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, knobs as kb, prm_class as pc, nonlin_dyn as nld_class
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()

prm_range = {
    "phi":       [  0.0,  2.0],
    "phi_rbend": [ -0.5,  0.0],
    "b_2":       [-10.0, 10.0],
    "b_2_bend":  [ -1.5,  1.5]
}

design_val = {
    "eps_x_des" : 50e-12,
    "nu_uc_des" : np.array([3.0/7.0, 1.0/7.0]),
    "eta_x_des" : 1e-4
}

L_min        = 0.10
L_max        = 0.50
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0


class opt_uc_class:
    # Private.

    def __init__(self, prm_class: pc.prm_class) -> None:
        self._first        = True

        self._lat_prop     = prm_class.lattice[0]
        self._nld          = prm_class.lattice[1]
        self._uc_centre      = prm_class.s_loc
        self._des_val_list = prm_class.design_vals
        self._weights      = prm_class.weights
        self._bend_list    = prm_class.dipoles[0]
        self._rbend_list   = prm_class.dipoles[1]
        self._prm_list     = prm_class.params[0]
        self._dprm_list    = prm_class.params[1]

        self._phi_bend     = np.zeros(len(self._bend_list))
        self._phi_rbend    = np.zeros(len(self._bend_list))

        self._b_2          = np.nan

        self._nu_uc        = np.nan
        self._eta_centre   = np.nan
        self._xi           = np.nan

        self._constr       = {}

        self._chi_2_min    = 1e30
        self._n_iter       = -1
        self._file_name    = "opt_uc.txt"

    # Public.

    def compute_lat_prop(self) -> bool:
        # Compute Twiss parameters along the lattice.
        stable = self._lat_prop.comp_per_sol()
        if not stable:
            print("\ncomp_per_sol: unstable")
            return False

        # Compute radiation properties.
        stable, stable_rad = self._lat_prop.compute_radiation()
        if not stable or not stable_rad:
            print("\ncompute_radiation: unstable")
            return False

        # Compute linear chromaticity.
        self._nld.compute_map()
        stable, _, self._xi = lo.compute_nu_xi(
            self._lat_prop.desc, self._lat_prop._no, self._nld._M)
        if not stable:
            print("\ncompute_nu_xi: unstable")
            return False
        return True

    def compute_prop(self) -> bool:
        self._phi_sp = self._lat_prop.compute_phi_lat()
        self._dphi = self._phi_sp - self._phi_sp_0

        for k, bend in enumerate(self._bend_list):
            self._phi_bend[k] = bend.compute_bend_phi()

        for k, bend in enumerate(self._rbend_list):
            self._phi_rbend[k] = self._lat_prop.get_phi_elem(bend, 0)

        lat_stable = self.compute_lat_prop()

        return lat_stable

    def compute_constr(self) -> bool:
        self._eta_centre, _, _, _ = self._lat_prop.get_Twiss(self._uc_centre)
        Twiss_sp = self._lat_prop.get_Twiss(-1)
        _, _, _, self._nu_uc = Twiss_sp

        self._constr = {
            "eps_x":   (self._lat_prop._eps[ind.X]
                        -self._des_val_list["eps_x_des"])**2,
            "nu_uc_x": (self._nu_uc[ind.X]
                        -self._des_val_list["nu_uc_des"][ind.X])**2,
            "nu_uc_y": (self._nu_uc[ind.Y]
                        -self._des_val_list["nu_uc_des"][ind.Y])**2,
            "eta_x":   self._eta_centre[ind.x]**2,
            "xi":      self._xi[ind.X]**2 + self._xi[ind.Y]**2
        }

    def compute_chi_2(self) -> float:
        self.compute_constr()

        prt_len = 11
        prt = not False

        if self._first:
            constr_len = len(self._constr)
            weights_len = len(self._weights)
            if constr_len-weights_len != 0:
                print(f"\ncompute_chi_2() - len(constr) != "
                      f"len(weights): {constr_len:d} {weights_len:d}")
                raise ValueError("Number of constraints != weights.")
            self._first = False

        chi_2 = 0e0
        if prt:
            print()
        for key in self._weights:
            dchi_2 = self._weights[key]*self._constr[key]
            chi_2 += dchi_2
            if prt:
                print(f"  dchi_2({key:s})", end="")
                for k in range(prt_len-len(key)):
                    print(" ", end="")
                print(f" = {dchi_2:9.3e}")

        return chi_2

    def prt_iter(self, prm, chi_2) -> None:
        print(f"\n{self._n_iter:3d} chi_2 = {self._chi_2_min:9.3e}",
              f"({chi_2-self._chi_2_min:9.3e})")
        print(f"    eps_x [pm.rad] = "
              f"{1e12*self._lat_prop._eps[ind.X]:5.3f} "
              f"({1e12*self._des_val_list["eps_x_des"]:5.3f})")

        print(f"    nu_uc          = [{self._nu_uc[ind.X]:7.5f}, "
              f"{self._nu_uc[ind.Y]:7.5f}] "
              f"([{self._des_val_list["nu_uc_des"][ind.X]:7.5f}, "
              f"{self._des_val_list["nu_uc_des"][ind.Y]:7.5f}])")
        print(f"    eta_x          = {self._eta_centre[ind.x]:9.3e}")
        print(f"    xi             = [{self._xi[ind.X]:5.3f}, "
              f"{self._xi[ind.Y]:5.3f}]")

        print(f"\n    phi_sp         = {self._phi_sp:8.5f}")
        print(f"    C [m]          = "
              f"{self._lat_prop.compute_circ():8.5f}")

        print()
        for k, phi in enumerate(self._phi_bend):
                b_2 = \
                    self._bend_list[k].compute_bend_b_2xL() \
                    /self._bend_list[k].compute_bend_L_tot()
                print(f"    phi_bend_{k+1:1d}     = "
                      f"[{phi:8.5f}, {b_2:8.5f}]")
        print()
        for k, rbend in enumerate(self._rbend_list):
            print(f"    {rbend:10s}     = {self._phi_rbend[k]:8.5f}")

        print(f"\n    b_3            =", self._nld._b_3_list)

        self._lat_prop.prt_rad()
        self._prm_list.prt_prm(prm)

    def f_sp(self, prm: np.ndarray) -> float:
        self._n_iter += 1
        self._prm_list.set_prm(prm)

        if not self.compute_prop():
            chi_2 = 1e30
            if not False:
                print(f"\n{self._n_iter:3d} chi_2 = {chi_2:11.5e}",
                      f"({self._chi_2_min:11.5e})")
                self._prm_list.prt_prm(prm)
            return chi_2

        chi_2 = self.compute_chi_2()

        if chi_2 < self._chi_2_min:
            self.prt_iter(prm, chi_2)
            if False:
                self._lat_prop.prt_Twiss("twiss.txt")
            pc.prt_lat(
                self._lat_prop, self._file_name, self._prm_list,
                bend_list=self._bend_list)
            self._chi_2_min = min(chi_2, self._chi_2_min)
        else:
            print(f"\n{self._n_iter:3d} chi_2 = {self._chi_2_min:9.3e}",
                  f"({chi_2-self._chi_2_min:9.3e})")
            if False:
                print(f"\n{self._n_iter:3d}",
                      f"dchi_2 = {chi_2-self._chi_2_min:9.3e}")
                # self._prm_list.prt_prm(prm)

        return chi_2

    def opt_uc(self) -> opt.OptimizeResult:
        """Use Case: optimise super period.
        """

        max_iter = 10000
        f_tol    = 1e-7
        x_tol    = 1e-7
        g_tol    = 1e-7

        prm, bounds = self._prm_list.get_prm()
        if True:
            self._phi_sp_0 = 18.0
        else:
            self._phi_sp_0 = self._lat_prop.compute_phi_lat()
        self.f_sp(prm)

        # Methods:
        #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
        #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.
        # Methods with boundaries:
        #   L-BFGS-B, TNC, and SLSQP.
        # Methods with constraints:
        #   SLSQP.

        opt_dict = {
            "CG": {"options": {
                "gtol": g_tol, "maxiter": max_iter, "eps": self._dprm_list}},
            "TNC": {"options": {
                "ftol": f_tol, "gtol": g_tol, "xtol": x_tol, "maxfun": max_iter,
                "eps": 1e-3}},
            "BFGS": {"options": {
                "ftol": f_tol, "gtol": g_tol, "maxiter": max_iter,
                "eps": 1e-5}},
            "SLSQP": {"options": {
                "ftol": f_tol, "maxiter": max_iter, "eps": 1e-4}}
            }

        method = "SLSQP"
        
        minimum = opt.minimize(
            self.f_sp,
            prm,
            method=method,
            # callback=prt_iter,
            jac=None,
            bounds = bounds,
            options=opt_dict[method]["options"]
        )

        print("\n".join(minimum))


def get_bends():
    bend_lists = {
        "bend_list": [
            ["d2_h2_sl_d0a", "d2_h2_sl_d0b", "d2_h2_sl_d0c", "d2_h2_sl_df1",
             "d2_h2_sl_df2", "d2_h2_sl_df3", "d2_h2_sl_df4",
             "d2_h2_sl_df5"]
        ],
        "rbend_list": ["r2_h2"]
    }

    bend_list = []
    for b_list in bend_lists["bend_list"]:
        bend_list.append(pc.bend_class(lat_prop, b_list))

    return bend_list, bend_lists["rbend_list"]


def get_prms(bend_list, eps):
    prm = [
        (bend_list[0], "b_2_bend", prm_range["b_2_bend"]),
        ("r2_h2", "b_2", prm_range["b_2"]),
        ("r2_h2", "phi", [-0.6,  0.0])
    ]

    prm_list = pc.prm_class(lat_prop, prm)

    dprm_list = np.full(len(prm), eps)

    return prm_list, dprm_list


def get_weights():
    weights = {
        "eps_x"   : 1e16,
        "nu_uc_x" : 1e0,
        "nu_uc_y" : 1e0,
        "eta_x"   : 1e-1,
        "xi"      : 1e-6
    }
    return weights


cod_eps   = 1e-15
E_0       = 3.0e9

A_max     = np.array([6e-3, 3e-3])
delta_max = 3e-2
beta_inj  = np.array([3.0, 3.0])

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
file_name = os.path.join(home_dir, sys.argv[1]+".lat")

lat_prop = lp.lattice_properties_class(file_name, E_0, cod_eps, 2)

lat_prop.prt_lat("max_4u_uc_lat.txt")

uc_centre = lat_prop._lattice.find("d2_h2_sl_d0a", 0).index
print("\nunit cell centre {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_centre].name, uc_centre))

np.set_printoptions(formatter={"float": "{:10.3e}".format})

b_3_list = ["s1_h3", "s2_h3", "s3_h3", "s4_h3"]
nld = nld_class.nonlin_dyn_class(lat_prop, A_max, beta_inj, delta_max, b_3_list)

nld.zero_mult(3)
nld.zero_mult(4)

# Compute Twiss parameters along lattice.
if not lat_prop.comp_per_sol():
    print("\ncomp_per_sol: unstable")
    exit(1)
# Compute radiation properties.
if not lat_prop.compute_radiation():
    print("\ncompute_radiation: unstable")
    exit(1)
lat_prop.prt_lat_param()
lat_prop.prt_rad()
lat_prop.prt_M()
lat_prop.prt_M_rad()

weight_list = get_weights()
bend_list, rbend_list = get_bends()
prm_list, dprm_list = get_prms(bend_list, 1e-4)

@dataclass
class prm_class:
    lattice:     ClassVar[list] = [lat_prop, nld]
    s_loc:       ClassVar[list] = uc_centre
    design_vals: ClassVar[dict] = design_val
    weights:     ClassVar[list] = weight_list
    dipoles:     ClassVar[list] = [bend_list, rbend_list]
    params:      ClassVar[list] = [prm_list, dprm_list]

opt_uc = opt_uc_class(prm_class)
opt_uc.opt_uc()
