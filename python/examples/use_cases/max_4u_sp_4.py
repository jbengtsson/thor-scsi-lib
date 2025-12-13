"""Use Case:
     Optimisation of a higher-order-achromat for negative linear momentum
     compaction.
     Use one super period.
"""


import os
import sys
from dataclasses import dataclass
from typing import ClassVar

import math
import numpy as np
from scipy import optimize as opt

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, prm_class as pc, nonlin_dyn as nld_class
from thor_scsi.utils.output import mat2txt, vec2txt


ind = ind.index_class()

prm_range = {
    "phi":       [  0.0,  2.0],
    "phi_rbend": [ -0.5,  0.0],
    "b_2":       [-10.0, 10.0],
    "b_2_bend":  [ -1.5,  1.5]
}

design_val = {
    "eps_x_des"    : 35e-12,
    "phi_1_des"    : 1.3,
    "phi_rb_1_des" : -0.2,
    "b_2_des"      : 2.0,
    "nu_uc_des"    : np.array([3.0/7.0, 1.0/7.0]),
    # "eta_x_uc_des" : 0.0,
    "eta_x_uc_des" : -1e-3,
    "nu_sp_des"    : np.array([58.20/20.0,  17.28/20.0]),
    # "nu_sp_des"    : np.array([58.14/20.0,  17.27/20.0]),
    # "nu_sp_des"    : np.array([57.14/20.0,  20.27/20.0]),
    # "nu_sp_des"    : np.array([58.17/20.0,  17.22/20.0]),
    "beta_des"     : [5.0, 3.0],
    # Phase advance across the straight.
    "dnu_des"      : [0.75, 0.25]
}


class opt_sp_class:
    # Private.

    def __init__(self, prm_class: pc.prm_class) -> None:
        self._first          = True

        self._lat_prop       = prm_class.lattice[0]
        self._nld            = prm_class.lattice[1]
        self._uc_list        = prm_class.s_loc[0]
        self._sp_list        = prm_class.s_loc[1]
        self._des_val_list   = prm_class.design_vals
        self._weights        = prm_class.weights
        self._bend_list      = prm_class.dipoles[0]
        self._rbend_list     = prm_class.dipoles[1]
        self._prm_list       = prm_class.params[0]
        self._dprm_list      = prm_class.params[1]

        self._phi_bend       = np.zeros(len(self._bend_list))
        self._phi_rbend      = np.zeros(len(self._bend_list))

        self._b_2            = np.nan
        self._phi_sp_0       = np.nan
        self._phi_sp         = np.nan
        self._dphi           = np.nan

        self._Twiss_sp       = None
        self._eta_list       = np.zeros((2, 4))
        self._alpha_list     = np.zeros((2, 2))

        self._nu_uc          = np.nan
        self._alpha_c        = np.nan
        self._nu_0           = np.nan
        self._dnu            = np.nan
        self._xi             = np.nan

        self._eta_nl         = lat_prop.tpsa()

        self._constr         = {}

        self._chi_2_min      = 1e30
        self._n_iter         = -1
        self._file_name      = "opt_sp.txt"

        self._h_im_scl_rms_0 = np.nan
        self._h_im_scl_rms_1 = np.nan
        self._K_re_scl_rms_0 = np.nan
        self._K_re_scl_rms_1 = np.nan

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

        # Compute 2nd order dispersion
        self._eta_nl = nld.compute_eta(self._lat_prop._M)
        return True

    def compute_prop(self) -> bool:
        self._phi_sp = self._lat_prop.compute_phi_lat()
        self._dphi = self._phi_sp - self._phi_sp_0

        for k, bend in enumerate(self._bend_list):
            self._phi_bend[k] = bend.compute_bend_phi()

        for k, bend in enumerate(self._rbend_list):
            self._phi_rbend[k] = self._lat_prop.get_phi_elem(bend, 0)

        # self._b_2 = \
        #     self._bend_list[0].compute_bend_b_2xL() \
        #     /self._bend_list[0].compute_bend_L_tot()

        self._b_2 = self._lat_prop.get_b_n_elem("s1_h2", 0, 2)

        lat_stable = self.compute_lat_prop()

        return lat_stable

    def compute_constr(self) -> bool:
        self._alpha_c = self._lat_prop._alpha_c
        self._eta_list[0], self._alpha_list[0], _, self._nu_0 = \
            self._lat_prop.get_Twiss(self._uc_list[0])
        self._eta_list[1], self._alpha_list[1], _, self._nu_1 = \
            self._lat_prop.get_Twiss(self._uc_list[1])
        self._nu_uc = self._nu_1 - self._nu_0
        self._eta_list[1], self._alpha_list[1], _, _ = \
            self._lat_prop.get_Twiss(self._uc_list[1])
        self._Twiss_sp = self._lat_prop.get_Twiss(-1)
        self._dnu = 2e0*self._lat_prop.get_Twiss(sp_list[0])[3][:]

        eta, _, beta_sp, nu_sp = self._Twiss_sp

        self._constr = {
            "eps_x":       (self._lat_prop._eps[ind.X]
                            -self._des_val_list["eps_x_des"])**2,
            "dphi":        self._dphi**2,
            "phi_1":       (self._phi_bend[0]
                            -self._des_val_list["phi_1_des"])**2,
            "phi_rb":      (self._phi_rbend[0]
                            -self._des_val_list["phi_rb_1_des"])**2,
            "b_2":         (self._b_2-self._des_val_list["b_2_des"])**2,
            "alpha^(1)_c": 1e0/self._alpha_c[1]**2,
            "alpha^(2)_c": self._alpha_c[2]**2,
            "U_0":         self._lat_prop._U_0**2,
            "eta_x_uc":    ((self._eta_list[0][ind.x]
                             -self._des_val_list["eta_x_uc_des"])**2
                            +(self._eta_list[1][ind.x]
                              -self._des_val_list["eta_x_uc_des"]))**2,
            "eta'_x_uc":   (self._eta_list[0][ind.px]**2
                            +self._eta_list[1][ind.px]**2),
            "alpha_uc":    (self._alpha_list[0][ind.X]**2
                            + self._alpha_list[0][ind.Y]**2
                            + self._alpha_list[1][ind.X]**2
                            + self._alpha_list[1][ind.Y]**2),
            "nu_uc_x":    (self._nu_uc[ind.X]
                           -self._des_val_list["nu_uc_des"][ind.X])**2,
            "nu_uc_y":    (self._nu_uc[ind.Y]
                           -self._des_val_list["nu_uc_des"][ind.Y])**2,
            "eta_x":      eta[ind.x]**2,
            "nu_sp_x":    (nu_sp[ind.X]
                           -self._des_val_list["nu_sp_des"][ind.X])**2,
            "nu_sp_y":    (nu_sp[ind.Y]
                           -self._des_val_list["nu_sp_des"][ind.Y])**2,
            "beta_x":     (beta_sp[ind.X]
                           -self._des_val_list["beta_des"][ind.X])**2,
            "beta_y":     (beta_sp[ind.Y]
                           -self._des_val_list["beta_des"][ind.Y])**2,
            "dnu_x":      (self._dnu[ind.X]
                           -self._des_val_list["dnu_des"][ind.X])**2,
            "dnu_y":      (self._dnu[ind.Y]
                           -self._des_val_list["dnu_des"][ind.Y])**2,
            "xi":         self._xi[ind.X]**2 + self._xi[ind.Y]**2,
            "eta^(2)_x":  self._eta_nl.get([0, 0, 0, 0, 2])**2,
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
        eta, alpha, beta, nu_sp = self._Twiss_sp

        print(f"\n{self._n_iter:3d} chi_2 = {self._chi_2_min:9.3e}",
              f"({chi_2-self._chi_2_min:9.3e})")
        print(f"    eps_x [pm.rad] = "
              f"{1e12*self._lat_prop._eps[ind.X]:5.3f} "
              f"({1e12*self._des_val_list["eps_x_des"]:5.3f})")

        print(f"\n    dphi [deg]     = {self._dphi:9.3e}")

        print("\n    alpha_c (multipoles zeroed)\n"
              f"                   = [{self._alpha_c[1]:9.3e}, "
              f"{self._alpha_c[2]:9.3e}]")
        print(f"    nu_uc          = [{self._nu_uc[ind.X]:7.5f}, "
              f"{self._nu_uc[ind.Y]:7.5f}] "
              f"([{self._des_val_list["nu_uc_des"][ind.X]:7.5f}, "
              f"{self._des_val_list["nu_uc_des"][ind.Y]:7.5f}])")
        print(f"    nu_sp          = [{nu_sp[ind.X]:7.5f}, "
              f"{nu_sp[ind.Y]:7.5f}] "
              f"([{self._des_val_list["nu_sp_des"][ind.X]:7.5f}, "
              f"{self._des_val_list["nu_sp_des"][ind.Y]:7.5f}])")
        print(f"    dnu            = [{self._dnu[ind.X]:7.5f}, "
              f"{self._dnu[ind.Y]:7.5f}] "
              f"([{self._des_val_list["dnu_des"][ind.X]:7.5f}, "
              f"{self._des_val_list["dnu_des"][ind.Y]:7.5f}])")
        print(f"    xi             = [{self._xi[ind.X]:5.3f}, "
              f"{self._xi[ind.Y]:5.3f}]")
        print(f"    eta^(2)_x      = "
              f"{self._eta_nl.get([0, 0, 0, 0, 2]):9.3e}")

        print(f"\n    eta_uc         = "
              f"{self._eta_list[0][ind.x]:9.3e}"
              f" {self._eta_list[1][ind.x]:9.3e}")
        print(f"    eta'_uc        = "
              f"{self._eta_list[0][ind.px]:9.3e}"
              f" {self._eta_list[1][ind.px]:9.3e}")
        print(f"    alpha_uc       = "
              f"[{self._alpha_list[0][ind.X]:9.3e}, "
              f"{self._alpha_list[0][ind.Y]:9.3e}] "
              f"[{self._alpha_list[1][ind.X]:9.3e}, "
              f"{self._alpha_list[1][ind.Y]:9.3e}]")
        print(f"    eta_x          = {eta[ind.x]:9.3e}")
        print(f"    beta           = "
              f"[{beta[ind.X]:7.5f}, {beta[ind.Y]:7.5f}] "
              f"([{self._des_val_list["beta_des"][ind.X]:7.5f}, "
              f"{self._des_val_list["beta_des"][ind.Y]:7.5f}])")

        print(f"\n    phi_sp         = {self._phi_sp:8.5f}")
        print(f"    C [m]          = "
              f"{self._lat_prop.compute_circ():8.5f}")

        print()
        for k, phi in enumerate(self._phi_bend):
                self._b_2 = \
                    self._bend_list[k].compute_bend_b_2xL() \
                    /self._bend_list[k].compute_bend_L_tot()
                print(f"    phi_bend_{k+1:1d}     = "
                      f"[{phi:8.5f}, {self._b_2:8.5f}]")
        print()
        for k, rbend in enumerate(self._rbend_list):
            print(f"    {rbend:10s}     = {self._phi_rbend[k]:8.5f}")

        print(f"\n    b_3            =", self._nld._b_3_list)
        print(f"\n    b_2            = {self._b_2:8.5f}")


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

    def opt_sp(self) -> opt.OptimizeResult:
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


def get_bends(dip_type, lat_prop):
    if dip_type == 1:
        # max_4u/m4U_250316_h03_01_01_01_tracy-2.
        bend_lists = {
            "bend_list": [
                ["d1_h3_sl_dm5", "d1_h3_sl_dm4", "d1_h3_sl_dm3", "d1_h3_sl_dm2",
                 "d1_h3_sl_dm1", "d1_h3_sl_d0", "d1_h3_sl_ds1", "d1_h3_sl_ds2",
                 "d1_h3_sl_ds3", "d1_h3_sl_ds4", "d1_h3_sl_ds5"],
                ["d2_h3_sl_df0", "d2_h3_sl_df1", "d2_h3_sl_df2", "d2_h3_sl_df3",
                 "d2_h3_sl_df4", "d2_h3_sl_df5", "d2_h3_sl_df6"],
                ["d3_h3_sl_df0", "d3_h3_sl_df1", "d3_h3_sl_df2", "d3_h3_sl_df3",
                 "d3_h3_sl_df4", "d3_h3_sl_df5", "d3_h3_sl_df6"]
            ],
            "rbend_list": ["r1_h3", "r2_h3", "r3_h3"]
        }
    elif dip_type == 2:
        # m4U_241223_h02_01_01_01_tracy_2.
        bend_lists = {
            "bend_list": [
                ["d1_h3_sl_dm5", "d1_h3_sl_dm4", "d1_h3_sl_dm3", "d1_h3_sl_dm2",
                 "d1_h3_sl_dm1", "d1_h3_sl_ds0", "d1_h3_sl_ds1", "d1_h3_sl_ds2",
                 "d1_h3_sl_ds3", "d1_h3_sl_ds4", "d1_h3_sl_ds5", "d1_h3_sl_ds6"
                 ],
                ["d2_h3_sl_d0a", "d2_h3_sl_d0b", "d2_h3_sl_d0c", "d2_h3_sl_df1",
                 "d2_h3_sl_df2", "d2_h3_sl_df3", "d2_h3_sl_df4",
                 "d2_h3_sl_df5"],
            ],
            "rbend_list": ["r1_h3", "r2_h3"]
        }
    elif dip_type == 3:
        # m4U_250527_h02_13_02_01_tracy-2.
        if True:
            bend_lists = {
                "bend_list": [
                    ["d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3",
                     "d1_h2_sl_dm2", "d1_h2_sl_dm1", "d1_h2_sl_ds0",
                     "d1_h2_sl_ds1", "d1_h2_sl_ds2", "d1_h2_sl_ds3",
                     "d1_h2_sl_ds4", "d1_h2_sl_ds5", "d1_h2_sl_ds6"],
                    ["d2_h2_sl_d0a", "d2_h2_sl_d0b", "d2_h2_sl_d0c",
                     "d2_h2_sl_df1", "d2_h2_sl_df2", "d2_h2_sl_df3",
                     "d2_h2_sl_df4", "d2_h2_sl_df5"]
                ],
                "rbend_list": ["r1_h2", "r2_h2"]
            }
        else:
            bend_lists = {
                "bend_list": [
                    ["d1_h2_sl_dm5", "d1_h2_sl_dm3", "d1_h2_sl_dm2",
                     "d1_h2_sl_ds0", "d1_h2_sl_ds1", "d1_h2_sl_ds2",
                     "d1_h2_sl_ds3", "d1_h2_sl_ds5"],
                    ["d2_h2_sl_d0a", "d2_h2_sl_d0c", "d2_h2_sl_df1",
                     "d2_h2_sl_df3", "d2_h2_sl_df4", "d2_h2_sl_df5"]
                ],
                "rbend_list": ["r1_h2", "r2_h2"]
            }
    elif dip_type == 4:
        # m4U_250610_h02_16_02_01_tracy_2.
        bend_lists = {
            "bend_list": [
                ["d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3", "d1_h2_sl_dm2",
                 "d1_h2_sl_dm1", "d1_h2_sl_ds0", "d1_h2_sl_ds1", "d1_h2_sl_ds2",
                 "d1_h2_sl_ds3", "d1_h2_sl_ds4", "d1_h2_sl_ds5",
                 "d1_h2_sl_ds6"],
                ["d2_h2_sl_d0a", "d2_h2_sl_d0b", "d2_h2_sl_d0c", "d2_h2_sl_df1",
                 "d2_h2_sl_df2", "d2_h2_sl_df3", "d2_h2_sl_df4",
                 "d2_h2_sl_df5"],
                ["d3_h2_sl_d0a", "d3_h2_sl_d0b", "d3_h2_sl_d0c", "d3_h2_sl_df1",
                 "d3_h2_sl_df2", "d3_h2_sl_df3", "d3_h2_sl_df4", "d3_h2_sl_df5"]
            ],
            "rbend_list": ["r1_h2", "r2_h2", "r3_h2"]
        }

    bend_list = []
    for b_list in bend_lists["bend_list"]:
        bend_list.append(pc.bend_class(lat_prop, b_list))

    return bend_list, bend_lists["rbend_list"]


def get_prms(prm_type, lat_prop, bend_list, eps):
    if prm_type == 1:
        prm = [
            ("q1_h2",        "b_2",      prm_range["b_2"]),
            ("q2_h2",        "b_2",      prm_range["b_2"]),

            (bend_list[0],   "b_2_bend", prm_range["b_2_bend"]),
            (bend_list[1],   "b_2_bend", prm_range["b_2_bend"]),
            (bend_list[2],   "b_2_bend", prm_range["b_2_bend"]),

            (bend_list[0],   "phi_bend", prm_range["phi"]),
            (bend_list[1],   "phi_bend", prm_range["phi"]),
            (bend_list[2],   "phi_bend", prm_range["phi"]),

            ("r1_h2",        "b_2",      prm_range["b_2"]),
            ("r2_h2",        "b_2",      prm_range["b_2"]),
            ("r3_h2",        "b_2",      prm_range["b_2"]),

            ("r1_h2",        "phi",      prm_range["phi_rbend"]),
            ("r2_h2",        "phi",      prm_range["phi_rbend"]),
            ("r3_h2",        "phi",      prm_range["phi_rbend"])
        ]
    elif prm_type == 2:
        prm = [
            ("q1_h2",        "b_2",      prm_range["b_2"]),
            ("q2_h2",        "b_2",      prm_range["b_2"]),
            ("r1_h2",        "b_2",      prm_range["b_2"]),
            ("lego",         "phi",      [-1.0, 1.0]),
            ("lego",         "b_2",      prm_range["b_2"]),

            (bend_list[0],   "b_2_bend", prm_range["b_2_bend"]),
            (bend_list[1],   "b_2_bend", prm_range["b_2_bend"]),

            (bend_list[0],   "phi_bend", [0.5, 1.5]),
            (bend_list[1],   "phi_bend", [1.5, 3.0]),

            ("r2_h2",        "b_2",      prm_range["b_2"]),

            ("r2_h2",        "phi",      prm_range["phi_rbend"])
        ]
    elif prm_type == 3:
        prm = [
            ("q1_h2",        "b_2", prm_range["b_2"]),
            ("q2_h2",        "b_2", prm_range["b_2"]),
            ("q1_h2",        "phi", [-1.0, 1.0]),
            ("q2_h2",        "phi", [-1.0, 1.0]),

            ("r1_h2",        "b_2", prm_range["b_2"]),

            ("lego",         "phi", [-1.0, 1.0]),
            ("lego",         "b_2", prm_range["b_2"]),

            ("d1_h2_sl_dm5", "b_2", prm_range["b_2"]),
            ("d1_h2_sl_dm3", "b_2", prm_range["b_2"]),
            ("d1_h2_sl_dm2", "b_2", prm_range["b_2"]),
            ("d1_h2_sl_ds0", "b_2", prm_range["b_2"]),
            ("d1_h2_sl_ds1", "b_2", prm_range["b_2"]),
            ("d1_h2_sl_ds2", "b_2", prm_range["b_2"]),
            ("d1_h2_sl_ds3", "b_2", prm_range["b_2"]),
            ("d1_h2_sl_ds5", "b_2", prm_range["b_2"]),

            ("d1_h2_sl_dm5", "phi", [-1.5, 1.5]),
            ("d1_h2_sl_dm3", "phi", [-1.5, 1.5]),
            ("d1_h2_sl_dm2", "phi", [-1.5, 1.5]),
            ("d1_h2_sl_ds0", "phi", [-1.5, 1.5]),
            ("d1_h2_sl_ds1", "phi", [-1.5, 1.5]),
            ("d1_h2_sl_ds2", "phi", [-1.5, 1.5]),
            ("d1_h2_sl_ds3", "phi", [-1.5, 1.5]),
            ("d1_h2_sl_ds5", "phi", [-1.5, 1.5]),
 
            ("d2_h2_sl_d0a", "b_2", prm_range["b_2"]),
            ("d2_h2_sl_d0c", "b_2", prm_range["b_2"]),
            ("d2_h2_sl_df1", "b_2", prm_range["b_2"]),
            ("d2_h2_sl_df3", "b_2", prm_range["b_2"]),
            ("d2_h2_sl_df4", "b_2", prm_range["b_2"]),
            ("d2_h2_sl_df5", "b_2", prm_range["b_2"]),

            ("d2_h2_sl_d0a", "phi", [-1.5, 1.5]),
            ("d2_h2_sl_d0c", "phi", [-1.5, 1.5]),
            ("d2_h2_sl_df1", "phi", [-1.5, 1.5]),
            ("d2_h2_sl_df3", "phi", [-1.5, 1.5]),
            ("d2_h2_sl_df4", "phi", [-1.5, 1.5]),
            ("d2_h2_sl_df5", "phi", [-1.5, 1.5]),

            ("r2_h2",        "b_2", prm_range["b_2"]),
            ("r2_h2",        "phi", prm_range["phi_rbend"])
        ]
    elif prm_type == 4:
        prm = [
            ("q1_h2",        "b_2",      prm_range["b_2"]),
            ("q2_h2",        "b_2",      prm_range["b_2"]),

            (bend_list[0],   "b_2_bend", prm_range["b_2_bend"]),
            (bend_list[1],   "b_2_bend", prm_range["b_2_bend"]),
            (bend_list[2],   "b_2_bend", prm_range["b_2_bend"]),

            (bend_list[0],   "phi_bend", [1.4,  1.5]),
            (bend_list[1],   "phi_bend", [1.5,  3.0]),
            (bend_list[2],   "phi_bend", [1.5,  2.5]),

            ("r1_h2",        "b_2",      prm_range["b_2"]),
            ("r2_h2",        "b_2",      prm_range["b_2"]),

            ("r1_h2",        "phi",      prm_range["phi_rbend"]),
            ("r2_h2",        "phi",      prm_range["phi_rbend"]),
        ]
    elif prm_type == 5:
        prm = [
            ("q1_h3",        "b_2",      prm_range["b_2"]),
            ("q2_h3",        "b_2",      prm_range["b_2"]),

            (bend_list[0],   "b_2_bend", prm_range["b_2_bend"]),
            (bend_list[1],   "b_2_bend", prm_range["b_2_bend"]),
            (bend_list[2],   "b_2_bend", prm_range["b_2_bend"]),

            ("d1_h3_sl_dm5", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_dm4", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_dm3", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_dm2", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_dm1", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_d0",  "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds1", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds2", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds3", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds4", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds5", "phi",      prm_range["b_2_bend"]),

            ("d2_h3_sl_df0", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df1", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df2", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df3", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df4", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df5", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df6", "phi",      prm_range["b_2_bend"]),

            ("d3_h3_sl_df0", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df1", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df2", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df3", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df4", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df5", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df6", "phi",      prm_range["b_2_bend"]),

            ("r1_h3",        "b_2",      prm_range["b_2"]),
            ("r2_h3",        "b_2",      prm_range["b_2"]),
            ("r3_h3",        "b_2",      prm_range["b_2"]),

            ("r1_h3",        "phi",      [-0.2,  0.2]),
            ("r2_h3",        "phi",      [-0.2,  0.0]),
            ("r3_h3",        "phi",      [-0.2,  0.0])
        ]
    elif prm_type == 6:
        prm = [
            ("q1_h3",        "b_2",      prm_range["b_2"]),
            ("q2_h3",        "b_2",      prm_range["b_2"]),

            ("d1_h3_sl_dm5", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_dm4", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_dm3", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_dm2", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_dm1", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_d0",  "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_ds1", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_ds2", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_ds3", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_ds4", "b_2",      prm_range["b_2"]),
            ("d1_h3_sl_ds5", "b_2",      prm_range["b_2"]),

            ("d2_h3_sl_df0", "b_2",      prm_range["b_2"]),
            ("d2_h3_sl_df1", "b_2",      prm_range["b_2"]),
            ("d2_h3_sl_df2", "b_2",      prm_range["b_2"]),
            ("d2_h3_sl_df3", "b_2",      prm_range["b_2"]),
            ("d2_h3_sl_df4", "b_2",      prm_range["b_2"]),
            ("d2_h3_sl_df5", "b_2",      prm_range["b_2"]),
            ("d2_h3_sl_df6", "b_2",      prm_range["b_2"]),

            ("d3_h3_sl_df0", "b_2",      prm_range["b_2"]),
            ("d3_h3_sl_df1", "b_2",      prm_range["b_2"]),
            ("d3_h3_sl_df2", "b_2",      prm_range["b_2"]),
            ("d3_h3_sl_df3", "b_2",      prm_range["b_2"]),
            ("d3_h3_sl_df4", "b_2",      prm_range["b_2"]),
            ("d3_h3_sl_df5", "b_2",      prm_range["b_2"]),
            ("d3_h3_sl_df6", "b_2",      prm_range["b_2"]),

            ("d1_h3_sl_dm5", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_dm4", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_dm3", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_dm2", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_dm1", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_d0",  "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds1", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds2", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds3", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds4", "phi",      prm_range["b_2_bend"]),
            ("d1_h3_sl_ds5", "phi",      prm_range["b_2_bend"]),

            ("d2_h3_sl_df0", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df1", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df2", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df3", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df4", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df5", "phi",      prm_range["b_2_bend"]),
            ("d2_h3_sl_df6", "phi",      prm_range["b_2_bend"]),

            ("d3_h3_sl_df0", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df1", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df2", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df3", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df4", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df5", "phi",      prm_range["b_2_bend"]),
            ("d3_h3_sl_df6", "phi",      prm_range["b_2_bend"]),

            ("r1_h3",        "b_2",      prm_range["b_2"]),
            ("r2_h3",        "b_2",      prm_range["b_2"]),
            ("r3_h3",        "b_2",      prm_range["b_2"]),

            ("r1_h3",        "phi",      [-0.2,  0.2]),
            ("r2_h3",        "phi",      [-0.2,  0.0]),
            ("r3_h3",        "phi",      [-0.2,  0.0])
        ]

    prm_list = pc.prm_class(lat_prop, prm)

    dprm_list = np.full(len(prm), eps)

    return prm_list, dprm_list


def get_weights():
    weights = {
        "eps_x"       : 1e1*1e17,
        "dphi"        : 1e-1, 
        "phi_1"       : 0e-2,  
        "phi_rb"      : 0e-3,  
        "b_2"         : 0e-3, 
        "alpha^(1)_c" : 1e1*1e-14,  
        "alpha^(2)_c" : 1e1,
        "U_0"         : 1e-15,
        "eta_x_uc"    : 1e1*1e2, 
        "eta'_x_uc"   : 1e2, 
        "alpha_uc"    : 1e-1,
        "nu_uc_x"     : 1e-2,
        "nu_uc_y"     : 1e-2,
        "eta_x"       : 1e2,
        "nu_sp_x"     : 0e0, 
        "nu_sp_y"     : 0e0,
        "beta_x"      : 0e-6,
        "beta_y"      : 0e-6,
        "dnu_x"       : 1e-3,
        "dnu_y"       : 1e-3,
        "xi"          : 1e-1*1e-6,
        "eta^(2)_x"   : 1e-6 
    }
    return weights


no = 2

cod_eps = 1e-10
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
delta_max = 3e-2
beta_inj  = np.array([3.0, 3.0])

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")

file_name = os.path.join(home_dir, sys.argv[1]+".lat")

lat_prop = lp.lattice_properties_class(file_name, E_0, cod_eps, no)
lat_prop.prt_lat("lat_prop_lat.txt")

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

dip_type = 3
prm_type = 3

# uc_list = np.array(lat_prop._lattice.find("d2_h3_sl_df0", 0).index)
# uc_list = np.append(uc_list, lat_prop._lattice.find("d3_h3_sl_df0", 1).index)
if dip_type == 3:
    uc_list = np.array(lat_prop._lattice.find("d2_h2_sl_d0a", 0).index)
    uc_list = np.append(uc_list, lat_prop._lattice.find("d2_h2_sl_d0a", 2)
                        .index)
else:
    uc_list = np.array(lat_prop._lattice.find("d2_h2_sl_d0a", 0).index)
    uc_list = np.append(uc_list, lat_prop._lattice.find("d3_h2_sl_d0a", 1)
                        .index)

sp_list = np.array(lat_prop._lattice.find("lsborder", 0).index)
sp_list = np.append(sp_list, lat_prop._lattice.find("lsborder", 1).index)

print(f"\nunit cell entrance          ",
      f"{lat_prop._lattice[uc_list[0]].name:15s} loc = {uc_list[0]:d}")
print(f"unit cell exit               {lat_prop._lattice[uc_list[1]].name:15s}",
      f"loc = {uc_list[1]:d}")
print(f"super period first sextupole {lat_prop._lattice[sp_list[0]].name:15s}",
      f"loc = {sp_list[0]:d}")
print(f"super period last sextupole  {lat_prop._lattice[sp_list[1]].name:15s}",
      f"loc = {sp_list[1]:d}")

weight_list = get_weights()
bend_list, rbend_list = get_bends(dip_type, lat_prop)
prm_list, dprm_list = get_prms(prm_type, lat_prop, bend_list, 1e-4)

@dataclass
class prm_class:
    lattice:     ClassVar[list] = [lat_prop, nld]
    s_loc:       ClassVar[list] = [uc_list, sp_list]
    design_vals: ClassVar[dict] = design_val
    weights:     ClassVar[list] = weight_list
    dipoles:     ClassVar[list] = [bend_list, rbend_list]
    params:      ClassVar[list] = [prm_list, dprm_list]

opt_sp = opt_sp_class(prm_class)
opt_sp.opt_sp()
