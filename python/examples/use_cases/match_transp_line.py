"""Use Case:
     Transport line linear optics matching.
"""


import os
import sys
import copy as _copy
from dataclasses import dataclass
from typing import ClassVar

import gtpsa

import math
import numpy as np
from scipy import optimize as opt

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, prm_class as pc, courant_snyder as cs, \
    nonlin_dyn as nld_class
from thor_scsi.utils.output import mat2txt, vec2txt


ind = ind.index_class()


# Desired linear optics at the exit: alpha, beta, [eta_x, eta'_x].
Twiss_design = np.array([[0.0, 0.0], [0.0, 0.0], [10.0, 10.0]])

# Linear optics at the entrance.
Twiss_entrance = np.array([[5.19744e-02, 0.0], [0.0, 0.0], [4.28596, 2.43637]])

weights = np.array([[1e6, 1e7], [1e2, 1e2], [1e-4, 1e-4]]) 


class transp_line_class:
    # Private.

    def __init__(self, prm_class: pc.prm_class) -> None:
        self._first       = True

        self._lat_prop    = prm_class.lattice
        self._A_0         = prm_class.A_0
        self._weights     = prm_class.weights
        self._prm_list    = prm_class.params[0]
        self._dprm_list   = prm_class.params[1]

        self._b_2         = np.nan

        self._eta_1       = np.zeros((4))
        self._alpha_1     = np.zeros((2))
        self._beta_1      = np.zeros((2))

        self._Twiss_k_len = np.nan
        self._Twiss_k     = {}

        self._chi_2_min   = 1e30
        self._n_iter      = -1
        self._file_name   = "match_prm.txt"

    # Public.

    def compute_Twiss(self):
        # Compute Twiss parameters along the lattice.
        A_0 = _copy.copy(self._A_0)
        lat_prop._lattice.propagate(
            lat_prop._model_state, A_0, 0, len(lat_prop._lattice))
        # Returns Twiss parameters [eta, alpha, beta] and dnu.
        self._eta_1, self._alpha_1, self._beta_1, _ = \
            cs.compute_Twiss_A(A_0.jacobian())
        self._Twiss_k = np.array([self._eta_1[:2], self._alpha_1, self._beta_1])

    def compute_chi_2(self) -> float:
        prt_len = 11
        prt = not False

        if self._first:
            self._Twiss_k_len = len(self._Twiss_k)
            weights_len = len(self._weights)
            if self._Twiss_k_len-weights_len != 0:
                print(f"\ncompute_chi_2() - len(constr) != "
                      f"len(weights): {self._constr_len:d} {weights_len:d}")
                raise ValueError(
                    f"\ncompute_chi_2() - Number of constraints"
                    f" ({self._Twiss_k_len:d}) !="
                    f" number of weights ({weights_len:d}).")
            self._first = False

        chi_2 = 0e0
        if prt:
            print()
        for k in range(self._Twiss_k_len):
            dchi_2 = np.sum(self._weights[k]*
                (self._Twiss_k[k]-Twiss_design[k])**2)
            chi_2 += dchi_2
            if prt:
                print(f"  dchi_2 = {dchi_2:9.3e}")

        return chi_2

    def prt_iter(self, prm, chi_2) -> None:
        print(f"\n{self._n_iter:3d} chi_2 = {self._chi_2_min:9.3e}",
              f"({chi_2-self._chi_2_min:9.3e})")
        print(f"    [eta_x, eta'_x] = "
              f"[{self._eta_1[ind.x]:9.3e}, {self._eta_1[ind.px]:9.3e}]")
        print(f"    alpha           = "
              f"[{self._alpha_1[ind.X]:7.5f}, {self._alpha_1[ind.Y]:7.5f}]")
        print(f"    beta            = "
              f"[{self._beta_1[ind.X]:7.5f}, {self._beta_1[ind.Y]:7.5f}]")

        self._prm_list.prt_prm(prm)

    def f_match(self, prm: np.ndarray) -> float:
        self._n_iter += 1
        self._prm_list.set_prm(prm)

        self.compute_Twiss()

        chi_2 = self.compute_chi_2()

        if chi_2 < self._chi_2_min:
            self.prt_iter(prm, chi_2)
            if False:
                self._lat_prop.prt_Twiss("twiss.txt")
            pc.prt_lat(self._lat_prop, self._file_name, self._prm_list)
            self._chi_2_min = min(chi_2, self._chi_2_min)
        else:
            print(f"\n{self._n_iter:3d} chi_2 = {self._chi_2_min:9.3e}",
                  f"({chi_2-self._chi_2_min:9.3e})")
            if False:
                print(f"\n{self._n_iter:3d}",
                      f"dchi_2 = {chi_2-self._chi_2_min:9.3e}")
                # self._prm_list.prt_prm(prm)

        return chi_2

    def match(self) -> opt.OptimizeResult:
        max_iter = 10000
        f_tol    = 1e-5
        x_tol    = 1e-7
        g_tol    = 1e-7

        prm, bounds = self._prm_list.get_prm()
        if True:
            self._phi_sp_0 = 18.0
        else:
            self._phi_sp_0 = self._lat_prop.compute_phi_lat()
        self.f_match(prm)

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
            self.f_match,
            prm,
            method=method,
            # callback=prt_iter,
            jac=None,
            bounds = bounds,
            options=opt_dict[method]["options"]
        )

        print("\n".join(minimum))


def get_prms(eps):
    prm = [
        ("D2",  "L",   [0.1-0.05, 0.1+0.1]),
        ("D3",  "L",   [0.1-0.07, 0.1+0.05]),
        ("D4",  "L",   [0.25-0.1, 0.25+0.1]),

        # ("B2",  "phi", [2.17-1.0, 2.17+0.7]),
        # ("B3",  "phi", [0.59-0.59, 0.59+0.8]),

        ("QF2", "b_2", [-10.0, 10.0]),
        ("B2",  "b_2", [-10.0, 10.0]),
        ("B3",  "b_2", [-10.0, 10.0]),
        ("QD1", "b_2", [-10.0, 10.0]),
        ("QF3", "b_2", [-10.0, 10.0])
    ]

    prm_list = pc.prm_class(lat_prop, prm)

    dprm_list = np.full(len(prm), eps)

    return prm_list, dprm_list


no      = 1
cod_eps = 1e-10
E_0     = 2.5e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB")

file_name = os.path.join(home_dir, sys.argv[1]+".lat")

lat_prop = lp.lattice_properties_class(file_name, E_0, cod_eps, no)
lat_prop.prt_lat("lat_prop_lat.txt")

A_0 = lo.compute_A(
    Twiss_entrance[0], Twiss_entrance[1], Twiss_entrance[2], lat_prop._desc)

prm_list, dprm_list = get_prms(1e-4)

@dataclass
class prm_class:
    lattice:      ClassVar[list] = lat_prop
    A_0:          ClassVar["gtpsa._gtpsa.ss_vect_tpsa"] = A_0
    Twiss_design: ClassVar["numpy.ndarray"] = Twiss_design
    weights:      ClassVar[list] = weights
    params:       ClassVar[list] = [prm_list, dprm_list]

transp_line = transp_line_class(prm_class)

transp_line.match()
