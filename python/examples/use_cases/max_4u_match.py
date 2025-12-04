"""Use Case:
     Module for general linear optics matching to obtain a periodic solution.
"""


import os
import sys
import copy as _copy
from dataclasses import dataclass
from typing import ClassVar
from typing import Tuple

import numpy as np
from scipy import optimize as opt

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, get_set_mpole as gs, \
    index_class as ind, linear_optics as lo, prm_class as pc, \
    courant_snyder as cs, tpsa_class as tpsa
from thor_scsi.utils.output import mat2txt, vec2txt


ind = ind.index_class()

prm_range = {
    "phi"       : [  0.0,  2.0],
    "phi_rbend" : [ -0.5,  0.0],
    "b_2"       : [-10.0, 10.0],
    "b_2_bend"  : [ -1.5,  1.5]
}

design_val = {
    "eta_des"   : [0.0, 0.0], # [eta_x, eta'_x] at the centre of the straight.
    "alpha_des" : [0.0, 0.0], # Alpha_x,y at the centre of the straight.
    "beta_des"  : [5.0, 3.0], # Beta_x,y at the centre of the straight.
    "dnu_des"   : [0.5, 0.25] # Phase advance across the straight.
}


class opt_match_class:
    # Private.

    def __init__(self, prm_class: pc.prm_class) -> None:
        self._first        = True

        self._A0           = lat_prop.ss_vect_tpsa()
        self._A_7x7        = lat_prop.ss_vect_tpsa()
        self._A            = np.nan
        self._lat_prop     = prm_class.lattice
        self._uc_centre    = prm_class.s_loc
        self._des_val_list = prm_class.design_vals
        self._weights      = prm_class.weights
        self._bend_list    = prm_class.dipoles[0]
        self._rbend_list   = prm_class.dipoles[1]
        self._prm_list     = prm_class.params[0]
        self._dprm_list    = prm_class.params[1]

        self._phi_bend     = np.zeros(len(self._bend_list))
        self._phi_rbend    = np.zeros(len(self._bend_list))

        self._phi_sp       = np.nan
        self._dphi         = np.nan

        self._b_2          = np.nan

        self._nu           = np.nan
        self._Twiss_k      = np.nan
        self._eta_centre   = np.nan
        self._xi           = np.nan

        self._constr       = {}

        self._chi_2_min    = 1e30
        self._n_iter       = -1
        self._file_name    = "opt_match.txt"

    # Public.

    def compute_periodic_cell(self, lat_prop, uc_list):
        M = lat_prop.ss_vect_tpsa()
        M.set_identity()
        # The 3rd argument is the 1st element index & the 4th the number of
        # elements to propagate across.
        lat_prop._lattice.propagate(
            lat_prop._model_state, M, uc_list[0], uc_list[1]-uc_list[0]+1)
        stable, nu, A, A_inv, _ = lo.compute_M_diag(2, M.jacobian())
        eta, Twiss = lo.transform_matrix_extract_twiss(A)
        alpha = np.array((Twiss[ind.X][0], Twiss[ind.Y][0]))
        beta = np.array((Twiss[ind.X][1], Twiss[ind.Y][1]))
        Twiss = eta, alpha, beta
        print("\ncompute_periodic_cell:")
        lat_prop.prt_Twiss_param(Twiss)
        return Twiss, A

    def compute_prop(self) -> bool:
        self._phi_sp = self._lat_prop.compute_phi_lat()
        self._dphi = self._phi_sp - self._phi_sp_0

        A1 = _copy.copy(self._A0)
        lat_prop._lattice.propagate(
            lat_prop._model_state, A1, uc_list[1]+1, sp_list[1]-uc_list[1])
        self._Twiss_k = cs.compute_Twiss_A(A1.jacobian())

        A1 = _copy.copy(self._A0)
        lat_prop._lattice.propagate(
            lat_prop._model_state, A1, uc_list[1]+1, sp_list[0]-uc_list[1])
        _, _, _, self._dnu = cs.compute_Twiss_A(A1.jacobian())
        self._Twiss_k[3][:] = 2e0*(self._Twiss_k[3][:]-self._dnu[:])

    def compute_constr(self) -> bool:
        self._constr = {
            "eta_x":   (self._Twiss_k[0][ind.x]
                        -self._des_val_list["eta_des"][ind.x])**2,
            "eta'_x":  (self._Twiss_k[0][ind.px]
                        -self._des_val_list["eta_des"][ind.px])**2,
            "alpha_x": (self._Twiss_k[1][ind.X]
                        -self._des_val_list["alpha_des"][ind.X])**2,
            "alpha_y": (self._Twiss_k[1][ind.Y]
                        -self._des_val_list["alpha_des"][ind.Y])**2,
            "beta_x":  (self._Twiss_k[2][ind.X]
                        -self._des_val_list["beta_des"][ind.X])**2,
            "beta_y":  (self._Twiss_k[2][ind.Y]
                        -self._des_val_list["beta_des"][ind.Y])**2,
            "dnu_x":   (self._Twiss_k[3][ind.X]
                        -self._des_val_list["dnu_des"][ind.X])**2,
            "dnu_y":   (self._Twiss_k[3][ind.Y]
                        -self._des_val_list["dnu_des"][ind.Y])**2
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
        print(f"    [eta_x, eta'_x]    = "
              f"[{self._Twiss_k[0][ind.x]:10.3e}"
              f", {self._Twiss_k[0][ind.px]:10.3e}] "
              f"([{self._des_val_list["eta_des"][ind.x]:10.3e}"
              f", {self._des_val_list["eta_des"][ind.px]:10.3e}])")
        print(f"    [alpha_x, alpha_y] = "
              f"[{self._Twiss_k[1][ind.x]:10.3e}"
              f", {self._Twiss_k[1][ind.px]:10.3e}] "
              f"([{self._des_val_list["alpha_des"][ind.x]:10.3e}"
              f", {self._des_val_list["alpha_des"][ind.px]:10.3e}])")
        print(f"    [beta_x, beta_y]   = "
              f"[{self._Twiss_k[2][ind.x]:5.3f}"
              f", {self._Twiss_k[2][ind.px]:5.3f}] "
              f"([{self._des_val_list["beta_des"][ind.x]:5.3f}"
              f", {self._des_val_list["beta_des"][ind.px]:5.3f}])")
        print(f"    [dnu_x, dnu_y]     = "
              f"[{self._Twiss_k[3][ind.x]:5.3f}"
              f", {self._Twiss_k[3][ind.px]:5.3f}] "
              f"([{self._des_val_list["dnu_des"][ind.x]:5.3f}"
              f", {self._des_val_list["dnu_des"][ind.px]:5.3f}])")

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

        self._prm_list.prt_prm(prm)

    def f_match(self, prm: np.ndarray) -> float:
        self._n_iter += 1
        self._prm_list.set_prm(prm)

        self.compute_prop()

        chi_2 = self.compute_chi_2()

        if chi_2 < self._chi_2_min:
            self.prt_iter(prm, chi_2)
            if False:
                self._lat_prop.prt_Twiss("twiss.txt")
            pc.prt_lat(
                self._lat_prop, self._file_name, self._prm_list,
                bend_list=self._bend_list)
            self._chi_2_min = min(chi_2, self._chi_2_min)
 
            # Problematic => system not time invariant.
            _, A = self.compute_periodic_cell(lat_prop, uc_list)
            self._A0.set_jacobian(self._A_7x7)
        else:
            print(f"\n{self._n_iter:3d} chi_2 = {self._chi_2_min:9.3e}",
                  f"({chi_2-self._chi_2_min:9.3e})")
            if False:
                print(f"\n{self._n_iter:3d}",
                      f"dchi_2 = {chi_2-self._chi_2_min:9.3e}")
                # self._prm_list.prt_prm(prm)

        return chi_2

    def opt_match(self) -> opt.OptimizeResult:
        """Use Case: optimise super period.
        """

        max_iter = 10000
        f_tol    = 1e-7
        x_tol    = 1e-7
        g_tol    = 1e-7

        Twiss_0, self._A = self.compute_periodic_cell(lat_prop, uc_list)

        print("\nmatch_straight:\n\nEntrance:")
        lat_prop.prt_Twiss_param(Twiss_0)

        self._A0.set_zero()
        self._A_7x7 = np.zeros((7, 7))
        self._A_7x7[:6, :6] = self._A
        self._A0.set_jacobian(self._A_7x7)

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


def get_bends():
    bend_lists = {
        "bend_list": [
            ["d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3", "d1_h2_sl_dm2",
             "d1_h2_sl_dm1", "d1_h2_sl_ds0", "d1_h2_sl_ds1", "d1_h2_sl_ds2",
             "d1_h2_sl_ds3", "d1_h2_sl_ds4", "d1_h2_sl_ds5",
             "d1_h2_sl_ds6"],
        ],
        "rbend_list": ["r2_h2"]
    }

    bend_list = []
    for b_list in bend_lists["bend_list"]:
        bend_list.append(pc.bend_class(lat_prop, b_list))

    return bend_list, bend_lists["rbend_list"]


def get_prms(bend_list, eps):
    prm = [
        ("q1_h2", "b_2", prm_range["b_2"]),
        ("q2_h2", "b_2", prm_range["b_2"]),
        ("r1_h2", "b_2", prm_range["b_2"]),
        ("lego", "b_2", prm_range["b_2"]),
        ("lego", "phi", [-0.5, 0.5]),
        (bend_list[0], "b_2_bend", prm_range["b_2_bend"]),
        (bend_list[0], "phi_bend", [0.5,  1.5]),
   ]

    prm_list = pc.prm_class(lat_prop, prm)

    dprm_list = np.full(len(prm), eps)

    return prm_list, dprm_list


def get_weights():
    weights = {
        "eta_x"   : 1e3,  # eta_x at the centre of the straight.
        "eta'_x"  : 1e3,  # eta'_x at the centre of the straight.
        "alpha_x" : 1e0,  # alpha_x at the centre of the straight.
        "alpha_y" : 1e0,  # alpha_y at the centre of the straight.
        "beta_x"  : 1e-3, # beta_x at the centre of the straight.
        "beta_y"  : 1e-3, # beta_y at the centre of the straight.
        "dnu_x"   : 0e-4, # dnu_x across the straight.
        "dnu_y"   : 0e-4  # dnu_y across the straight.
    }
    return weights


no = 1

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
delta_max = 3e-2
beta_inj  = np.array([3.0, 3.0])

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = lp.lattice_properties_class(file_name, E_0, cod_eps, no)

lat_prop.prt_lat("max_4u_match_lat.txt")

uc_list = np.zeros(2, dtype=int)
uc_list[0] = lat_prop._lattice.find("ucborder", 0).index
uc_list[1] = lat_prop._lattice.find("ucborder", 1).index

sp_list = np.zeros(2, dtype=int)
sp_list[0] = lat_prop._lattice.find("lsborder", 1).index
sp_list[1] = lat_prop._lattice.find("cav", 0).index

print("\nunit cell entrance           {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_list[0]].name, uc_list[0]))
print("unit cell exit               {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_list[1]].name, uc_list[1]))
print("super period last sextupole  {:5s} loc = {:d}".
      format(lat_prop._lattice[sp_list[0]].name, sp_list[0]))
print("super period exit            {:5s} loc = {:d}".
      format(lat_prop._lattice[sp_list[1]].name, sp_list[1]))

np.set_printoptions(formatter={"float": "{:10.3e}".format})

weight_list = get_weights()
bend_list, rbend_list = get_bends()
prm_list, dprm_list = get_prms(bend_list, 1e-4)

@dataclass
class prm_class:
    lattice:     ClassVar[list] = lat_prop
    s_loc:       ClassVar[list] = [uc_list, sp_list]
    design_vals: ClassVar[dict] = design_val
    weights:     ClassVar[list] = weight_list
    dipoles:     ClassVar[list] = [bend_list, rbend_list]
    params:      ClassVar[list] = [prm_list, dprm_list]

opt_match = opt_match_class(prm_class)
opt_match.opt_match()
