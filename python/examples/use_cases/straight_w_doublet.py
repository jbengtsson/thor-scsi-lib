"""Use Case:
     Optimising a straight with a matching section.
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


class opt_straight_class:
    # Private.

    def __init__(self, prm_class: pc.prm_class) -> None:
        self._first       = True

        self._lat_prop    = prm_class.lattice[0]
        self._nld         = prm_class.lattice[1]
        self._str_start   = prm_class.s_loc[0]
        self._str_centre  = prm_class.s_loc[1]
        self._str_end     = prm_class.s_loc[2]
        self._des_val     = prm_class.des_val
        self._max_iter    = prm_class.opt_prm[0]
        self._prm_tol     = prm_class.opt_prm[1]
        self._f_tol       = prm_class.opt_prm[2]
        self._eps_Jacob   = prm_class.opt_prm[3]
        self._weights     = prm_class.weights
        self._bend_list   = prm_class.dipoles[0]
        self._prm_list    = prm_class.params

        self._phi_bend    = np.zeros(len(self._bend_list))

        self._b_2         = None
        self._phi_tot_0   = None
        self._phi_tot     = None
        self._dphi        = None

        self._nu          = None
        self._length      = None
        self._alpha_c     = None
        self._eta_entr    = None
        self._alpha_entr  = None
        self._beta_entr   = None
        self._eta_centre  = None
        self._beta_centre = None
        self._xi          = None

        self._constr      = {}

        self._chi_2_min   = 1e30
        self._n_iter      = -1
        self._file_name   = "opt_straight.txt"

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
        # Compute lattice properties for the constraints.
        self._phi_tot = self._lat_prop.compute_phi_lat()
        self._dphi = self._phi_tot - self._phi_tot_0

        for k, bend in enumerate(self._bend_list):
            self._phi_bend[k] = bend.compute_bend_phi()

        lat_stable = self.compute_lat_prop()

        return lat_stable

    def compute_constr(self):
        # Compute the least-square constraints.
        # 1st parameter is length
        prm_name = self._prm_list._prm_list[0][0]
        L = self._lat_prop._lattice.find(prm_name, 0).get_length()
        self._length = 2e0*L
        self._alpha_c = self._lat_prop._alpha_c
        self._, _, _, nu_0 = self._lat_prop.get_Twiss(self._str_start)
        self._, _, _, nu_1 = self._lat_prop.get_Twiss(self._str_end)
        self._nu =  nu_1 - nu_0
        self._eta_centre, _, self._beta_centre, _ = \
            self._lat_prop.get_Twiss(self._str_centre)
        self._eta_entr, self._alpha_entr, self._beta_entr, _ = \
            self._lat_prop.get_Twiss(-1)

        self._constr = {
            "nu_x"         : (self._nu[ind.X]
                              -self._des_val['nu_des'][ind.X])**2,
            "nu_y"         : (self._nu[ind.Y]
                              -self._des_val['nu_des'][ind.Y])**2,
            "length"       : (self._length-self._des_val['length'])**2,
            "eps_x"        : (self._lat_prop._eps[ind.X]
                              -self._des_val['eps_x_des'])**2,
            "dphi"         : self._dphi**2,
            "alpha_c"      : 1e0/self._alpha_c[1]**2,
            "eta_x_entr"   : (self._eta_entr[ind.x]
                              -self._des_val['eta_x_entr_des'])**2,
            "beta_x_entr"  : (self._beta_entr[ind.X]
                              -self._des_val['beta_entr_des']
                              [ind.X])**2,
            "beta_y_entr"  : (self._beta_entr[ind.Y]
                              -self._des_val['beta_entr_des']
                              [ind.Y])**2,
            "eta_x_centre" : (self._eta_centre[ind.x]
                              -self._des_val['eta_x_centre_des'])**2,
            "beta_x_centre" : (self._beta_centre[ind.X]
                               -self._des_val['beta_centre_des']
                               [ind.X])**2,
            "beta_y_centre" : (self._beta_centre[ind.Y]
                               -self._des_val['beta_centre_des']
                               [ind.Y])**2,
            "xi"            : self._xi[ind.X]**2 + self._xi[ind.Y]**2
        }

    def compute_chi_2(self) -> float:
        self.compute_constr()

        prt_len = 13
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
        # Compute the least-square.
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
        # Print the results for the iteration.
        print(f"\n{self._n_iter:3d} chi_2 = {self._chi_2_min:9.3e}",
              f"({chi_2-self._chi_2_min:9.3e})")
        print(f"    eps_x [pm.rad] = "
              f"{1e12*self._lat_prop._eps[ind.X]:5.3f} "
              f"({1e12*self._des_val['eps_x_des']:5.3f})")

        print(f"\n    dphi [deg]     = {self._dphi:9.3e}")

        print("\n    alpha_c (multipoles zeroed)\n"
              f"                   = [{self._alpha_c[1]:9.3e}, "
              f"{self._alpha_c[2]:9.3e}]")

        print(f"\n    length         = {self._length:6.3f}                 "
              f" ({self._des_val['length']:5.3f})")
        print(f"    nu             = [{self._nu[ind.X]:7.5f}, "
              f"{self._nu[ind.Y]:7.5f}]      "
              f"([{self._des_val['nu_des'][ind.X]:7.5f}, "
              f"{self._des_val['nu_des'][ind.Y]:7.5f}])")
        print(f"    eta_x_entr     = {self._eta_entr[ind.x]:10.3e}"
              f"             "
              f" ({self._des_val['eta_x_entr_des']:9.3e})")
        print(f"    beta_entr      = "
              f"[{self._beta_entr[ind.X]:5.3f}"
              f", {self._beta_entr[ind.Y]:5.3f}]          "
              f"([{self._des_val['beta_entr_des'][ind.X]:5.3f}, "
              f"{self._des_val['beta_entr_des'][ind.Y]:5.3f}])")
        print(f"    eta_x_centre   = {self._eta_centre[ind.x]:10.3e}"
              f"              "
              f"({self._des_val['eta_x_centre_des']:9.3e})")
        print(f"    beta_centre    = "
              f"[{self._beta_centre[ind.X]:6.3f}"
              f", {self._beta_centre[ind.Y]:6.3f}]        "
              f"([{self._des_val['beta_centre_des'][ind.X]:6.3f}, "
              f"{self._des_val['beta_centre_des'][ind.Y]:6.3f}])")
        print(f"    xi             = [{self._xi[ind.X]:5.3f}, "
              f"{self._xi[ind.Y]:5.3f}]")

        print(f"\n    phi_tot        = {self._phi_tot:8.5f}")
        print(f"    C [m]          = "
              f"{self._lat_prop.compute_circ():8.5f}")

        print()
        for k, phi in enumerate(self._phi_bend):
                self._b_2 = \
                    self._bend_list[k].compute_bend_b_2xL() \
                    /self._bend_list[k].compute_bend_L_tot()
                print(f"    phi_bend_{k+1:1d}     = "
                      f"[{phi:8.5f}, {self._b_2:8.5f}]")

        print(f"\n    b_3            =", self._nld._b_3_list)

        self._lat_prop.prt_rad()
        self._prm_list.prt_prm(prm)

    def f_straight(self, prm: np.ndarray) -> float:
        # Least-square function.
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

    def opt_straight(self) -> opt.OptimizeResult:
        # Optimiser.

        prm, bounds = self._prm_list.get_prm()
        # Initial total bend angle.
        self._phi_tot_0 = self._lat_prop.compute_phi_lat()
        self.f_straight(prm)

        # Methods:
        #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
        #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.
        # Methods with boundaries:
        #   L-BFGS-B, TNC, and SLSQP.
        # Methods with constraints:
        #   SLSQP.

        opt_dict = {
            # Without boundaries.
            "CG": {"options": {
                "maxiter": self._max_iter, "gtol": self._f_tol,
                "eps": self._eps_Jacob}},
            "BFGS": {"options": {
                "maxiter": self._max_iter,
                "gtol": self._f_tol, "xrtol": self._prm_tol,
                "eps": self._prm_tol}},
            # With boundaries.
            "L-BFGS-B": {"options": {
                "maxiter": self._max_iter,
                "ftol": self._f_tol, "gtol": self._f_tol,
                "eps": self._prm_tol}},
            "SLSQP": {"options": {
                "ftol": self._f_tol, "maxiter": self._max_iter,
                "eps": self._eps_Jacob}},
            "TNC": {"options": {
                "maxfun": self._max_iter, "ftol": self._f_tol,
                "gtol": self._f_tol, "xtol": self._prm_tol,
                "eps": self._eps_Jacob}}
        }

        method = "TNC"
        
        minimum = opt.minimize(
            self.f_straight, prm, method=method,
            # callback=prt_iter,
            jac=None, bounds = bounds, options=opt_dict[method]["options"]
        )

        print("\n".join(minimum))
        # print(minimum)
        # for key, value in minimum.items():
        #     print(f"{key}: {value}")


def init():
    # TPSA order - i.e., include linear chromaticity.
    no = 2

    cod_eps = 1e-15
    E_0     = 2.5e9

    home_dir = os.path.join(
        os.environ["HOME"], "Nextcloud", "thor_scsi", "JB")
    file_name = os.path.join(home_dir, sys.argv[1]+".lat")

    lat_prop = lp.lattice_properties_class(file_name, E_0, cod_eps, no)

    lat_prop.prt_lat("Twiss_lat.txt")

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

    return lat_prop


def define_system(lat_prop):
    A_max     = np.array([6e-3, 3e-3])
    delta_max = 3e-2
    beta_inj  = np.array([3.0, 3.0])

    str_start = lat_prop._lattice.find("QF2", 0).index-1
    str_centre = lat_prop._lattice.find("D5", 0).index
    str_end = lat_prop._lattice.find("QF2", 1).index
    print("Straight start {:5s} loc = {:d}".
          format(lat_prop._lattice[str_start].name, str_start))
    print("Straight centre {:5s} loc = {:d}".
          format(lat_prop._lattice[str_centre].name, str_centre))
    print("Straigth end {:5s} loc = {:d}".
          format(lat_prop._lattice[str_end].name, str_end))

    b_3_list = []
    nld = nld_class.nonlin_dyn_class(
        lat_prop, A_max, beta_inj, delta_max, b_3_list)

    # Parameter ranges.
    prms_range = {
        "phi":       [ -5.0,  5.0],
        "b_2":       [-10.0, 10.0],
        "b_2_bend":  [-10.0, 10.0]
    }

    # Design values.
    design_values = {
        "nu_des"           : [0.9, 0.3],
        "length"           : 5.0,
        "eps_x_des"        : 150e-12,
        "eta_x_entr_des"   : 3.88155e-02,
        "beta_entr_des"    : [4.26217, 2.43390],
        "eta_x_centre_des" : 0.0,
        "beta_centre_des"  : [10.0, 2.0]
    }

    # Varying bend radius dipole.
    bend_lists = [["B2_1", "B2_2", "B2_3", "B2_4", "B2_5"]]

    bend_list = []
    bend_list.append(pc.bend_class(lat_prop, bend_lists[0]))

    # Parameters.
    if False:
        prms = [
            ("D5",         "L",        [0.4, 3.0]),
            (bend_list[0], "phi_bend", prms_range["phi"]),
            (bend_list[0], "b_2_bend", prms_range["b_2_bend"]),
            ("QF2",        "b_2",      prms_range["b_2"]),
            ("QF3",        "b_2",      prms_range["b_2"]),
            ("QD1",        "b_2",      prms_range["b_2"])
        ]
    else:
        # Straight should be 1st parameter.
        prms = [
            ("D5",   "L",   [0.4,  3.0]),
            # ("D2",   "L",   [0.10, 0.25]),
            # ("D3",   "L",   [0.10, 0.25]),
            # ("D4",   "L",   [0.10, 0.25]),

            ("B2_1", "phi", prms_range["phi"]),
            ("B2_2", "phi", prms_range["phi"]),
            ("B2_3", "phi", prms_range["phi"]),
            ("B2_4", "phi", prms_range["phi"]),
            ("B2_5", "phi", prms_range["phi"]),

            ("B2_1", "b_2", prms_range["b_2_bend"]),
            ("B2_2", "b_2", prms_range["b_2_bend"]),
            ("B2_3", "b_2", prms_range["b_2_bend"]),
            ("B2_4", "b_2", prms_range["b_2_bend"]),
            ("B2_5", "b_2", prms_range["b_2_bend"]),

            ("QF2",  "b_2", prms_range["b_2"]),
            ("QD1",  "b_2", prms_range["b_2"]),
            ("QF3",  "b_2", prms_range["b_2"])

            # ("QF2",  "L",   [0.2, 0.3]),
            # ("QF3",  "L",   [0.2, 0.3]),
            # ("QD1",  "L",   [0.2, 0.3])
        ]

    # Weights for least-square minimisation.
    weight_list = {
        "length"        : 1e-3,
        "eps_x"         : 1e14,
        "dphi"          : 0e0,
        "alpha_c"       : 1e-13,
        "nu_x"          : 0e-3,
        "nu_y"          : 0e-3,
        "eta_x_entr"    : 1e3,
        "beta_x_entr"   : 1e-3,
        "beta_y_entr"   : 1e-3,
        "eta_x_centre"  : 1e1,
        "beta_x_centre" : 1e-7,
        "beta_y_centre" : 1e-7,
        "xi"            : 1e-7
    }

    # Epsilon for numerical evalutation of the Jacobian. 
    eps_Jacob = 1e-4

    prm_list = pc.prm_class(lat_prop, prms)

   # Package the system constrants & parameters.
    @dataclass
    class prm_class:
        # Max number of iterations, parameter Precision, and optimal function
        # precision.
        opt_prm: ClassVar[list] = [10000, 1e-5, 1e-7, eps_Jacob]
        lattice: ClassVar[list] = [lat_prop, nld]
        s_loc:   ClassVar[list] = [str_start, str_centre, str_end]
        des_val: ClassVar[dict] = design_values
        weights: ClassVar[list] = weight_list
        dipoles: ClassVar[list] = [bend_list]
        params:  ClassVar[list] = prm_list

    # Generate the corresponding object.
    opt_straight = opt_straight_class(prm_class)

    return opt_straight


# Main program.
lat_prop = init()
opt_straight = define_system(lat_prop)

# Optimise.
opt_straight.opt_straight()
