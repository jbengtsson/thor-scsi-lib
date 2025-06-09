"""Use Case:
     Implementation and general optimisation of a higher-order-achromat.
     Use one unit cell.
"""


import os
import sys
from dataclasses import dataclass
from typing import ClassVar

import math
import numpy as np
from scipy import optimize as opt
from scipy.interpolate import CubicSpline

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, prm_class as pc, nonlin_dyn as nld_class
from thor_scsi.utils.output import mat2txt, vec2txt

import gtpsa


ind = ind.index_class()

b_2_bend_max = 1.0
b_2_max      = 10.0

eps_x_des    = 60e-12
phi_1_des    =  1.4
phi_rb_1_des = -0.2
b_2_des      =  2.0
alpha_c_des  = 0.5e-4
nu_uc_des    = [0.4, 0.1]
nu_sp_des    = np.array([57.2/20.0,  20.75/20.0])
# nu_sp_des    = np.array([57.2/20.0,  17.75/20.0])
# nu_sp_des    = np.array([57.47/20.0, 17.68/20.0])
if False:
    # With 4 unit dipole cells.
    nu_sp_des -= nu_uc_des
beta_des     = [5.0, 3.0]
dnu_des      = [0.75, 0.25]     # Phase advance across the straight.


@dataclass
class gtpsa_prop:
    # GTPSA properties.
    # Number of variables - phase-space coordinates & 1 for parameter
    #dependence.
    nv: ClassVar[int] = 6 + 1
    # Max order.
    no: ClassVar[int] = 1
    # Truncation order.
    to: ClassVar[int] = 1
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]

    def tpsa():
        return gtpsa.tpsa(gtpsa_prop.desc, gtpsa_prop.no)
    def ctpsa():
        return gtpsa.ctpsa(gtpsa_prop.desc, gtpsa_prop.no)
    def ss_vect_tpsa():
        return gtpsa.ss_vect_tpsa(gtpsa_prop.desc, gtpsa_prop.no)


def get_phi(el):
    if (type(el) == ts.Bending) or (type(el) == ts.Quadrupole):
        if el.get_curvature() is not None:
            L = el.get_length()
            phi = math.degrees(L*el.get_curvature())
        else:
            phi = 0e0
    else:
        phi = 0e0
    return phi


def compute_alpha_c(map):
    C = lat_prop.compute_circ()
    index = np.zeros(gtpsa_prop.nv, int)
    alpha_c = np.zeros(gtpsa_prop.no+1)
    for k in range(1, gtpsa_prop.no+1):
        alpha_c[k] = map.ct.get([0, 0, 0, 0, k])/C
    return alpha_c


class opt_sp_class:
    # Private

    def __init__(
            self, lat_prop, prm_list, dprm_list, uc_list, sp_list, weight,
            bend_list, rbend_list, eps_x_des, phi_1_des, phi_rb_1_des, b_2_des,
            nu_uc_des, nu_sp_des, beta_des, dnu_des, nld):

        self._lat_prop       = lat_prop
        self._nld            = nld
        self._prm_list       = prm_list
        self._uc_list        = uc_list
        self._sp_list        = sp_list
        self._weight         = weight
        self._bend_list      = bend_list
        self._rbend_list     = rbend_list
        self._eps_x_des      = eps_x_des
        self._phi_1_des      = phi_1_des
        self._phi_rb_1_des   = phi_rb_1_des
        self._b_2_des        = b_2_des
        self._nu_uc_des      = nu_uc_des
        self._nu_sp_des      = nu_sp_des
        self._beta_des       = beta_des
        self._dnu_des        = dnu_des

        self._phi_bend       = np.zeros(3)
        self._phi_rbend      = np.zeros(3)
        self._b_2            = np.nan
        self._phi_sp_0       = np.nan
        self._phi_sp         = np.nan
        self._dphi           = np.nan
        self._Twiss_sp       = np.nan
        self._eta_list       = np.zeros((2, 4))
        self._alpha_list     = np.zeros((2, 2))
        self._nu_uc          = np.nan
        self._alpha_c        = np.nan
        self._nu_0           = np.nan
        self._dnu            = np.nan
        self._xi             = np.nan
        self._eta_nl         = gtpsa_prop.tpsa()

        self._chi_2_min      = 1e30
        self._n_iter         = -1
        self._file_name      = "opt_sp.txt"

        self._h_im_scl_rms_0 = np.nan
        self._h_im_scl_rms_1 = np.nan
        self._K_re_scl_rms_0 = np.nan
        self._K_re_scl_rms_1 = np.nan

    # Public.

    def prt_iter(self, prm, chi_2):
        eta, alpha, beta, nu_sp = self._Twiss_sp

        print(f"\n{self._n_iter:3d} chi_2 = {self._chi_2_min:9.3e}",
              f"({chi_2-self._chi_2_min:9.3e})")
        print(f"    eps_x [pm.rad] = {1e12*self._lat_prop._eps[ind.X]:5.3f}",
              f"({1e12*self._eps_x_des:5.3f})")

        print(f"\n    dphi [deg]     = {self._dphi:9.3e}")

        print(f"\n    alpha_c        = [{self._alpha_c[1]:9.3e},",
              f"{self._alpha_c[2]:9.3e}]")
        print(f"    nu_uc          = [{self._nu_uc[ind.X]:7.5f},",
              f"{self._nu_uc[ind.Y]:7.5f}] ([{self._nu_uc_des[ind.X]:7.5f},",
              f"{self._nu_uc_des[ind.Y]:7.5f}])")
        print(f"    nu_sp          = [{nu_sp[ind.X]:7.5f},",
              f"{nu_sp[ind.Y]:7.5f}] ([{self._nu_sp_des[ind.X]:7.5f},",
              f"{self._nu_sp_des[ind.Y]:7.5f}])")
        print(f"    dnu            = [{self._dnu[ind.X]:7.5f},",
              f"{self._dnu[ind.Y]:7.5f}] ([{self._dnu_des[ind.X]:7.5f},",
              f"{self._dnu_des[ind.Y]:7.5f}])")
        print(f"    xi             = [{self._xi[ind.X]:5.3f},",
              f"{self._xi[ind.Y]:5.3f}]")
        print(f"    eta^(2)_x      = {self._eta_nl.get([0, 0, 0, 0, 2]):9.3e}")

        print(f"\n    eta'_uc        = {self._eta_list[0][ind.px]:9.3e}",
              f" {self._eta_list[1][ind.px]:9.3e}")
        print(f"    alpha_uc       = [{self._alpha_list[0][ind.X]:9.3e},",
              f"{self._alpha_list[0][ind.Y]:9.3e}]",
              f"[{self._alpha_list[1][ind.X]:9.3e},",
              f"{self._alpha_list[1][ind.Y]:9.3e}]")
        print(f"    eta_x          = {eta[ind.x]:9.3e}")
        print(f"    beta           = [{beta[ind.X]:7.5f}, {beta[ind.Y]:7.5f}]",
              f"([{self._beta_des[ind.X]:7.5f}, {self._beta_des[ind.Y]:7.5f}])")

        print(f"\n    phi_sp         = {self._phi_sp:8.5f}")
        print(f"    C [m]          = {self._lat_prop.compute_circ():8.5f}")

        print()
        for phi in self._phi_bend:
            print(f"    phi_bend       = {phi:8.5f}")
        print()
        for k, rbend in enumerate(self._rbend_list):
            print(f"    {rbend:10s}     = {self._phi_rbend[k]:8.5f}")

        print(f"\n    b_2            = {self._b_2:8.5f}")

        self._lat_prop.prt_rad()
        self._prm_list.prt_prm(prm)

    def compute_chi_2(self):
        prt = not False

        eta, alpha, beta, nu_sp = self._Twiss_sp

        dchi_2 = weight[0]*(self._lat_prop._eps[ind.X]-self._eps_x_des)**2
        chi_2 = dchi_2
        if prt:
            print(f"\n  dchi2(eps_x)        = {dchi_2:9.3e}")

        dchi_2 = weight[1]*(self._phi_bend[0]-self._phi_1_des)**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(phi_1)        = {dchi_2:9.3e}")

        dchi_2 = weight[2]*(self._phi_rbend[0]-self._phi_rb_1_des)**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(phi_rb_1)     = {dchi_2:9.3e}")

        dchi_2 = weight[3]*(self._b_2-self._b_2_des)**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(b_2)          = {dchi_2:9.3e}")

        dchi_2 = weight[4]*(self._dphi)**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(dphi)         = {dchi_2:9.3e}")

        dchi_2 = weight[5]*1e0/self._alpha_c[1]**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(alpha^(1)_c)  = {dchi_2:9.3e}")

        dchi_2 = weight[6]*(self._alpha_c[2])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(alpha^(2)_c)  = {dchi_2:9.3e}")

        dchi_2 = weight[7]*self._lat_prop._U_0**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(U_0)          = {dchi_2:9.3e}")

        dchi_2 = weight[8]*(self._eta_list[0][ind.px]**2
                            +self._eta_list[1][ind.px]**2)
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(eta'_uc)      = {dchi_2:9.3e}")

        dchi_2 = weight[9]*(
            self._alpha_list[0][ind.X]**2+self._alpha_list[0][ind.Y]**2
            +self._alpha_list[1][ind.X]**2+self._alpha_list[1][ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(alpha_uc)     = {dchi_2:9.3e}")

        dchi_2 = weight[10]*(self._nu_uc[ind.X]-self._nu_uc_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(nu_uc_x)      = {dchi_2:9.3e}")

        dchi_2 = weight[11]*(self._nu_uc[ind.Y]-self._nu_uc_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(nu_uc_y)      = {dchi_2:9.3e}")

        dchi_2 = weight[12]*eta[ind.x]**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(eta_x)        = {dchi_2:9.3e}")

        dchi_2 = weight[13]*(nu_sp[ind.X]-self._nu_sp_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(nu_sp_x)      = {dchi_2:9.3e}")

        dchi_2 = weight[14]*(nu_sp[ind.Y]-self._nu_sp_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(nu_sp_y)      = {dchi_2:9.3e}")

        dchi_2 = weight[15]*(beta[ind.X]-self._beta_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(beta_x)       = {dchi_2:9.3e}")

        dchi_2 = weight[16]*(beta[ind.Y]-self._beta_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(beta_y)       = {dchi_2:9.3e}")

        dchi_2 = weight[17]*(self._dnu[ind.X]-self._dnu_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(dnu_x)        = {dchi_2:9.3e}")

        dchi_2 = weight[18]*(self._dnu[ind.Y]-self._dnu_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(dnu_y)        = {dchi_2:9.3e}")

        dchi_2 = weight[19]*(self._xi[ind.X]**2+self._xi[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(xi)           = {dchi_2:9.3e}")

        dchi_2 = weight[20]*self._eta_nl.get([0, 0, 0, 0, 2])**2
        chi_2 += dchi_2
        if prt:
            print(f"  dchi2(eta^(2)_x)    = {dchi_2:9.3e}")

        return chi_2

    def opt_sp(self):
        """Use Case: optimise super period.
        """

        def f_sp(prm):
            self._n_iter += 1
            self._prm_list.set_prm(prm)

            self._phi_sp = self._lat_prop.compute_phi_lat()
            self._dphi = self._phi_sp - self._phi_sp_0

            for k, bend in enumerate(self._bend_list):
                self._phi_bend[k] = bend.compute_bend_phi()

            for k, bend in enumerate(self._rbend_list):
                self._phi_rbend[k] = self._lat_prop.get_phi_elem(bend, 0)

            self._b_2 = \
                self._bend_list[0].compute_bend_b_2xL() \
                /self._bend_list[0].compute_bend_L_tot()

            chi_2 = 0e0
            # Compute Twiss parameters along the lattice.
            stable = self._lat_prop.comp_per_sol(gtpsa_prop)
            if not stable:
                print("\ncomp_per_sol: unstable")
                chi_2 = 1e30
            # Compute radiation properties.
            stable, stable_rad = self._lat_prop.compute_radiation()
            if not stable or not stable_rad:
                print("\ncompute_radiation: unstable")
                chi_2 = 1e30
            # Compute linear chromaticity.
            self._nld.compute_map(gtpsa_prop, lat_prop)
            stable, _, self._xi = lo.compute_nu_xi(
                gtpsa_prop, self._nld._M)
            if not stable:
                print("\ncompute_nu_xi: unstable")
                chi_2 = 1e30
             # Compute 2nd order dispersion
            self._eta_nl = nld.compute_eta(
                gtpsa_prop, lat_prop, lat_prop._M)

            if chi_2 != 1e30:
                self._alpha_c = compute_alpha_c(self._nld._M)
                self._eta_list[0], self._alpha_list[0], _, self._nu_0 = \
                    self._lat_prop.get_Twiss(uc_list[0])
                self._eta_list[1], self._alpha_list[1], _, self._nu_1 = \
                    self._lat_prop.get_Twiss(uc_list[1])
                self._nu_uc = self._nu_1 - self._nu_0
                if True or (len(uc_list) == 3):
                    self._eta_list[1], self._alpha_list[1], _, _ = \
                        self._lat_prop.get_Twiss(self._uc_list[1])
                self._Twiss_sp = self._lat_prop.get_Twiss(-1)
                self._dnu = 2e0*self._lat_prop.get_Twiss(sp_list[0])[3][:]

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
            else:
                if not False:
                    print(f"\n{self._n_iter:3d} chi_2 = {chi_2:11.5e}",
                          f"({self._chi_2_min:11.5e})")
                    self._prm_list.prt_prm(prm)

            return chi_2

        max_iter = 10000
        f_tol    = 1e-7
        x_tol    = 1e-7
        g_tol    = 1e-7

        prm, bounds = self._prm_list.get_prm()
        if True:
            self._phi_sp_0 = 18.0
        else:
            self._phi_sp_0 = lat_prop.compute_phi_lat()
        f_sp(prm)

        # Methods:
        #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
        #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.
        # Methods with boundaries:
        #   L-BFGS-B, TNC, and SLSQP.
        # Methods with constraints:
        #   SLSQP.

        optimiser = "CG"
        if optimiser == "TNC":
            minimum = opt.minimize(
                f_sp,
                prm,
                method="TNC",
                # callback=prt_iter,
                bounds = bounds,
                options={"ftol": f_tol, "gtol": g_tol, "xtol": x_tol,
                         "maxfun": max_iter, "eps": 1e-3}
            )
        if optimiser == "BFGS":
            minimum = opt.minimize(
                f_sp,
                prm,
                method="BFGS",
                # callback=prt_iter,
                bounds = bounds,
                options={"ftol": f_tol, "gtol": g_tol, "maxiter": max_iter,
                         "eps": 1e-5}
            )
        if optimiser == "SLSQP":
            # eps - step size for numerical evalution of the Jacobian.
            minimum = opt.minimize(
                f_sp,
                prm,
                method="SLSQP",
                # callback=prt_iter,
                bounds = bounds,
                options={"ftol": f_tol, "maxiter": max_iter, "eps": 1e-4}
            )
        elif optimiser == "CG":
            # eps - step size.

            minimum = opt.minimize(
                f_sp,
                prm,
                method="CG",
                # callback=prt_iter,
                jac=None,
                options={"gtol": g_tol, "maxiter": max_iter, "eps": dprm_list}
            )

        print("\n".join(minimum))


def get_bends(lat):
    if lat == 1:
        # max_4u/m4U_250316_h03_01_01_01_tracy-2.

        d1_list = [
            "d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3", "d1_h2_sl_dm2",
            "d1_h2_sl_dm1",
            "d1_h2_sl_d0", "d1_h2_sl_ds1", "d1_h2_sl_ds2", "d1_h2_sl_ds3",
            "d1_h2_sl_ds4", "d1_h2_sl_ds5"
        ]
        d2_list = [
            "d2_h2_sl_df0", "d2_h2_sl_df1", "d2_h2_sl_df2", "d2_h2_sl_df3",
            "d2_h2_sl_df4", "d2_h2_sl_df5", "d2_h2_sl_df6"
        ]
        d3_list = [
            "d3_h2_sl_df0", "d3_h2_sl_df1", "d3_h2_sl_df2", "d3_h2_sl_df3",
            "d3_h2_sl_df4", "d3_h2_sl_df5", "d3_h2_sl_df6"
        ]

        d1_bend = pc.bend_class(lat_prop, d1_list)
        d2_bend = pc.bend_class(lat_prop, d2_list)
        d3_bend = pc.bend_class(lat_prop, d3_list)
        bend_list = [d1_bend, d2_bend, d3_bend]

        rbend_list = ["r1_h2", "r2_h2", "r3_h2"]

        b_3_list = ["s3_h2", "s4_h2"]
    elif lat == 2:
        # m4U_241223_h02_01_01_01_tracy_2.

        d1_list = [
            "d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3", "d1_h2_sl_dm2",
            "d1_h2_sl_dm1",
            "d1_h2_sl_ds0", "d1_h2_sl_ds1", "d1_h2_sl_ds2", "d1_h2_sl_ds3",
            "d1_h2_sl_ds4", "d1_h2_sl_ds5", "d1_h2_sl_ds6"
        ]
        d2_list = [
            "d2_h2_sl_d0a", "d2_h2_sl_d0b", "d2_h2_sl_d0c", "d2_h2_sl_df1",
            "d2_h2_sl_df2", "d2_h2_sl_df3", "d2_h2_sl_df4", "d2_h2_sl_df5"
        ]

        d1_bend = pc.bend_class(lat_prop, d1_list)
        d2_bend = pc.bend_class(lat_prop, d2_list)
        bend_list = [d1_bend, d2_bend]

        rbend_list = ["r1_h2", "r2_h2"]
        b_3_list = ["s3_h2", "s4_h2"]

    return bend_list, rbend_list, b_3_list


def get_prms(set, bend_list, eps):
    if set == 1:
        prm = [
            ("q1_h2",      "b_2",     -10.0, 10.0),
            ("q2_h2",      "b_2",     -10.0, 10.0),

            ("r1_h2",      "b_2",     -10.0, 10.0),
            ("r2_h2",      "b_2",     -10.0, 10.0),
            ("r3_h2",      "b_2",     -10.0, 10.0),

            (bend_list[0], "b_2_bend", -1.5,  1.5),
            (bend_list[1], "b_2_bend", -1.5,  1.5),
            (bend_list[2], "b_2_bend", -1.5,  1.5),

            ("r1_h2",      "phi",      -0.5,  0.5),
            ("r2_h2",      "phi",      -0.2,  0.2),
            ("r3_h2",      "phi",      -0.2,  0.2),

            (bend_list[0], "phi_bend",  1.4,  1.5),
            (bend_list[1], "phi_bend",  1.5,  3.0),
            (bend_list[2], "phi_bend",  1.5,  2.5)
        ]
    elif set == 2:
        prm = [
            ("q1_h2",        "b_2",      -10.0, 10.0),
            ("q2_h2",        "b_2",      -10.0, 10.0),

            ("r1_h2",        "b_2",      -10.0, 10.0),
            ("r2_h2",        "b_2",      -10.0, 10.0),
            ("r3_h2",        "b_2",      -10.0, 10.0),

            (bend_list[0],   "b_2_bend", -1.5,  1.5),
            (bend_list[1],   "b_2_bend", -1.5,  1.5),
            (bend_list[2],   "b_2_bend", -1.5,  1.5),

            ("r1_h2",        "phi",      -0.2,  0.2),
            ("r2_h2",        "phi",      -0.2,  0.0),
            ("r3_h2",        "phi",      -0.2,  0.0),

            ("d1_h2_sl_dm5", "phi",      -1.5,  1.5),
            ("d1_h2_sl_dm4", "phi",      -1.5,  1.5),
            ("d1_h2_sl_dm3", "phi",      -1.5,  1.5),
            ("d1_h2_sl_dm2", "phi",      -1.5,  1.5),
            ("d1_h2_sl_dm1", "phi",      -1.5,  1.5),
            ("d1_h2_sl_d0",  "phi",      -1.5,  1.5),
            ("d1_h2_sl_ds1", "phi",      -1.5,  1.5),
            ("d1_h2_sl_ds2", "phi",      -1.5,  1.5),
            ("d1_h2_sl_ds3", "phi",      -1.5,  1.5),
            ("d1_h2_sl_ds4", "phi",      -1.5,  1.5),
            ("d1_h2_sl_ds5", "phi",      -1.5,  1.5),

            ("d2_h2_sl_df0", "phi",      -1.5,  1.5),
            ("d2_h2_sl_df1", "phi",      -1.5,  1.5),
            ("d2_h2_sl_df2", "phi",      -1.5,  1.5),
            ("d2_h2_sl_df3", "phi",      -1.5,  1.5),
            ("d2_h2_sl_df4", "phi",      -1.5,  1.5),
            ("d2_h2_sl_df5", "phi",      -1.5,  1.5),
            ("d2_h2_sl_df6", "phi",      -1.5,  1.5),

            ("d3_h2_sl_df0", "phi",      -1.5,  1.5),
            ("d3_h2_sl_df1", "phi",      -1.5,  1.5),
            ("d3_h2_sl_df2", "phi",      -1.5,  1.5),
            ("d3_h2_sl_df3", "phi",      -1.5,  1.5),
            ("d3_h2_sl_df4", "phi",      -1.5,  1.5),
            ("d3_h2_sl_df5", "phi",      -1.5,  1.5),
            ("d3_h2_sl_df6", "phi",      -1.5,  1.5)
        ]
    elif set == 3:
        prm = [
            ("q1_h2",        "b_2",    -10.0, 10.0),
            ("q2_h2",        "b_2",    -10.0, 10.0),

            ("r1_h2",        "b_2",    -10.0, 10.0),
            ("r2_h2",        "b_2",    -10.0, 10.0),
            ("r3_h2",        "b_2",    -10.0, 10.0),

            ("d1_h2_sl_dm5", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_dm4", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_dm3", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_dm2", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_dm1", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_d0",  "b_2",    -10.0, 10.0),
            ("d1_h2_sl_ds1", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_ds2", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_ds3", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_ds4", "b_2",    -10.0, 10.0),
            ("d1_h2_sl_ds5", "b_2",    -10.0, 10.0),

            ("d2_h2_sl_df0", "b_2",    -10.0, 10.0),
            ("d2_h2_sl_df1", "b_2",    -10.0, 10.0),
            ("d2_h2_sl_df2", "b_2",    -10.0, 10.0),
            ("d2_h2_sl_df3", "b_2",    -10.0, 10.0),
            ("d2_h2_sl_df4", "b_2",    -10.0, 10.0),
            ("d2_h2_sl_df5", "b_2",    -10.0, 10.0),
            ("d2_h2_sl_df6", "b_2",    -10.0, 10.0),

            ("d3_h2_sl_df0", "b_2",    -10.0, 10.0),
            ("d3_h2_sl_df1", "b_2",    -10.0, 10.0),
            ("d3_h2_sl_df2", "b_2",    -10.0, 10.0),
            ("d3_h2_sl_df3", "b_2",    -10.0, 10.0),
            ("d3_h2_sl_df4", "b_2",    -10.0, 10.0),
            ("d3_h2_sl_df5", "b_2",    -10.0, 10.0),
            ("d3_h2_sl_df6", "b_2",    -10.0, 10.0),

            ("r1_h2",        "phi",     -0.2,  0.2),
            ("r2_h2",        "phi",     -0.2,  0.0),
            ("r3_h2",        "phi",     -0.2,  0.0),

            ("d1_h2_sl_dm5", "phi",     -1.5,  1.5),
            ("d1_h2_sl_dm4", "phi",     -1.5,  1.5),
            ("d1_h2_sl_dm3", "phi",     -1.5,  1.5),
            ("d1_h2_sl_dm2", "phi",     -1.5,  1.5),
            ("d1_h2_sl_dm1", "phi",     -1.5,  1.5),
            ("d1_h2_sl_d0",  "phi",     -1.5,  1.5),
            ("d1_h2_sl_ds1", "phi",     -1.5,  1.5),
            ("d1_h2_sl_ds2", "phi",     -1.5,  1.5),
            ("d1_h2_sl_ds3", "phi",     -1.5,  1.5),
            ("d1_h2_sl_ds4", "phi",     -1.5,  1.5),
            ("d1_h2_sl_ds5", "phi",     -1.5,  1.5),

            ("d2_h2_sl_df0", "phi",     -1.5,  1.5),
            ("d2_h2_sl_df1", "phi",     -1.5,  1.5),
            ("d2_h2_sl_df2", "phi",     -1.5,  1.5),
            ("d2_h2_sl_df3", "phi",     -1.5,  1.5),
            ("d2_h2_sl_df4", "phi",     -1.5,  1.5),
            ("d2_h2_sl_df5", "phi",     -1.5,  1.5),
            ("d2_h2_sl_df6", "phi",     -1.5,  1.5),

            ("d3_h2_sl_df0", "phi",     -1.5,  1.5),
            ("d3_h2_sl_df1", "phi",     -1.5,  1.5),
            ("d3_h2_sl_df2", "phi",     -1.5,  1.5),
            ("d3_h2_sl_df3", "phi",     -1.5,  1.5),
            ("d3_h2_sl_df4", "phi",     -1.5,  1.5),
            ("d3_h2_sl_df5", "phi",     -1.5,  1.5),
            ("d3_h2_sl_df6", "phi",     -1.5,  1.5)
        ]

    prm_list = pc.prm_class(lat_prop, prm)

    dprm_list = np.full(len(prm), eps)

    return prm_list, dprm_list


def get_weights():
    weight = np.array([
        1e17,  # 0,  eps_x.
        1e-1,  # 1,  phi_1.
        0e-6,  # 2,  phi_rb_1.
        1e-6,  # 3,  b_2.
        1e-1,  # 4,  dphi.
        0e-14, # 5,  alpha^(1)_c.
        0e1,   # 6,  alpha^(2)_c.
        1e-15, # 7,  U_0.
        1e2,   # 8,  etap_x_uc.
        1e-2,  # 9,  alpha_uc.
        1e-2,  # 10, nu_uc_x.
        1e-2,  # 11, nu_uc_y.
        1e2,   # 12, eta_x.
        1e0,   # 13, nu_sp_x.
        1e0,   # 14, nu_sp_y.
        0e-6,  # 15, beta_x.
        0e-6,  # 16, beta_y.
        1e-3,  # 17, dnu_x.
        1e-3,  # 18, dnu_y.
        1e-7,  # 19, xi.
        1e-6   # 20, eta^(2)_x.
    ])

    return weight


# TPSA max order.
gtpsa_prop.no = 2
gtpsa_prop.to = gtpsa_prop.no
gtpsa_prop.desc = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)

cod_eps = 1e-10
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
delta_max = 3e-2
beta_inj  = np.array([3.0, 3.0])

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")

file_name = os.path.join(home_dir, sys.argv[1]+".lat")

lat_prop = lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)
lat_prop.prt_lat("lat_prop_lat.txt")

if False:
    # Check.
    d1_list = [
        "d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3", "d1_h2_sl_dm2",
        "d1_h2_sl_dm1",
        "d1_h2_sl_d0", "d1_h2_sl_ds1", "d1_h2_sl_ds2", "d1_h2_sl_ds3",
        "d1_h2_sl_ds4", "d1_h2_sl_ds5"
    ]

    d1_bend = pc.bend_class(lat_prop, d1_list)

    prm = [(d1_bend, "b_2_bend", -1.5,  1.5)]

    prm_list = pc.prm_class(lat_prop, prm)

    prm_list._prm_list[0][0].print()
    prm_list.set_prm([2.0*prm_list.get_prm()[0][0]])
    print()
    prm_list._prm_list[0][0].print()

    assert False, "\nCheck completed"

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol(gtpsa_prop):
        print("\ncomp_per_sol: unstable")
        raise ValueError

    # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation: unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop.prt_lat_param(gtpsa_prop)
    lat_prop.prt_rad()
    lat_prop.prt_M()
    lat_prop.prt_M_rad()

uc_list = []
uc_list.append(lat_prop._lattice.find("d2_h2_sl_df0", 0).index)
uc_list.append(lat_prop._lattice.find("d3_h2_sl_df0", 1).index)
uc_list = np.array(uc_list)

sp_list = []
sp_list.append(lat_prop._lattice.find("lsborder", 0).index)
sp_list.append(lat_prop._lattice.find("lsborder", 1).index)
sp_list = np.array(sp_list)

print(f"\nunit cell entrance           ",
      f"{lat_prop._lattice[uc_list[0]].name:15s} loc = {uc_list[0]:d}")
print(f"unit cell exit               {lat_prop._lattice[uc_list[1]].name:15s}",
      f"loc = {uc_list[1]:d}")
print(f"super period first sextupole {lat_prop._lattice[sp_list[0]].name:15s}",
      f"loc = {sp_list[0]:d}")
print(f"super period last sextupole  {lat_prop._lattice[sp_list[1]].name:15s}",
      f" loc = {sp_list[1]:d}")

bend_list, rbend_list, b_3_list = get_bends(1)

nld = nld_class.nonlin_dyn_class(
    gtpsa_prop, lat_prop, A_max, beta_inj, delta_max, b_3_list)

nld.zero_mult(lat_prop, 3)
nld.zero_mult(lat_prop, 4)

weight = get_weights()

prm_list, dprm_list = get_prms(1, bend_list, 1e-4)

opt_sp = opt_sp_class(
    lat_prop, prm_list, dprm_list, uc_list, sp_list, weight, bend_list,
    rbend_list, eps_x_des, phi_1_des, phi_rb_1_des, b_2_des, nu_uc_des,
    nu_sp_des, beta_des, dnu_des, nld)

opt_sp.opt_sp()
