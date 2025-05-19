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

eps_x_des   = 49e-12
alpha_c_des = -1e-4
nu_uc_des   = [0.4, 0.1]
nu_sp_des   = [2.89, 0.91]
beta_des    = [5.0, 3.0]
dnu_des     = [0.75, 0.25] # Phase advance across the straight.
dx_max      = 10e-3        # Max orbit.


@dataclass
class gtpsa_prop:
    # GTPSA properties.
    # Number of variables - phase-space coordinates & 1 for parameter
    #dependence.
    nv: ClassVar[int] = 6 + 1
    # Max order.
    no: ClassVar[int] = 1
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]


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


def compute_layout(lat_prop):
    s_buf = []
    X_buf = []
    Y_buf = []
    p_x_buf = []
    s = 0e0
    s_ref = s
    X = 0e0
    Y = 0e0
    p_x = 0e0
    s_buf.append(s)
    X_buf.append(X)
    Y_buf.append(Y)
    p_x_buf.append(p_x)
    for k, elem in enumerate(lat_prop._lattice):
        elem = lat_prop._lattice[k]
        L = elem.get_length()
        phi = math.radians(get_phi(elem))
        if phi == 0e0:
            X += L*np.cos(p_x)
            Y += L*np.sin(p_x)
        else:
            rho = L/phi
            X += rho*(np.sin(p_x+phi)-np.sin(p_x))
            Y += rho*(np.cos(p_x)-np.cos(p_x+phi))
        s += L
        p_x += phi
        # Cubic spline fit requires strictly monotonic series.
        if s > s_ref:
            s_buf.append(s)
            X_buf.append(X)
            Y_buf.append(Y)
            p_x_buf.append(p_x)
            s_ref = s
    s_buf = np.array(s_buf)
    X_buf = np.array(X_buf)
    Y_buf = np.array(Y_buf)
    p_x_buf = np.array(p_x_buf)
    X_cs = CubicSpline(s_buf, X_buf, bc_type="natural")
    Y_cs = CubicSpline(s_buf, Y_buf, bc_type="natural")
    return s_buf, X_cs, Y_cs, p_x_buf


def compute_alpha_c(map):
    C = lat_prop.compute_circ()
    index = np.zeros(gtpsa_prop.nv, int)
    alpha_c = np.zeros(gtpsa_prop.no+1)
    for k in range(1, gtpsa_prop.no+1):
        index[ind.delta] = k
        alpha_c[k] = map.ct.get(index)/C
    return alpha_c


class opt_sp_class:
    # Private

    def __init__(
            self, lat_ref, lat_prop, prm_list, uc_list, sp_list, weight,
            d1_bend, d2_bend, d2_2_bend, rb_list, eps_x_des, nu_uc_des,
            nu_sp_des, beta_des, dnu_des, nld):

        self._lat_prop       = lat_prop
        self._rf_cav_name    = "cav"
        self._cav_loc        = lat_prop._lattice.find(self._rf_cav_name, 0)
        self._nld            = nld
        self._prm_list       = prm_list
        self._uc_list        = uc_list
        self._sp_list        = sp_list
        self._weight         = weight
        self._d1_bend        = d1_bend
        self._d2_bend        = d2_bend
        self._d2_2_bend      = d2_2_bend
        self._rb_list        = rb_list
        self._eps_x_des      = eps_x_des
        self._nu_uc_des      = nu_uc_des
        self._nu_sp_des      = nu_sp_des
        self._beta_des       = beta_des
        self._dnu_des        = dnu_des

        self._Twiss_sp       = np.nan
        self._eta_list       = np.zeros((2, 4))
        self._alpha_list     = np.zeros((2, 2))
        self._nu_uc          = np.nan
        self._alpha_c        = np.nan
        self._nu_0           = np.nan
        self._dnu            = np.nan
        self._xi             = np.nan

        self._chi_2_min      = 1e30
        self._n_iter         = -1
        self._file_name      = "opt_sp.txt"

        self._h_im_scl_rms_0 = np.nan
        self._h_im_scl_rms_1 = np.nan
        self._K_re_scl_rms_0 = np.nan
        self._K_re_scl_rms_1 = np.nan

        self._s_ref          = np.nan
        self._X_ref          = np.nan
        self._Y_ref          = np.nan
        self._p_x_ref        = np.nan
        self._X              = np.nan
        self._Y              = np.nan
        self._DX             = np.nan
        self._DY             = np.nan
        self._dx             = np.nan

    # Public.

    def compute_phi_bend(self, bend_list):
        phi = 0e0
        for k in range(len(bend_list)):
            phi += self._lat_prop.get_phi_elem(bend_list[k], 0)
        return phi

    def compute_orbit(self, s_ref, X_ref, Y_ref, p_x_ref, X, Y, k):
        self._Dx = X(s_ref[k]) - X_ref(s_ref[k])
        self._Dy = Y(s_ref[k]) - Y_ref(s_ref[k])
        self._dx = np.sqrt(
            (self._Dx*np.sin(p_x_ref[k]))**2+(self._Dy*np.cos(p_x_ref[k]))**2)

    def prt_orbit(self):
        file_name = "orbit.txt"
        file = open(file_name, "w")

        print("# k      s            DX            DY            dx\n"
              "#                     [m]           [m]           [m]",
              file=file)
        for k in range(len(self._s_ref)):
            self.compute_orbit(
                self._s_ref, self._X_ref, self._Y_ref, self._p_x_ref, self._X,
                self._Y, k)
            print(f"{k:4d}  {self._s_ref[k]:9.5f}  {self._Dx:12.5e}"
                  f"  {self._Dy:12.5e}  {self._dx:12.5e}", file=file)

    def prt_phi_rb(self):
        for k, rb in enumerate(self._rb_list):
            print(f"    {rb:9s}      ="
                  f" {self._lat_prop.get_phi_elem(rb, 0):8.5f}")

    def prt_dip(self):
        phi = self._lat_prop.compute_phi_lat()
        phi_d1 = self.compute_phi_bend(self._d1_bend._bend_list)
        phi_d2 = self.compute_phi_bend(self._d2_bend._bend_list)
        phi_d2_2 = self.compute_phi_bend(self._d2_2_bend._bend_list)
        print(f"\n    phi_sp         = {phi:8.5f}")
        print(f"    phi_d1         = {phi_d1:8.5f}")
        print(f"    phi_d2         = {phi_d2:8.5f}")
        print(f"    phi_d2_2       = {phi_d2_2:8.5f}")
        self.prt_phi_rb()

    def prt_iter(self, prm, chi_2):
        eta, alpha, beta, nu_sp = self._Twiss_sp
        delta_hat = self._alpha_c[1]/(2e0*self._alpha_c[2])

        print("\n{:3d} chi_2 = {:9.3e} ({:9.3e})".format(
            self._n_iter, self._chi_2_min, chi_2-self._chi_2_min))
        print("    eps_x [pm.rad] = {:5.3f} [{:5.3f}]".
              format(1e12*self._lat_prop._eps[ind.X], 1e12*self._eps_x_des))

        print("\n    dx [mm]        = {:3.1f}".format(1e3*self._dx))

        print("\n    alpha_c        = [{:9.3e}, {:9.3e}] delta^ [%] = {:4.2f}".
              format(self._alpha_c[1], self._alpha_c[2], 1e2*delta_hat))
        print("    nu_uc          = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(self._nu_uc[ind.X], self._nu_uc[ind.Y],
                     self._nu_uc_des[ind.X], self._nu_uc_des[ind.Y]))
        print("    nu_sp          = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu_sp[ind.X], nu_sp[ind.Y], nu_sp_des[ind.X],
                     nu_sp_des[ind.Y]))
        print("    dnu            = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(self._dnu[ind.X], self._dnu[ind.Y],
                     self._dnu_des[ind.X], dnu_des[ind.Y]))
        print("    xi             = [{:5.3f}, {:5.3f}]".
              format(self._xi[ind.X], self._xi[ind.Y]))

        if len(uc_list) == 2:
            print("\n    eta'_uc        = {:9.3e}".
                  format(self._eta_list[0][ind.px]))
            print("    alpha_uc       = [{:9.3e}, {:9.3e}]".
                  format(self._alpha_list[0][ind.X],
                         self._alpha_list[0][ind.Y]))
        else:
            print("\n    eta'_uc        = {:9.3e} {:9.3e}".
                  format(self._eta_list[0][ind.px], self._eta_list[1][ind.px]))
            print("    alpha_uc       = [{:9.3e}, {:9.3e}] [{:9.3e}, {:9.3e}]".
                  format(self._alpha_list[0][ind.X], self._alpha_list[0][ind.Y],
                         self._alpha_list[1][ind.X], self._alpha_list[1][ind.Y]
                         ))
        print("    eta_x          = {:9.3e}".format(eta[ind.x]))
        print("    beta           = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(beta[ind.X], beta[ind.Y], self._beta_des[ind.X],
                    self._beta_des[ind.Y]))

        print("\n    C [m]          = {:8.5f}".
              format(self._lat_prop.compute_circ()))

        self.prt_dip()

        self._lat_prop.prt_rad()
        self._prm_list.prt_prm(prm)
        self.prt_orbit()

    def compute_chi_2(self):
        prt = not False

        eta, alpha, nu_sp, beta = self._Twiss_sp

        dchi_2 = weight[0]*(self._lat_prop._eps[ind.X]-self._eps_x_des)**2
        chi_2 = dchi_2
        if prt:
            print("\n  dchi2(eps_x)        = {:9.3e}".format(dchi_2))

        dchi_2 = weight[1]*(self._lat_prop._alpha_c-alpha_c_des)**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(alpha_c)      = {:9.3e}".format(dchi_2))

        dchi_2 = weight[2]*self._lat_prop._U_0**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(U_0)          = {:9.3e}".format(dchi_2))

        dchi_2 = weight[3]*(self._eta_list[0][ind.px]**2
                            +self._eta_list[1][ind.px]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(eta'_uc)      = {:9.3e}".format(dchi_2))

        dchi_2 = weight[4]*(
            self._alpha_list[0][ind.X]**2+self._alpha_list[0][ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(alpha_uc)     = {:9.3e}".format(dchi_2))

        dchi_2 = weight[4]*(
            self._alpha_list[1][ind.X]**2+self._alpha_list[1][ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(alpha_uc)     = {:9.3e}".format(dchi_2))

        dchi_2 = weight[5]*eta[ind.x]**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(eta_x)        = {:9.3e}".format(dchi_2))

        dchi_2 = weight[6]*(self._nu_uc[ind.X]-self._nu_uc_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_x)      = {:9.3e}".format(dchi_2))

        dchi_2 = weight[7]*(self._nu_uc[ind.Y]-self._nu_uc_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_y)      = {:9.3e}".format(dchi_2))

        dchi_2 = weight[8]*(self._dnu[ind.X]-self._dnu_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(dnu_x)        = {:9.3e}".format(dchi_2))

        dchi_2 = weight[8]*(self._dnu[ind.Y]-self._dnu_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(dnu_y)        = {:9.3e}".format(dchi_2))

        dchi_2 = weight[10]*(nu_sp[ind.X]-self._nu_sp_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_sp_x)      = {:9.3e}".format(dchi_2))

        dchi_2 = weight[11]*(nu_sp[ind.Y]-self._nu_sp_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_sp_y)      = {:9.3e}".format(dchi_2))

        dchi_2 = weight[12]*(beta[ind.X]-self._beta_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(beta_x)       = {:9.3e}".format(dchi_2))

        dchi_2 = weight[13]*(beta[ind.Y]-self._beta_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(beta_y)       = {:9.3e}".format(dchi_2))

        dchi_2 = weight[14]*(self._xi[ind.X]**2+self._xi[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(xi)           = {:9.3e}".format(dchi_2))

        dchi_2 = weight[15]*(self._dx-dx_max)**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(dx)           = {:9.3e}".format(dchi_2))

        if False:
            self._phi_1 = self.compute_phi_bend(self._d1_bend._bend_list)
            dchi_2 = 1e-2*(self._phi_d1-1.2)**2
            chi_2 += dchi_2
            if prt:
                print("  dchi2(phi_d1)       = {:9.3e}".format(dchi_2))

        return chi_2

    def opt_sp(self):
        """Use Case: optimise super period.
        """

        def f_sp(prm):
            self._n_iter += 1
            self._prm_list.set_prm(prm)

            try:
                # Compute Twiss parameters along the lattice.
                if not self._lat_prop.comp_per_sol():
                    print("\ncomp_per_sol: unstable")
                    raise ValueError

                # Adjust RF phase for sign of alpha_c.
                if self._lat_prop._alpha_c >= 0e0:
                    self._cav_loc.set_phase(0e0)
                else:
                    self._cav_loc.set_phase(180e0)
                    print("\nopt_sp:  alpha_c = {:10.3e} => phi_rf = 180 deg".
                          format(self._lat_prop._alpha_c))

                # Compute radiation properties.
                stable, stable_rad = self._lat_prop.compute_radiation()
                if not stable or not stable_rad:
                    print("\ncompute_radiation: unstable")
                    raise ValueError
 
                # Compute linear chromaticity.
                lat_prop._desc.truncate(2)
                self._nld.compute_map(lat_prop, 2)
                stable, _, self._xi = lo.compute_nu_xi(
                    lat_prop._desc, lat_prop._no, self._nld._M)
                if not stable:
                    print("\ncompute_nu_xi: unstable")
                    raise ValueError
            except ValueError:
                chi_2 = 1e30
                if not False:
                    print("\n{:3d} chi_2 = {:11.5e} ({:11.5e})".
                          format(self._n_iter, chi_2, self._chi_2_min))
                    if False:
                        self.prt_dip()
                        self._prm_list.prt_prm(prm)
            else:
                self._alpha_c = compute_alpha_c(self._nld._M)
                lat_prop._desc.truncate(1)
                _, _, _, self._nu_0 = self._lat_prop.get_Twiss(uc_list[0])
                self._eta_list[0], self._alpha_list[0], _, self._nu_1 = \
                    self._lat_prop.get_Twiss(uc_list[1])
                self._nu_uc = self._nu_1 - self._nu_0
                if len(uc_list) == 3:
                    self._eta_list[1], self._alpha_list[1], _, _ = \
                        self._lat_prop.get_Twiss(self._uc_list[1])
                self._Twiss_sp = self._lat_prop.get_Twiss(-1)
                self._dnu = 2e0*self._lat_prop.get_Twiss(sp_list[0])[3][:]

                _, self._X, self._Y, _ = compute_layout(lat_prop)
                self.compute_orbit(
                    self._s_ref, self._X_ref, self._Y_ref, self._p_x_ref,
                    self._X, self._Y, self._sp_list[2])

                chi_2 = self.compute_chi_2()

                if chi_2 < self._chi_2_min:
                    self.prt_iter(prm, chi_2)
                    if False:
                        self._lat_prop.prt_Twiss("twiss.txt")
                    pc.prt_lat(
                        self._lat_prop, self._file_name, self._prm_list,
                        d1_bend=self._d1_bend)
                    self._chi_2_min = min(chi_2, self._chi_2_min)
                else:
                    if not False:
                        print("\n{:3d} dchi_2 = {:9.3e} ({:9.3e})".format(
                            self._n_iter, chi_2-self._chi_2_min,
                            self._chi_2_min))
                        if not False:
                            self.prt_dip()
                        if False:
                            self._prm_list.prt_prm(prm)

            return chi_2

        max_iter = 1000

        self._s_ref, self._X_ref, self._Y_ref, self._p_x_ref = compute_layout(
            lat_ref)

        prm, bounds = self._prm_list.get_prm()
        f_sp(prm)

        def constraint_phi_lat(prm):
            phi = self._lat_prop.compute_phi_lat() - 18.0
            print(f"\nconstraint_phi_lat: phi = {phi:10.3e}")
            return phi

        def constraint_d1_min(prm):
            phi_min = 1.1
            phi = self.compute_phi_bend(self._d1_bend._bend_list) - phi_min
            print(f"\nconstraint_d1_min:  phi = {phi:10.3e}")
            return phi

        def constraint_d1_max(prm):         
            phi_max = 1.5
            phi = phi_max - self.compute_phi_bend(self._d1_bend._bend_list)
            print(f"\nconstraint_d1_max:  phi = {phi:10.3e}")
            return phi

        if True:
            constraints = ({"type": "eq",   "fun": constraint_phi_lat},
                           {"type": "ineq", "fun": constraint_d1_min},
                           {"type": "ineq", "fun": constraint_d1_max})
        else:
            constraints = ({"type": "eq",   "fun": constraint_phi_lat})

        f_tol = 1e-10
        g_tol = 1e-10
        x_tol = 1e-7
        eps   = 1e-4 # Step size for Jacobian.

        print("\nbounds =")
        for b in bounds:
            print(f"  ({b[0]:6.2f}, {b[1]:6.2f})")

        # Methods:
        #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
        #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.
        # Methods with boundaries:
        #   L-BFGS-B, TNC, and SLSQP.

        optimiser = "SLSQP"

        if optimiser == "SLSQP":
            minimum = opt.minimize(
                f_sp,
                prm,
                method="SLSQP",
                # callback=prt_iter,
                bounds = bounds,
                options={"ftol": f_tol, "maxiter": max_iter, "eps": 1e-3},
                constraints=constraints
            )
        elif optimiser == "CG":
            minimum = opt.minimize(
                f_sp,
                prm,
                method="CG",
                # callback=prt_iter,
                jac=None,
                options={"gtol": g_tol, "maxiter": max_iter,
                         # "eps": dprm_list
                         }
            )

        print("\n".join(minimum))


# TPSA max order.
gtpsa_prop.no = 2

cod_eps = 1e-10
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
delta_max = 3e-2
beta_inj  = np.array([3.0, 3.0])

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")

file_name_ref = os.path.join(home_dir, "max_iv", "max_iv_baseline_2.lat")
file_name = os.path.join(home_dir, sys.argv[1]+".lat")

lat_ref = lp.lattice_properties_class(gtpsa_prop, file_name_ref, E_0, cod_eps)

lat_prop = lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)
lat_prop.prt_lat("lat_prop_lat.txt")

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol: unstable")
        raise ValueError

    # Adjust RF phase for sign of alpha_c.
    cav_loc = lat_prop._lattice.find("cav", 0)
    if lat_prop._alpha_c >= 0e0:
        cav_loc.set_phase(0.0)
    else:
        print("  alpha_c = {:10.3e} phi_rf = 180 deg".
              format(lat_prop._alpha_c))
        cav_loc.set_phase(180.0)

   # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation: unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop.prt_lat_param()
    lat_prop.prt_rad()
    lat_prop.prt_M()
    lat_prop.prt_M_rad()

uc_list = []
uc_list.append(lat_prop._lattice.find("ucborder", 0).index)
uc_list.append(lat_prop._lattice.find("ucborder", 1).index)
uc_list.append(lat_prop._lattice.find("ucborder", 3).index)
uc_list = np.array(uc_list)

sp_list = np.zeros(3, dtype=int)
sp_list[0] = lat_prop._lattice.find("lsborder", 0).index
sp_list[1] = lat_prop._lattice.find("lsborder", 1).index
# Achromat centre - from cubic spline list.
sp_list[2] = lat_prop._lattice.find("d2_h2_sl_d0a", 5).index

print("\nunit cell entrance           {:15s} loc = {:d}".
      format(lat_prop._lattice[uc_list[0]].name, uc_list[0]))
print("unit cell exit               {:15s} loc = {:d}".
      format(lat_prop._lattice[uc_list[1]].name, uc_list[1]))
if len(uc_list) == 3:
    print("next unit cell exit          {:15s} loc = {:d}".
          format(lat_prop._lattice[uc_list[2]].name, uc_list[2]))
print("super period first sextupole {:15s} loc = {:d}".
      format(lat_prop._lattice[sp_list[0]].name, sp_list[0]))
print("super period last sextupole  {:15s} loc = {:d}".
      format(lat_prop._lattice[sp_list[1]].name, sp_list[1]))
print("max orbit location                           loc = {:d}".
      format(sp_list[2]))

d1_list = [
    "d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3", "d1_h2_sl_dm2",
    "d1_h2_sl_dm1", "d1_h2_sl_d0", "d1_h2_sl_ds1", "d1_h2_sl_ds2",
    "d1_h2_sl_ds3", "d1_h2_sl_ds4", "d1_h2_sl_ds5"
]
d2_list = ["d2_h2_sl_df0", "d2_h2_sl_df1", "d2_h2_sl_df2", "d2_h2_sl_df3",
           "d2_h2_sl_df4", "d2_h2_sl_df5", "d2_h2_sl_df6"]
d2_2_list = ["d2_2_h2_sl_df0", "d2_2_h2_sl_df1", "d2_2_h2_sl_df2",
             "d2_2_h2_sl_df3", "d2_2_h2_sl_df4", "d2_2_h2_sl_df5",
             "d2_2_h2_sl_df6"]

d1_bend   = pc.bend_class(lat_prop, d1_list)
d2_bend   = pc.bend_class(lat_prop, d2_list)
d2_2_bend = pc.bend_class(lat_prop, d2_2_list)

b_3_list = ["s3_h2", "s4_h2"]

nld = nld_class.nonlin_dyn_class(lat_prop, A_max, beta_inj, delta_max, b_3_list)
nld.zero_mult(lat_prop, 3)
nld.zero_mult(lat_prop, 4)

step = 1;

if step == 1:
    weight = np.array([
        1e15,  # eps_x.
        1e2,   # alpha_c.
        0e-15, # U_0.
        1e2,   # etap_x_uc.
        1e-2,  # alpha_uc.
        1e2,   # eta_x.
        1e1,   # nu_uc_x.
        1e1,   # nu_uc_y.
        1e-3,  # dnu_x.
        1e-3,  # dnu_y.
        0e-7,  # nu_sp_x.
        0e-7,  # nu_sp_y.
        0e-6,  # beta_x.
        0e-6,  # beta_y.
        1e-9,  # xi.
        0e-3   # orbit.
    ])

    prms = [
        ("q1_h2",   "b_2",      -10.0, 10.0),
        ("q2_h2",   "b_2",      -10.0, 10.0),
        ("r1_h2",   "b_2",      -10.0, 10.0),
        ("r2_h2",   "b_2",      -10.0, 10.0),
        ("r3_h2",   "b_2",      -10.0, 10.0),
        (d1_bend,   "b_2_bend", -1.5,   1.5),
        (d2_bend,   "b_2_bend", -1.5,   1.5),
        (d2_2_bend, "b_2_bend", -1.5,   1.5),

        # ("r1_h2",   "phi",      -0.5,   0.5),
        # ("r2_h2",   "phi",      -0.3,   0.0),
        # ("r3_h2",   "phi",      -0.3,   0.0),
        (d1_bend,   "phi_bend",  1.1,   1.5),
        (d2_bend,   "phi_bend",  1.5,   1.75),
        (d2_2_bend, "phi_bend",  1.5,   2.5)
   ]

    rb_list = ["r1_h2", "r2_h2", "r3_h2"]
if step == 2:
    weight = np.array([
        1e15,  # eps_x.
        1e3,   # alpha_c.
        0e-15, # U_0.
        1e2,   # etap_x_uc.
        1e-2,  # alpha_uc.
        1e2,   # eta_x.
        1e1,   # nu_uc_x.
        1e1,   # nu_uc_y.
        1e-3,  # dnu_x.
        1e-3,  # dnu_y.
        0e-7,  # nu_sp_x.
        0e-7,  # nu_sp_y.
        0e-6,  # beta_x.
        0e-6,  # beta_y.
        1e-7,  # xi.
        0e-3   # orbit.
    ])

    prms = [
        ("q1_h2",          "b_2", -10.0, 10.0),
        ("q2_h2",          "b_2", -10.0, 10.0),
        ("r1_h2",          "b_2", -10.0, 10.0),
        ("r2_h2",          "b_2", -10.0, 10.0),
        ("r3_h2",          "b_2", -10.0, 10.0),

        ("d1_h2_sl_dm5",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_dm4",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_dm3",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_dm2",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_dm1",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_d0",    "b_2",  -1.5,  1.5),
        ("d1_h2_sl_ds1",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_ds2",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_ds3",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_ds4",   "b_2",  -1.5,  1.5),
        ("d1_h2_sl_ds5",   "b_2",  -1.5,  1.5),

        ("d2_h2_sl_df0",   "b_2",  -1.5,  1.5),
        ("d2_h2_sl_df1",   "b_2",  -1.5,  1.5),
        ("d2_h2_sl_df2",   "b_2",  -1.5,  1.5),
        ("d2_h2_sl_df3",   "b_2",  -1.5,  1.5),
        ("d2_h2_sl_df4",   "b_2",  -1.5,  1.5),
        ("d2_h2_sl_df5",   "b_2",  -1.5,  1.5),
        ("d2_h2_sl_df6",   "b_2",  -1.5,  1.5),

        ("d2_2_h2_sl_df0", "b_2",  -1.5,  1.5),
        ("d2_2_h2_sl_df1", "b_2",  -1.5,  1.5),
        ("d2_2_h2_sl_df2", "b_2",  -1.5,  1.5),
        ("d2_2_h2_sl_df3", "b_2",  -1.5,  1.5),
        ("d2_2_h2_sl_df4", "b_2",  -1.5,  1.5),
        ("d2_2_h2_sl_df5", "b_2",  -1.5,  1.5),
        ("d2_2_h2_sl_df6", "b_2",  -1.5,  1.5),

        ("r1_h2",          "phi",  -0.5,   0.5),
        ("r2_h2",          "phi",  -0.3,   0.0),
        ("r3_h2",          "phi",  -0.3,   0.0),

        ("d1_h2_sl_dm5",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_dm4",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_dm3",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_dm2",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_dm1",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_d0",    "phi",  -1.5,  1.5),
        ("d1_h2_sl_ds1",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_ds2",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_ds3",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_ds4",   "phi",  -1.5,  1.5),
        ("d1_h2_sl_ds5",   "phi",  -1.5,  1.5),

        ("d2_h2_sl_df0",   "phi",  -1.5,  1.5),
        ("d2_h2_sl_df1",   "phi",  -1.5,  1.5),
        ("d2_h2_sl_df2",   "phi",  -1.5,  1.5),
        ("d2_h2_sl_df3",   "phi",  -1.5,  1.5),
        ("d2_h2_sl_df4",   "phi",  -1.5,  1.5),
        ("d2_h2_sl_df5",   "phi",  -1.5,  1.5),
        ("d2_h2_sl_df6",   "phi",  -1.5,  1.5),

        ("d2_2_h2_sl_df0", "phi",  -1.5,  1.5),
        ("d2_2_h2_sl_df1", "phi",  -1.5,  1.5),
        ("d2_2_h2_sl_df2", "phi",  -1.5,  1.5),
        ("d2_2_h2_sl_df3", "phi",  -1.5,  1.5),
        ("d2_2_h2_sl_df4", "phi",  -1.5,  1.5),
        ("d2_2_h2_sl_df5", "phi",  -1.5,  1.5),
        ("d2_2_h2_sl_df6", "phi",  -1.5,  1.5)
   ]

    rb_list = ["r1_h2", "r2_h2", "r3_h2"]

prm_list = pc.prm_class(lat_prop, prms)

opt_sp = opt_sp_class(
    lat_ref, lat_prop, prm_list, uc_list, sp_list, weight, d1_bend,
    d2_bend, d2_2_bend, rb_list, eps_x_des, nu_uc_des, nu_sp_des, beta_des,
    dnu_des, nld)

opt_sp.opt_sp()
