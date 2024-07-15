"""Use Case:
     Implementation and general optimisation of a higher-order-achromat.
     Use one unit cell.
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
    linear_optics as lo, prm_class as pc
from thor_scsi.utils.output import mat2txt, vec2txt

import gtpsa


ind = ind.index_class()

L_min        = 0.10
L_max        = 0.50
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0

eps_x_des    = 139e-12
nu_uc_des    = [2.0/6.0, 0.12]
nu_sp_des    = [2.25, 0.8]
beta_des     = [5.7, 2.0]
dnu_des      = [0.5, 0.25]     # Phase advance across the straight.


def opt_sp(
        lat_prop, prm_list, uc_0, uc_1, uc_2, sp_1, sp_2, weight, b1_list,
        b2_list, phi_lat, eps_x_des, nu_uc_des, nu_sp_des, beta_des, dnu_des):
    """Use Case: optimise super period.
    """

    chi_2_min    = 1e30
    n_iter       = -1
    file_name    = "opt_sp.txt"

    Twiss_uc     = np.nan
    Twiss_uc_ref = np.nan
    nu_sp        = np.nan

    def prt_iter(
            prm, chi_2, Twiss_sp, eta_uc_1, alpha_uc_1, eta_uc_2, alpha_uc_2,
            nu_uc, dnu, xi):
        nonlocal nu_uc_des, nu_sp_des, beta_des

        eta, alpha, beta, nu_sp = Twiss_sp

        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        rb = "qf1"

        phi = lat_prop.compute_phi_lat()
        phi_b1 = compute_phi_bend(lat_prop, b1_list)
        phi_b2 = compute_phi_bend(lat_prop, b2_list)
        phi_rb = lat_prop.get_phi_elem(rb, 0)

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("    eps_x [pm.rad] = {:5.3f} [{:5.3f}]".
              format(1e12*lat_prop._eps[ind.X], 1e12*eps_x_des))
        print("    U_0 [keV]      = {:5.3f}".format(1e-3*lat_prop._U_0))
        print("\n    eta'_uc        = {:9.3e} {:9.3e}".
              format(eta_uc_1[ind.px], eta_uc_2[ind.px]))
        print("    alpha_uc       = [{:9.3e}, {:9.3e}] [{:9.3e}, {:9.3e}]".
              format(alpha_uc_1[ind.X], alpha_uc_1[ind.Y], alpha_uc_2[ind.X],
                     alpha_uc_2[ind.Y]))
        print("    nu_uc          = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu_uc[ind.X], nu_uc[ind.Y], nu_uc_des[ind.X],
                     nu_uc_des[ind.Y]))
        print("\n    eta_x          = {:9.3e}".format(eta[ind.x]))
        print("\n    nu_sp          = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu_sp[ind.X], nu_sp[ind.Y], nu_sp_des[ind.X],
                     nu_sp_des[ind.Y]))
        print("    beta           = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(beta[ind.X], beta[ind.Y], beta_des[ind.X],
                     beta_des[ind.Y]))
        print("    dnu            = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(dnu[ind.X], dnu[ind.Y], dnu_des[ind.X], dnu_des[ind.Y]))
        print("    xi             = [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        print("\n    phi_sp         = {:8.5f}".format(phi))
        print("    C [m]          = {:8.5f}".format(lat_prop.compute_circ()))
        print("\n    phi_b1         = {:8.5f}".format(phi_b1))
        print("    phi_b2         = {:8.5f}".format(phi_b2))
        print("    phi_rb         = {:8.5f}".format(phi_rb))
        prm_list.prt_prm(prm)

    def compute_chi_2(
            Twiss_sp, eta_uc_1, alpha_uc_1, eta_uc_2, alpha_uc_2, nu_uc, dnu,
            xi):
        nonlocal nu_uc_des, nu_sp_des, beta_des

        prt = not False

        eta, alpha, nu_sp, beta = Twiss_sp

        dchi_2 = weight[0]*(lat_prop._eps[ind.X]-eps_x_des)**2
        chi_2 = dchi_2
        if prt:
            print("\n  dchi2(eps_x)    = {:9.3e}".format(dchi_2))

        dchi_2 = weight[1]*lat_prop._U_0**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(U_0)      = {:9.3e}".format(dchi_2))

        dchi_2 = weight[2]*(eta_uc_1[ind.px]**2+eta_uc_2[ind.px]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(eta'_uc)  = {:9.3e}".format(dchi_2))

        dchi_2 = weight[3]*(
            alpha_uc_1[ind.X]**2+alpha_uc_2[ind.X]**2+alpha_uc_1[ind.Y]**2
            +alpha_uc_2[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(alpha_uc) = {:9.3e}".format(dchi_2))

        dchi_2 = weight[4]*(nu_uc[ind.X]-nu_uc_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_x)  = {:9.3e}".format(dchi_2))

        dchi_2 = weight[5]*(nu_uc[ind.Y]-nu_uc_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_y)  = {:9.3e}".format(dchi_2))

        dchi_2 = weight[6]*eta[ind.x]**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(eta_x)    = {:9.3e}".format(dchi_2))

        dchi_2 = weight[7]*(nu_sp[ind.X]-nu_sp_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_sp_x)  = {:9.3e}".format(dchi_2))

        dchi_2 = weight[8]*(nu_sp[ind.Y]-nu_sp_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_sp_y)  = {:9.3e}".format(dchi_2))

        dchi_2 = weight[9]*(beta[ind.X]-beta_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(beta_x)   = {:9.3e}".format(dchi_2))

        dchi_2 = weight[10]*(beta[ind.Y]-beta_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(beta_y)   = {:9.3e}".format(dchi_2))

        dchi_2 = weight[11]*(dnu[ind.X]-dnu_des[ind.X])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(dnu_x)    = {:9.3e}".format(dchi_2))

        dchi_2 = weight[12]*(dnu[ind.Y]-dnu_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(dnu_y)    = {:9.3e}".format(dchi_2))

        dchi_2 = weight[13]*xi[ind.X]**2
        if prt:
            print("  dchi2(xi_x)     = {:9.3e}".format(dchi_2))

        dchi_2 = weight[14]*xi[ind.Y]**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(xi_y)     = {:9.3e}".format(dchi_2))

        return chi_2

    def f_sp(prm):
        nonlocal chi_2_min, n_iter, file_name

        n_iter += 1
        prm_list.set_prm(prm)
        phi_lat.set_phi_lat()

        try:
            # Compute Twiss parameters along the lattice.
            if not lat_prop.comp_per_sol():
                print("\ncomp_per_sol: unstable")
                raise ValueError

            # Compute radiation properties.
            if not lat_prop.compute_radiation():
                print("\ncompute_radiation: unstable")
                raise ValueError

            # Compute linear chromaticity.
            M = lo.compute_map(
                lat_prop._lattice, lat_prop._model_state,
                desc=lat_prop._desc, tpsa_order=2)
            stable, _, xi = \
                lo.compute_nu_xi(lat_prop._desc, lat_prop._no, M)
            if not stable:
                print("\ncompute_nu_xi: unstable")
                raise ValueError
        except ValueError:
            chi_2 = 1e30
            if False:
                print("\n{:3d} chi_2 = {:11.5e} ({:11.5e})".
                      format(n_iter, chi_2, chi_2_min))
                prm_list.prt_prm(prm)
        else:
            _, _, _, nu_0 = lat_prop.get_Twiss(uc_0-1)
            eta_uc_1, alpha_uc_1, _, nu_1 = lat_prop.get_Twiss(uc_1)
            nu_uc = nu_1 - nu_0
            eta_uc_2, alpha_uc_2, _, _ = lat_prop.get_Twiss(uc_2)
            Twiss_sp = lat_prop.get_Twiss(-1)

            dnu = \
                lat_prop.get_Twiss(-1)[3][:] - lat_prop.get_Twiss(sp_2)[3][:] \
                + lat_prop.get_Twiss(sp_1)[3][:]

            chi_2 = compute_chi_2(
                Twiss_sp, eta_uc_1, alpha_uc_1, eta_uc_2, alpha_uc_2, nu_uc,
                dnu, xi)

            if chi_2 < chi_2_min:
                prt_iter(
                    prm, chi_2, Twiss_sp, eta_uc_1, alpha_uc_1, eta_uc_2,
                    alpha_uc_2, nu_uc, dnu, xi)
                lat_prop.prt_Twiss("twiss.txt")
                pc.prt_lat(lat_prop, file_name, prm_list, phi_lat=phi_lat)
                chi_2_min = min(chi_2, chi_2_min)

        return chi_2

    max_iter = 1000
    f_tol    = 1e-10
    x_tol    = 1e-4
    g_tol    = 1e-5

    prm, bounds = prm_list.get_prm()
    f_sp(prm)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol, xtol
    # CG     gtol
    minimum = opt.minimize(
        f_sp,
        prm,
        method="CG",
        # callback=prt_iter,
        # bounds = bounds,
        # options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
        options={"gtol": g_tol, "maxiter": max_iter}
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
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_4u")
# lat_name = "max_iv_sp_jb"
lat_name = "max_4u_f_2"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

lat_prop.prt_lat(lat_name+"_lat.txt")

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol: unstable")
        raise ValueError

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

uc_0 = lat_prop._lattice.find("sf_h", 1).index
uc_1 = lat_prop._lattice.find("sf_h", 2).index
uc_2 = lat_prop._lattice.find("sf_h", 4).index
sp_1 = lat_prop._lattice.find("sd2", 0).index
sp_2 = lat_prop._lattice.find("sd2", 1).index
# For 2 sf [sf1, sf2] families.
# uc_0 = lat_prop._lattice.find("sf_h", 0).index
# uc_1 = lat_prop._lattice.find("sf_h", 1).index
print("\nunit cell entrance           {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_0].name, uc_0))
print("unit cell exit               {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_1].name, uc_1))
print("next unit cell exit          {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_2].name, uc_2))
print("super period first sextupole {:5s} loc = {:d}".
      format(lat_prop._lattice[sp_1].name, sp_1))
print("super period last sextupole  {:5s} loc = {:d}".
      format(lat_prop._lattice[sp_1].name, sp_2))

# Weights.
weight = np.array([
    1e17,  # eps_x.
    1e-17, # U_0.
    1e2,   # etap_x_uc.
    1e-1,  # alpha_uc.
    0e0,   # nu_uc_x.
    0e0,   # nu_uc_y.
    1e1,   # eta_x.
    0e-7,  # nu_sp_x.
    0e-3,  # nu_sp_y.
    0e-6,  # beta_x.
    0e-6,  # beta_y.
    1e-2,  # dnu_x.
    1e-2,  # dnu_y.
    1e-7,  # xi_x.
    1e-7   # xi_y.
])

b1_list = ["b1_0", "b1_1", "b1_2", "b1_3", "b1_4", "b1_5"]

b2_list = [
    "b2u_6", "b2u_5", "b2u_4", "b2u_3", "b2u_2", "b2u_1", "b2_0",
    "b2d_1", "b2d_2", "b2d_3", "b2d_4", "b2d_5"
]

b1_bend = pc.bend_class(lat_prop, b1_list, phi_max, b_2_max)
b2_bend = pc.bend_class(lat_prop, b2_list, phi_max, b_2_max)

if True:
    prms = [
        ("qf1",      "b_2"),
        ("qf1_e",    "b_2"),
        ("qd",       "b_2"),
        ("qf2",      "b_2"),

        ("b_2_bend", b1_bend),
        ("b_2_bend", b2_bend),

        ("phi_bend", b1_bend),
        ("qf1",      "phi")
    ]

    # To maintain the total bend angle.
    phi_lat = pc.phi_lat_class(lat_prop, 2, b2_bend)
else:
    prms = [
        ("qf1",      "b_2"),
        ("qf1_e",    "b_2"),
        ("qd",       "b_2"),
        ("qf2",      "b_2"),

        ("b1_0",     "b_2"),
        ("b1_1",     "b_2"),
        ("b1_2",     "b_2"),
        ("b1_3",     "b_2"),
        ("b1_4",     "b_2"),
        ("b1_5",     "b_2"),

        ("b2u_6",    "b_2"),
        ("b2u_5",    "b_2"),
        ("b2u_4",    "b_2"),
        ("b2u_3",    "b_2"),
        ("b2u_2",    "b_2"),
        ("b2u_1",    "b_2"),
        ("b2_0",     "b_2"),
        ("b2d_1",    "b_2"),
        ("b2d_2",    "b_2"),
        ("b2d_3",    "b_2"),
        ("b2d_4",    "b_2"),
        ("b2d_5",    "b_2"),

        ("b1_0",     "phi"),
        ("b1_1",     "phi"),
        ("b1_2",     "phi"),
        ("b1_3",     "phi"),
        ("b1_4",     "phi"),
        ("b1_5",     "phi"),

        ("qf1",      "phi")
    ]

    # To maintain the total bend angle.
    phi_lat = pc.phi_lat_class(lat_prop, 2, b2_bend)

prm_list = pc.prm_class(lat_prop, prms, b_2_max)

opt_sp(
    lat_prop, prm_list, uc_0, uc_1, uc_2, sp_1, sp_2, weight, b1_list, b2_list,
    phi_lat, eps_x_des, nu_uc_des, nu_sp_des, beta_des, dnu_des)
