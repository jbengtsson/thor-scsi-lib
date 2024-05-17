"""Use Case:
     Implementation and general optimisation of a higher-order-achromat.
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


def compute_periodic_cell(lat_prop, uc_0, uc_1):
    M = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    M.set_identity()
    # The 3rd argument is the 1st element index & the 4th the number of elements
    # to propagate across.
    lat_prop._lattice.propagate(
        lat_prop._model_state, M, uc_0, uc_1-uc_0+1)
    # Problem with nu.
    stable, nu, A, A_inv, _ = lo.compute_M_diag(2, M.jacobian())
    eta, Twiss = lo.transform_matrix_extract_twiss(A)
    alpha = np.array((Twiss[ind.X][0], Twiss[ind.Y][0]))
    beta = np.array((Twiss[ind.X][1], Twiss[ind.Y][1]))
    # Problem with nu/1-nu for compute_nus.
    # nu = lo.compute_nus(2, M.jacobian())
    nu = lo.compute_nu_symp(2, M.jacobian())
    Twiss = eta, alpha, beta, nu
    return stable, Twiss


def opt_sp(lat_prop, prm_list, uc_0, uc_1, weight, b1_list, b2_list):
    """Use Case: optimise super period.
    """

    nu_uc_des    = [1.0/3.0, 1.0/12.0]

    chi_2_min    = 1e30
    n_iter       = 0
    file_name    = "opt_sp.txt"

    Twiss_uc     = np.nan
    Twiss_uc_ref = np.nan
    nu_sp        = np.nan

    def prt_iter(prm, chi_2, Twiss_uc, nu_sp, xi):
        nonlocal Twiss_uc_ref, nu_uc_des

        eta_uc, alpha_uc, beta_uc, nu_uc = Twiss_uc
        eta_uc_ref, alpha_uc_ref, beta_uc_ref, nu_uc_ref = Twiss_uc_ref

        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        rb_1 = "qf1"
        rb_2 = "qf1_e"

        phi = lat_prop.compute_phi_lat()
        phi_b1 = compute_phi_bend(lat_prop, b1_list)
        phi_b2 = compute_phi_bend(lat_prop, b2_list)
        phi_rb_1 = lat_prop.get_phi_elem(rb_1, 0)
        phi_rb_2 = lat_prop.get_phi_elem(rb_2, 0)

        print("\n{:3d} chi_2 = {:16.10e}".format(n_iter, chi_2))
        print("    eps_x [pm.rad] = {:5.3f}".format(1e12*lat_prop._eps[ind.X]))
        print("    U_0 [keV]      = {:5.3f}".format(1e-3*lat_prop._U_0))
        print("\n    eta_x_uc       = [{:9.3e}, {:9.3e}] ([{:9.3e}, {:9.3e}])".
              format(eta_uc[ind.x], eta_uc[ind.px], eta_uc_ref[ind.x],
                     eta_uc_ref[ind.px]))
        print("    alpha_uc       = [{:9.3e}, {:9.3e}] ([{:9.3e}, {:9.3e}]))".
              format(alpha_uc[ind.X], alpha_uc[ind.Y], alpha_uc_ref[ind.X],
                     alpha_uc_ref[ind.Y]))
        print("    beta_uc        = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(beta_uc[ind.X], beta_uc[ind.Y], beta_uc_ref[ind.X],
                     beta_uc_ref[ind.Y]))
        print("    nu_uc_ref      = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu_uc_ref[ind.X], nu_uc_ref[ind.Y], nu_uc_des[ind.X],
                     nu_uc_des[ind.Y]))
        print("\n    nu_sp          = [{:7.5f}, {:7.5f}]".
              format(nu_sp[ind.X], nu_sp[ind.Y]))
        print("    xi             = [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        print("\n    phi_sp         = {:8.5f}".format(phi))
        print("    C [m]          = {:8.5f}".format(lat_prop.compute_circ()))
        print("\n    phi_b1         = {:8.5f}".format(phi_b1))
        print("    phi_b2         = {:8.5f}".format(phi_b2))
        print("    phi_rb_1       = {:8.5f}".format(phi_rb_1))
        print("    phi_rb_2       = {:8.5f}".format(phi_rb_2))
        prm_list.prt_prm(prm)

    def compute_chi_2(Twiss_uc, xi):
        nonlocal Twiss_uc_ref, nu_uc_des

        prt = not False

        eta_uc, alpha_uc, beta_uc, nu_uc = Twiss_uc
        eta_uc_ref, alpha_uc_ref, beta_uc_ref, nu_uc_ref = Twiss_uc_ref

        dchi_2 = weight[0]*lat_prop._eps[ind.X]**2
        chi_2 = dchi_2
        if prt:
            print("\n  dchi2(eps_x)     = {:10.3e}".format(dchi_2))

        dchi_2 = weight[1]*lat_prop._U_0**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(U_0)       = {:10.3e}".format(dchi_2))

        dchi_2 = \
            weight[2]*((eta_uc[ind.x]-eta_uc_ref[ind.x])**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(eta_x_uc)  = {:10.3e}".format(dchi_2))

        dchi_2 = \
            weight[3]*((eta_uc[ind.px]-eta_uc_ref[ind.px])**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(eta'_x_uc) = {:10.3e}".format(dchi_2))

        dchi_2 = \
            weight[4]*((alpha_uc[ind.X]-alpha_uc_ref[ind.X])**2
                       +(alpha_uc[ind.Y]-alpha_uc_ref[ind.Y])**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(alpha_uc)  = {:10.3e}".format(dchi_2))

        dchi_2 = \
            weight[5]*((beta_uc[ind.X]-beta_uc_ref[ind.X])**2
                       +(beta_uc[ind.Y]-beta_uc_ref[ind.Y])**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(beta_uc)   = {:10.3e}".format(dchi_2))

        dchi_2 = \
            weight[6]*(
                (nu_uc_ref[ind.X]-nu_uc_des[ind.X])**2
                +(nu_uc_ref[ind.Y]-nu_uc_des[ind.Y])**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc)     = {:10.3e}".format(dchi_2))

        dchi_2 = weight[7]*(xi[ind.X]**2+xi[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(xi)        = {:10.3e}".format(dchi_2))

        return chi_2

    def f_sp(prm):
        nonlocal chi_2_min, n_iter, Twiss_uc_ref, file_name

        n_iter += 1
        prm_list.set_prm(prm)

        try:
            # Compute Twiss parameters along lattice.
            if not lat_prop.comp_per_sol():
                print("\ncomp_per_sol : unstable")
                raise ValueError

            # Compute radiation properties.
            if not lat_prop.compute_radiation():
                print("\ncompute_radiation : unstable")
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
        else:
            eta_uc, alpha_uc, beta_uc, nu_0 = lat_prop.get_Twiss(uc_0-1)
            _, _, _, nu_1 = lat_prop.get_Twiss(uc_1)
            nu_uc = nu_1 - nu_0
            Twiss_uc = eta_uc, alpha_uc, beta_uc, nu_uc

            _, _, _, nu_sp = lat_prop.get_Twiss(-1)

            chi_2 = compute_chi_2(Twiss_uc, xi)

        if chi_2 < chi_2_min:
            prt_iter(prm, chi_2, Twiss_uc, nu_sp, xi)
            pc.prt_lat(lat_prop, file_name, prm_list)
            chi_2_min = min(chi_2, chi_2_min)

            # stable, Twiss_uc_ref = \
            #     compute_periodic_cell(lat_prop, uc_0, uc_1)
            # Twiss_uc_ref[0][ind.px] = 0e0
            # Twiss_uc_ref[1][ind.X], Twiss_uc_ref[1][ind.Y] = \
            #     [0e0, 0e0]

        return chi_2

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    stable, Twiss_uc_ref = compute_periodic_cell(lat_prop, uc_0, uc_1)
    Twiss_uc_ref[0][ind.px] = 0e0
    Twiss_uc_ref[1][ind.X], Twiss_uc_ref[1][ind.Y] = [0e0, 0e0]

    prm, bounds = prm_list.get_prm()

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
lat_name = "max_4u_sp_jb_2"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol - unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop.prt_lat_param()
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

uc_0 = lat_prop._lattice.find("sf_h", 1).index
uc_1 = lat_prop._lattice.find("sf_h", 2).index
print("\nunit cell entrance {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_0].name, uc_0))
print("unit cell exit     {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_1].name, uc_1))

# Weights.
weight = np.array([
    1e15,  # eps_x.
    1e-16, # U_0.
    1e-1,  # eta_x_uc.
    1e-1,  # eta'_x_uc.
    1e-6,  # alpha_uc.
    1e-5,  # beta_uc.
    1e1,   # nu_uc.
    1e-6   # xi.
])

b1_list = ["b1_0", "b1_1", "b1_2", "b1_3", "b1_4", "b1_5"]

b2_list = [
    "b2u_6", "b2u_5", "b2u_4", "b2u_3", "b2u_2", "b2u_1", "b2_0",
    "b2d_1", "b2d_2", "b2d_3", "b2d_4", "b2d_5"
]

b1_bend = pc.bend_class(lat_prop, b1_list, phi_max, b_2_max)
b2_bend = pc.bend_class(lat_prop, b2_list, phi_max, b_2_max)

# The end dipole bend angle is the parameter whereas the unit cell dipole bend
# angle is used for maintaining the total bend angle.
phi_n_list = [2, 2]
phi_list   = [b2_bend, b1_bend]

opt_phi = pc.phi_tot_class(lat_prop, phi_list, phi_n_list, phi_max)

prm_list = [
    ("qf1",      "b_2"),
    ("qf1_e",    "b_2"),
    ("qd",       "b_2"),
    ("qf2",      "b_2"),
    # ("b_2_bend", b1_bend),
    # ("b_2_bend", b2_bend),

    ("b1_0",     "b_2"),
    ("b1_1",     "b_2"),
    ("b1_2",     "b_2"),
    ("b1_3",     "b_2"),
    ("b1_4",     "b_2"),
    ("b1_5",     "b_2"),

    # ("b2u_6",    "b_2"),
    # ("b2u_5",    "b_2"),
    # ("b2u_4",    "b_2"),
    # ("b2u_3",    "b_2"),
    # ("b2u_2",    "b_2"),
    # ("b2u_1",    "b_2"),
    # ("b2_0",     "b_2"),
    # ("b2d_1",    "b_2"),
    # ("b2d_2",    "b_2"),
    # ("b2d_3",    "b_2"),
    # ("b2d_4",    "b_2"),
    # ("b2d_5",    "b_2")

    # ("qf1",      "phi"),
    # ("qf1_e",    "phi"),
    # ("phi_tot",  opt_phi)
]

prm_list = pc.prm_class(lat_prop, prm_list, b_2_max)

opt_sp(lat_prop, prm_list, uc_0, uc_1, weight, b1_list, b2_list)
