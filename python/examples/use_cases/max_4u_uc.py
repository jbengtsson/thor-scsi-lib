"""Use Case:
     Module for implementing & optimising a unit cell.
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

import gtpsa

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, index_class as ind, \
    linear_optics as lo, knobs as kb, prm_class as pc, nonlin_dyn as nld_class
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()

L_min        = 0.10
L_max        = 0.50
phi_max      = 0.85
b_2_bend_max = 1.0
b_2_max      = 10.0

eps_x_des    = 10e-12
nu_uc_des    = [0.49, 0.24]
eta_x_des    = 1e-4


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


def compute_linear_optics(lat_prop, nld):
    # Compute the beam dynamics properties.
    alpha_rad_z_min = -1e-8

    try:
        # Compute Twiss parameters along lattice.
        if not lat_prop.comp_per_sol():
            print("\ncomp_per_sol: unstable")
            raise ValueError
        else:
            _, _, beta, nu_uc = lat_prop.get_Twiss(len(lat_prop._lattice)-1)
            eta, _, _, _ = lat_prop.get_Twiss(uc_centre)

        # Adjust RF phase for sign of alpha_c.
        cav_loc = lat_prop._lattice.find("cav", 0)
        if lat_prop._alpha_c >= 0e0:
            cav_loc.set_phase(0e0)
        else:
            cav_loc.set_phase(180e0)
            print("\nopt_sp:  alpha_c = {:10.3e} => phi_rf = 180 deg".
                  format(lat_prop._alpha_c))

        # Compute radiation properties.
        stable, stable_rad = lat_prop.compute_radiation()
        if not stable:
            print("\ncompute_radiation: unstable")
            raise ValueError
        if lat_prop._alpha_rad[ind.Z] > alpha_rad_z_min:
            print("\ncompute_linear_optics:\n"+
                  " unstable in the longitudinal plane alpha_rad_z = {:9.3e}".
                  format(lat_prop._alpha_rad[ind.Z]))
            raise ValueError

        # Compute linear chromaticity.
        nld.compute_map(lat_prop, 2)
        stable, _, xi = lo.compute_nu_xi(lat_prop._desc, lat_prop._no, nld._M)
        if not stable:
            print("\nf_sp - compute_nu_xi: unstable")
            raise ValueError
    except ValueError:
        return False, np.nan, np.nan, np.nan, np.nan
    else:
        return True, eta, beta, nu_uc, xi


def compute_Twiss(M):
    alpha = np.zeros(2)
    beta = np.zeros(2)
    nu = np.zeros(2)
    for k in range(2):
        tr = np.trace(M[2*k:2*k+2, 2*k:2*k+2])
        nu[k] = np.arccos(tr/2e0)/(2e0*np.pi)
        if M[2*k, 2*k+1] < 0e0:
            nu[k] = 1e0 - nu[k]
        alpha[k] = (M[2*k, 2*k]-np.cos(2e0*np.pi*nu[k]))/np.sin(2e0*np.pi*nu[k])
        beta[k] = M[2*k, 2*k+1]/np.sin(2e0*np.pi*nu[k])
    return alpha, beta, nu


def opt_uc(lat_prop, prm_list, uc_centre, weight, sp_bend, rb_1, phi_b,
           eps_x_des, nu_uc_des, nld):
    """Use Case: optimise super period.
    """
    chi_2_min = 1e30
    n_iter    = 0
    file_name = "opt_uc.txt"

    def prt_iter(prm, chi_2, eta, beta, nu_uc, xi):
        nonlocal n_iter

        def compute_phi_bend(lat_prop, bend_list):
            phi = 0e0
            for k in range(len(bend_list)):
                phi += lat_prop.get_phi_elem(bend_list[k], 0)
            return phi

        phi = lat_prop.compute_phi_lat()
        phi_rb_1 = lat_prop.get_phi_elem(rb_1, 0)
        phi_bend = lat_prop.get_phi_elem(phi_b, 0)

        print("\n{:3d} chi_2 = {:11.5e}".format(n_iter, chi_2))
        print("    eps_x [pm.rad] = {:5.3f} [{:5.3f}]".
              format(1e12*lat_prop._eps[ind.X], 1e12*eps_x_des))
        print("    alpha_c        = {:9.3e}".format(lat_prop._alpha_c))
        print("    nu_uc          = [{:7.5f}, {:7.5f}] ([{:7.5f}, {:7.5f}])".
              format(nu_uc[ind.X], nu_uc[ind.Y], nu_uc_des[ind.X],
                     nu_uc_des[ind.Y]))
        print("    beta           = [{:5.3f}, {:5.3f}]".
              format(beta[ind.X], beta[ind.Y]))
        print("    eta_x          = {:9.3e}".format(eta[ind.X]))
        print("    xi             = [{:5.3f}, {:5.3f}]".
              format(xi[ind.X], xi[ind.Y]))
        print("\n    phi_sp         = {:8.5f}".format(phi))
        print("    C [m]          = {:8.5f}".format(lat_prop.compute_circ()))
        print("\n    phi_rb_1       = {:8.5f}".format(phi_rb_1))
        print("    phi_bend       = {:8.5f}".format(phi_bend))
        lat_prop.prt_rad()
        prm_list.prt_prm(prm)

    def compute_chi_2(eta, nu_uc, xi):
        prt = not False

        dchi_2 = weight[0]*(lat_prop._eps[ind.X]-eps_x_des)**2
        chi_2 = dchi_2
        if prt:
            print("\n  dchi2(eps_x)    = {:9.3e}".format(dchi_2))

        dchi_2 = weight[1]*(nu_uc[ind.X]-nu_uc_des[ind.X])**2            
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_x)  = {:9.3e}".format(dchi_2))

        dchi_2 = weight[2]*(nu_uc[ind.Y]-nu_uc_des[ind.Y])**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(nu_uc_y)  = {:9.3e}".format(dchi_2))

        dchi_2 = weight[3]*eta[ind.X]**2
        chi_2 += dchi_2
        if prt:
            print("  dchi2(eta_x)    = {:9.3e}".format(dchi_2))

        dchi_2 = weight[4]*(xi[ind.X]**2+xi[ind.Y]**2)
        chi_2 += dchi_2
        if prt:
            print("  dchi2(xi)       = {:9.3e}".format(dchi_2))

        return chi_2

    def f_uc(prm):
        nonlocal chi_2_min, n_iter

        n_iter += 1
        prm_list.set_prm(prm)
        if sp_bend is not None:
            sp_bend.set_phi_lat()

        stable, eta, beta, nu, xi = compute_linear_optics(lat_prop, nld)
        print(eta)

        if stable:
            chi_2 = compute_chi_2(eta, nu, xi)
            if chi_2 < chi_2_min:
                prt_iter(prm, chi_2, eta, beta, nu, xi)
                pc.prt_lat(lat_prop, "opt_uc.txt", prm_list)
                chi_2_min = min(chi_2, chi_2_min)
            return chi_2
        else:
            return 1e30


    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4
    g_tol    = 1e-5

    prm, bounds = prm_list.get_prm()

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    if True:
        minimum = opt.minimize(
            f_uc,
            prm,
            method="Powell",
            # callback=prt_iter,
            bounds = bounds,
            options={"ftol": f_tol, "xtol": x_tol, "maxiter": max_iter}
        )
    else:
        minimum = opt.minimize(
            f_uc,
            prm,
            method="CG",
            # callback=prt_iter,
            jac=None,
            bounds = bounds,
            options={"gtol": g_tol, "maxiter": max_iter}
        )

    print("\n".join(minimum))


# TPSA max order.
gtpsa_prop.no = 2

cod_eps   = 1e-15
E_0       = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

lat_prop.prt_lat("max_4u_uc_lat.txt")

uc_centre = lat_prop._lattice.find("dsim", 0).index
print("\nunit cell centre {:5s} loc = {:d}".
      format(lat_prop._lattice[uc_centre].name, uc_centre))

np.set_printoptions(formatter={"float": "{:10.3e}".format})

b_3_list = ["s1", "s2", "s3", "s4"]
nld = nld_class.nonlin_dyn_class(lat_prop, A_max, beta_inj, delta_max, b_3_list)
nld.zero_mult(lat_prop, 3)
nld.zero_mult(lat_prop, 4)

if not False:
    compute_linear_optics(lat_prop, nld)

    lat_prop.prt_lat_param()
    lat_prop.prt_rad()
    lat_prop.prt_M()
    lat_prop.prt_M_rad()

# Weights.
weight = np.array([
    1e17, # eps_x.
    0e0,  # nu_uc_x.
    0e0,  # nu_uc_y.
    1e0,  # eta_x.
    1e-9, # xi.
])

d2_list = []
# d2_list = ["d2_0", "d2_1", "d2_2", "d2_3", "d2_4", "d2_5"]

# d2_bend = pc.bend_class(lat_prop, d2_list, phi_max, b_2_max)

if True:
    prms = [
        ("r1_f1", "b_2"),
        ("dsim",  "b_2"),
        ("s3_f1", "b_2"),

        ("dsim",  "phi"),
]

    dprm_list = np.array([
        1e-3, 1e-3, 1e-3,
        1e-3, 1e-3
    ])

    el_name = "r1_f1"
    n = len(lat_prop._lattice.elements_with_name(el_name))
    sp_bend = pc.phi_lat_class(lat_prop, n, el_name)
else:
    prms = [
        ("qf1",      "b_2"),

        ("d2_0",     "b_2"),
        ("d2_1",     "b_2"),
        ("d2_2",     "b_2"),
        ("d2_3",     "b_2"),
        ("d2_4",     "b_2"),
        ("d2_5",     "b_2"),

        ("d2_0",     "phi"),
        ("d2_1",     "phi"),
        ("d2_2",     "phi"),
        ("d2_3",     "phi"),
        ("d2_4",     "phi"),
        # ("d2_5",     "phi"),

        ("qf1",     "phi"),
    ]

prm_list = pc.prm_class(lat_prop, prms, b_2_max)

rb_1    = "r1_f1"
phi_b   = "dsim"

opt_uc(lat_prop, prm_list, uc_centre, weight, sp_bend, rb_1, phi_b, eps_x_des,
       nu_uc_des, nld)
