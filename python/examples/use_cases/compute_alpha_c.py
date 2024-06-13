"""Use Case:
     Compute the:

       momentum compaction, alpha_c, to arbitrary order

       fixed points.

     Compute & graph the:

       linear optics

       2nd order dispersion.

       longitudinal phase-space portrait

"""


import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="ERROR")
logger = logging.getLogger("thor_scsi")


import math
import os
import sys

import copy

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import lattice_properties as lp, linear_optics as lo, \
    courant_snyder as cs, closed_orbit as co, index_class as ind

from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt


ind = ind.index_class()


def plot_D(lat_prop, s, disp, file_name):
    # Turn interactive mode off.
    plt.ioff()

    fig, gr_1 = plt.subplots(1)

    gr_1.set_title("Dispersion")
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel("")

    gr_1.plot(s, disp[1, :, ind.x], "g", label=r"$\eta_x$")
    gr_1.plot(s, disp[1, :, ind.px], "r", label=r"$\eta'_x$")
    gr_1.plot(s, disp[2, :, ind.x], "b", label=r"$\eta^{(2)}_x$")
    gr_1.legend()

    gr_1_r = gr_1.twinx()
    gr_1_r.set_ylim([-2.0, 20.0])
    gr_1_r.set_yticks([])
    gr_1_r.step(s, lat_prop._type_code, "k")

    fig.tight_layout()

    plt.savefig(file_name)
    print("\nPlot saved as:", file_name)

    plt.show()


def plot_D_alpha_1(lat_prop, s, disp, file_name):
    D_over_rho = compute_D_over_rho(lat_prop, disp)

    # Turn interactive mode off.
    plt.ioff()

    fig, gr_1 = plt.subplots(1)

    gr_1.set_title(r"Contribution to $\alpha^{(1)}_{\mathrm {c}}$")
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel(r"$\eta^{(1)}_x/\rho$")
    gr_1.step(s, D_over_rho, "b")

    gr_1_r = gr_1.twinx()
    gr_1_r.set_ylim([-2.0, 20.0])
    gr_1_r.set_yticks([])
    gr_1_r.step(s, lat_prop._type_code, "k")

    fig.tight_layout()

    plt.savefig(file_name)
    print("\nPlot saved as:", file_name)

    plt.show()


def plot_D_alpha_2(lat_prop, s, disp, file_name):
    D_der_sqr_over_2 = compute_D_der_sqr_over_2(lat_prop, disp)
    D_2_over_rho = compute_D_2_over_rho(lat_prop._lattice, disp)

    # Turn interactive mode off.
    plt.ioff()

    fig, gr_1 = plt.subplots(1)

    gr_1.set_title(r"Contribution to $\alpha^{(2)}_{\mathrm {c}}$")
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel("")
    gr_1.plot(s, D_der_sqr_over_2, "r", label=r"$(\eta'^{(1)})^2_x/2$")
    gr_1.legend(loc="upper left")

    y_min, y_max = gr_1.get_ylim();
    y_scl = (y_max-y_min)/20.0
    gr_1.step(s, y_scl*(lat_prop._type_code), "k")

    gr_1_r = gr_1.twinx()
    gr_1_r.plot(s, D_2_over_rho, "b", label=r"$\eta^{(2)}/\rho_x$")
    gr_1_r.legend(loc="upper right")

    fig.tight_layout()

    plt.savefig(file_name)
    print("\nPlot saved as:", file_name)

    plt.show()


def plot_H_long(phi, delta, H, file_name, title):
    # Turn interactive mode off.
    plt.ioff()

    fig, gr_1 = plt.subplots(1)

    gr_1.set_title(title)
    gr_1.set_xlabel(r"phi [$^\circ$]")
    gr_1.set_ylabel(r"$\delta$ [%]")
    gr_1.contour(phi, 1e2*delta, H, 30)

    fig.tight_layout()

    plt.savefig(file_name)
    print("\nPlot saved as:", file_name)

    plt.show()


def print_Twiss_param(str, Twiss):
    # eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    # that way I also check that Twiss has exactly three parameters
    eta, alpha, beta = Twiss
    print(str, end="")
    print(f"  eta    = [{eta[ind.X]:9.3e}, {eta[ind.Y]:9.3e}]")
    print(f"  alpha  = [{alpha[ind.X]:9.3e}, {alpha[ind.Y]:9.3e}]")
    print(f"  beta   = [{beta[ind.X]:5.3f}, {beta[ind.Y]:5.3f}]")


def compute_alpha_c(map):
    if not True:
        print("\nmap[ct]:\n")
        map.ct.print()

    C = lat_prop.compute_circ()
    print(f"\nC [m] = {C:5.3f}")
    index = np.zeros(nv, int)
    alpha_c = np.zeros(no+1)

    print("\nalpha_c:")
    for k in range(1, no+1):
        index[ind.delta] = k
        alpha_c[k] = map.ct.get(index)/C
        print("  {:1d} {:10.3e}".format(k, alpha_c[k]))

    n_alpha_c = len(alpha_c)
    if n_alpha_c >= 3:
        print("\nFixed points to O(3) [%]: {:5.2f}, {:5.2f}".
              format(0e0, -1e2*alpha_c[1]/alpha_c[2]))

    if n_alpha_c >= 4:
        po2 = alpha_c[2]/(2e0*alpha_c[3])
        q = alpha_c[1]/alpha_c[3]
        pm = np.sqrt(po2**2-q)
        print("Fixed points to O(4) [%]: {:5.2f}, {:5.2f}, {:5.2f}".
              format(0e0, -1e2*(po2+pm), -1e2*(po2-pm)))

    return alpha_c


def compute_D(map):
    A0 = gtpsa.ss_vect_tpsa(desc, no)
    # Compute canonical transformation to delta-dependent fixed point.
    ts.GoFix(map, A0)
    index = np.zeros(nv, int)
    print("\ndispersion:")
    for k in range(1, no+1):
        index[ind.delta] = k
        print("  {:1d} {:10.3e} {:10.3e}".
              format(k, A0.x.get(index), A0.px.get(index)))
    return A0


def compute_D_along_lattice(lat_prop, A0):
    A = A0
    ind1 = np.zeros(nv, int)
    ind1[ind.delta] = 1
    ind2 = np.zeros(nv, int)
    ind2[ind.delta] = 2
    n = len(lat_prop._lattice)
    s = np.zeros(n)
    disp = np.zeros((3, n, 2))
    s[0] = 0e0
    disp[1, 0, ind.x] = A0.x.get(ind1)
    disp[1, 0, ind.px] = A0.px.get(ind1)
    disp[2, 0, ind.x] = A0.x.get(ind2)
    disp[2, 0, ind.px] = A0.px.get(ind2)
    for k in range(n-1):
        lat_prop._lattice.propagate(lat_prop._model_state, A0, k+1, 1)
        s[k+1] = s[k]+ lat_prop._lattice[k+1].get_length()
        disp[1, k+1, ind.x] = A0.x.get(ind1)
        disp[1, k+1, ind.px] = A0.px.get(ind1)
        disp[2, k+1, ind.x] = A0.x.get(ind2)
        disp[2, k+1, ind.px] = A0.px.get(ind2)
    return s, disp


def print_D_along_lattice(lat, s, disp):
    print("\n                      eta_0,x   eta_0,p_x   eta_1,x   eta_1,p_x")
    for k in range(len(disp[1])):
        print("{:3d} {:6.3f} {:8s} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".
              format(k+1, s[k], lat[k].name, disp[1, k, ind.x],
                     disp[1, k, pind.x], disp[2, k, ind.x], disp[2, k, ind.px]))


def compute_D_over_rho(lat_prop, disp):
    n = len(lat_prop._lattice)
    D_over_rho = np.zeros(n)
    for k in range(n):
        # Compute average across element.
        if k > 0:
            D_avg = (disp[1, k, ind.x]+disp[1, k-1, ind.x])/2e0
        else:
            D_avg = disp[1, k, ind.x]
        if type(lat_prop._lattice[k]) == ts.Bending:
            rho_inv = lat_prop._lattice[k].get_curvature()
            D_over_rho[k] = D_avg*rho_inv
            # print("  {:3d} {:10s} {:5.3f} {:7.3f} {:10.3e} {:10.3e} {:10.3e}".
            #       format(k, lat[k].name, lat[k].get_length(), rho_inv,
            #              disp[1, k, ind.x], D_avg, D_over_rho[k]))
        else:
            D_over_rho[k] = 0e0
    return D_over_rho


def compute_D_der_sqr_over_2(lat_prop, disp):
    n = len(lat_prop._lattice)
    D_der_sqr_over_2 = np.zeros(n)
    for k in range(n):
        # Compute average across element.
        if k > 0:
            D_der_avg = (disp[1, k, ind.px]+disp[1, k-1, ind.px])/2e0
        else:
            D_der_avg = disp[1, k, ind.px]
        D_der_sqr_over_2[k] = D_der_avg**2/2e0
    return D_der_sqr_over_2


def compute_D_2_over_rho(lat, disp):
    n = len(lat)
    D_2_over_rho = np.zeros(n)
    for k in range(n):
        # Compute average across element.
        if k > 0:
            D_2_avg = (disp[2, k, ind.x]+disp[2, k-1, ind.x])/2e0
        else:
            D_2_avg = disp[1, k, ind.x]
        if type(lat[k]) == ts.Bending:
            rho_inv = lat[k].get_curvature()
            D_2_over_rho[k] = D_2_avg*rho_inv
            # print("  {:3d} {:10s} {:5.3f} {:7.3f} {:10.3e} {:10.3e} {:10.3e}".
            #       format(k, lat[k].name, lat[k].get_length(), rho_inv,
            #              disp[1, k, ind.x], D_2_avg, D_2_over_rho[k]))
        else:
            D_2_over_rho[k] = 0e0
    return D_2_over_rho


def print_map(str, map):
    n_dof = 3
    print(str)
    for k in range(2*n_dof):
        map.iloc[k].print()


def H_long(E0, phi, delta, h_rf, V_rf, phi0, alpha_c):
    H = V_rf/E0*(np.cos(phi+phi0)+phi*np.sin(phi0))
    for k in range(1, len(alpha_c)):
        H += 2e0*np.pi*h_rf*alpha_c[k]*delta**(k+1)/(k+1)
    return H


def compute_H_long(
        lat_prop, E0, alpha_c, n, phi_max, delta_max, U_0, neg_alpha_c):
    cav = lat_prop._lattice.find("cav", 0)
    h_rf = cav.get_harmonic_number()
    V_rf = cav.get_voltage()
    f_rf = cav.get_frequency()

    phi0 = - abs(np.arcsin(U_0/V_rf))
    if neg_alpha_c:
        phi0 += pi

    delta_rf = \
        np.sqrt(-V_rf*np.cos(np.pi+phi0)
                *(2e0-(np.pi-2e0*(np.pi+phi0))*np.tan(np.pi+phi0))
                /(alpha_c[1]*np.pi*h_rf*E0))

    print("\nh_rf                 = {:1d}".format(h_rf))
    print("V_rf [MV]            = {:3.1f}".format(1e-6*V_rf))
    print("f_rf [MHz]           = {:3.1f}".format(1e-6*f_rf))
    print("U_0 [keV]            = {:3.1f}".format(1e-3*U_0))

    if not neg_alpha_c:
        print("phi0 [deg]           = {:4.2f}".
              format(abs(phi0)*180e0/np.pi-180e0))
    else:
        print("phi0 [deg]           = 180 - {:4.2f}".
              format(abs(phi0)*180e0/np.py-180e0))

    print("RF bucket height [%] = {:3.1f}".format(1e2*delta_rf))

    phi = np.zeros(2*n+1)
    delta = np.zeros(2*n+1)
    H = np.zeros((2*n+1, 2*n+1))
    for i in range(-n, n+1):
        for j in range(-n, n+1):
            phi1 = i*phi_max*np.pi/(180e0*n)
            delta[j+n] = j*delta_max/n
            H[j+n, i+n] = \
                H_long(E0, phi1, delta[j+n], h_rf, V_rf, np.pi+phi0, alpha_c)
            phi[i+n] = (phi0+phi1)*180e0/np.pi

    return phi, delta, H


def print_H_long(file_name, phi, delta, H):
    with open(file_name, "w") as sys.stdout:
        for i in range(len(phi)):
            for j in range(len(delta)):
                print(" {:7.2f} {:7.3f} {:12.5e}".
                  format(phi[i], 1e2*delta[j], H[j, i]))
            print()


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 3
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 2.5e9
U_0     = 22.4e3

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

home_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB",
                     "BESSY-III", "ipac_2023")
lat_name = "b3_cf425cf_thor_scsi"
file_name = os.path.join(home_dir, lat_name+".lat")

print("\nlattice file:  \n", file_name)

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

print("\nTotal bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))
print("Circumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))

try:
    # Compute Twiss parameters along lattice.
    if not lat_prop.comp_per_sol():
        print("\ncomp_per_sol - unstable")
        raise ValueError

    # Compute radiation properties.
    if not lat_prop.compute_radiation():
        print("\ncompute_radiation - unstable")
        raise ValueError
except ValueError:
    exit
else:
    lat_prop.prt_lat_param()
    lat_prop.prt_rad()
    lat_prop.prt_M()
    lat_prop.prt_M_rad()
    lat_prop.prt_Twiss(lat_name+"_Twiss.txt")

if not True:
    print_Twiss(lat, data)

if True:
    lat_prop.plt_Twiss( "lin_opt.png", not False)

r = co.compute_closed_orbit(
    lat_prop._lattice, lat_prop._model_state, delta=0e0, eps=1e-10, desc=desc)

# Compute the Taylor map.

map = gtpsa.ss_vect_tpsa(desc, no)
map.set_identity()
r.x0.ct = 0e0
map += r.x0
lat_prop._lattice.propagate(lat_prop._model_state, map)
print("\nmap:\n" + mat2txt(map.jacobian()[:6, :6]))

alpha_c = compute_alpha_c(map)
A0 = compute_D(map)
s, disp = compute_D_along_lattice(lat_prop, A0)

if not True:
    print_D_along_lattice(lat_prop._lattice, s, disp, lat_prop._types)
    
if True:
    plot_D(lat_prop, s, disp, "D.png")
    plot_D_alpha_1(lat_prop, s, disp, "D_alpha.png")
    plot_D_alpha_2(lat_prop, s, disp, "D_alpha.png")

phi, delta, H = \
    compute_H_long(lat_prop, E_0, alpha_c, 10, 180e0, 20e-2, U_0, False)

if not True:
    print_H_long("H_long.dat", phi, delta, H)

if True:
    plot_H_long(phi, delta, H, "H_long.png", r"$H_\parallel$")
