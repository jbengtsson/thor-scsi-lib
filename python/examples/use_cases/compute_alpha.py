"""Use Case:
     Compute the momentum compaction, alpha_c, to arbitrary order.
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

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
    radiate as rad, closed_orbit as co

# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

# Configuration space coordinates.
X_, Y_, Z_ = [
    ts.spatial_index.X,
    ts.spatial_index.Y,
    ts.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    ts.phase_space_index_internal.x,
    ts.phase_space_index_internal.px,
    ts.phase_space_index_internal.y,
    ts.phase_space_index_internal.py,
    ts.phase_space_index_internal.ct,
    ts.phase_space_index_internal.delta,
]


def read_lattice(t_file):
    # Read in & parse lattice file.
    lat = accelerator_from_config(t_file)

    # Set lattice state (Rf cavity on/off, etc.)
    model_state = ts.ConfigType()

    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    return n_dof, lat, model_state


def plot_Twiss(ds, file_name, title):
    # Turn interactive mode off.
    plt.ioff()

    fig, (gr_1, gr_2) = plt.subplots(2)

    gr_1.set_title(title)
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel(r"$\beta_{x,y}$ [m]")
    gr_1.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), "b-",
              label=r"$\beta_x$")
    gr_1.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), "r-",
              label=r"$\beta_y$")
    gr_1.legend()

    gr_2.set_xlabel("s [m]")
    gr_2.set_ylabel(r"$\eta_x$ [m]")
    gr_2.plot(ds.s, ds.dispersion.sel(phase_coordinate="x"), label=r"$\eta_x$")
    fig.tight_layout()
    plt.savefig(file_name)


def print_Twiss(str, Twiss):
    """

    todo:
        rename str e.g. header? prefix?

    """
    # eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    # that way I also check that Twiss has exactly three parameters
    eta, alpha, beta = Twiss
    print(str, end="")
    print(f"  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
    print(f"  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
    print(f"  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")


def compute_periodic_solution(lat, model_state, named_index, desc):
    """
    Todo:
        model_state: rename to calculation_configuration or calc_config
    """
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    stable, M, A = lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc)
    print("\nM:", M)
    res = cs.compute_Twiss_A(A)
    Twiss = res[:3]
    print_Twiss("Twiss:\n", Twiss)
    A_map = gtpsa.ss_vect_tpsa(desc, no)
    A_map.set_jacobian(A)
    ds = \
        lo.compute_Twiss_along_lattice(
            n_dof, lat, model_state, A=A_map, desc=desc, mapping=named_index)

    return M, A, ds


def compute_alpha_c(map):
    if not True:
        print("\nmap[ct]:\n")
        map.ct.print()

    C = rad.compute_circ(lat)
    print(f"\nC [m] = {C:5.3f}")
    ind = np.zeros(nv, int)
    alpha_c = np.zeros(no+1)

    print("\nalpha_c:")
    for k in range(1, no+1):
        ind[delta_] = k
        alpha_c[k] = map.ct.get(ind)/C
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


def compute_disp(map):
    A0 = gtpsa.ss_vect_tpsa(desc, no)
    # Compute canonical transformation to delta-dependent fixed point.
    ts.GoFix(map, A0)
    ind = np.zeros(nv, int)
    print("\ndispersion:")
    for k in range(1, no+1):
        ind[delta_] = k
        print("  {:1d} {:10.3e} {:10.3e}".
              format(k, A0.x.get(ind), A0.px.get(ind)))
    return A0


def print_disp_along_lattice(lat, s, disp):
    print("\n                      eta_0,x   eta_0,p_x   eta_1,x   eta_1,p_x")
    for k in range(len(disp[1])):
        print("{:3d} {:6.3f} {:8s} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".
              format(k+1, s[k], lat[k].name, disp[1, k, x_], disp[1, k, px_],
                     disp[2, k, x_], disp[2, k, px_]))


def compute_disp_along_lattice(lat, model_state, A0):
    A = A0
    ind1 = np.zeros(nv, int)
    ind1[delta_] = 1
    ind2 = np.zeros(nv, int)
    ind2[delta_] = 2
    n = len(lat)
    s = np.zeros(n)
    disp = np.zeros((3, n, 2))
    s[0] = 0e0
    disp[1, 0, x_] = A0.x.get(ind1)
    disp[1, 0, px_] = A0.px.get(ind1)
    disp[2, 0, x_] = A0.x.get(ind2)
    disp[2, 0, px_] = A0.px.get(ind2)
    for k in range(n-1):
        lat.propagate(model_state, A0, k, 1)
        s[k+1] = s[k]+ lat[k].get_length()
        disp[1, k+1, x_] = A0.x.get(ind1)
        disp[1, k+1, px_] = A0.px.get(ind1)
        disp[2, k+1, x_] = A0.x.get(ind2)
        disp[2, k+1, px_] = A0.px.get(ind2)
    return s, disp


def plot_disp(s, disp, file_name, title):
    # Turn interactive mode off.
    plt.ioff()

    fig, gr_1 = plt.subplots(1)

    gr_1.set_title(title)
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel(r"$\eta_{x}$ [m]")
    gr_1.plot(s, disp[2, :, x_], "b-")
    gr_1.legend()

    fig.tight_layout()
    plt.savefig(file_name)


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


def compute_H_long(lat, E0, alpha_c, n, phi_max, delta_max, U0, neg_alpha_c):
    cav = lat.find("cav", 0)
    h_rf = cav.get_harmonic_number()
    V_rf = cav.get_voltage()
    f_rf = cav.get_frequency()

    phi0 = - abs(np.arcsin(U0/V_rf))
    if neg_alpha_c:
        phi0 += pi

    delta_rf = \
        np.sqrt(-V_rf*np.cos(np.pi+phi0)
                *(2e0-(np.pi-2e0*(np.pi+phi0))*np.tan(np.pi+phi0))
                /(alpha_c[1]*np.pi*h_rf*E0))

    print("\nh_rf                 = {:1d}".format(h_rf))
    print("V_rf [MV]            = {:3.1f}".format(1e-6*V_rf))
    print("f_rf [MHz]           = {:3.1f}".format(1e-6*f_rf))
    print("U0 [keV]             = {:3.1f}".format(1e-3*U0))

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


def plot_H_long(phi, delta, H, file_name, title):
    # Turn interactive mode off.
    plt.ioff()

    fig, gr_1 = plt.subplots(1)

    gr_1.set_title(title)
    gr_1.set_xlabel("phi [$^\circ$]")
    gr_1.set_ylabel(r"$\delta$ [%]")
    gr_1.contour(phi, 1e2*delta, H, 30)
    gr_1.legend()

    fig.tight_layout()
    plt.savefig(file_name)


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 3
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

E0 = 2.5e9
U0 = 22.4e3

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB",
                     "BESSY-III", "ipac_2023")
t_file = os.path.join(t_dir, "b3_cf425cf_thor_scsi.lat")

print("\nlattice file:  \n", t_file)

n_dof, lat, model_state = read_lattice(t_file)

M, A, data = \
    compute_periodic_solution(lat, model_state, named_index, desc)

if True:
    plot_Twiss(data, "lin_opt.png", "Linear Optics")
    plt.show()
    print("\nPlot saved as: lin_opt.png")

r = co.compute_closed_orbit(lat, model_state, delta=0e0, eps=1e-10, desc=desc)

# Compute the Taylor map.

map = gtpsa.ss_vect_tpsa(desc, no)
map.set_identity()
r.x0.ct = 0e0
map += r.x0
lat.propagate(model_state, map)
print("\nmap:", map)

alpha_c = compute_alpha_c(map)
A0 = compute_disp(map)
s, disp = compute_disp_along_lattice(lat, model_state, A0)

if not True:
    print_disp_along_lattice(lat, s, disp)
    
if True:
    plot_disp(s, disp, "disp_2.png", "2nd Order Dispersion")
    plt.show()
    print("\nPlot saved as: disp_2.png")

phi, delta, H = compute_H_long(lat, E0, alpha_c, 10, 180e0, 20e-2, U0, False)

if not True:
    print_H_long("H_long.dat", phi, delta, H)

if True:
    plot_H_long(phi, delta, H, "H_long.png", "H_Long")
    plt.show()
    print("\nPlot saved as: H_long.png")
