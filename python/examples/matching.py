"""Straight section matching for a super period with a triplet.
   Constraints:
    [alpha_x,y, beta_x,y, eta_x, eta'_x] at the centre of the straigth.
   Parameters:
     triplet gradients (3) & spacing (3).
"""

import logging
# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="INFO")
logger = logging.getLogger("thor_scsi")


from scipy import optimize

import os
import copy

import numpy as np
import matplotlib.pyplot as plt

import gtpsa

import thor_scsi.lib as tslib
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv
from thor_scsi.utils.linear_optics import propagate_and_find_fixed_point, \
    compute_Twiss_along_lattice, compute_dispersion
from thor_scsi.utils.courant_snyder import compute_A_CS, compute_A, compute_Twiss_A
from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(6, 1)

# Configuration space coordinates.
X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta,
]


def plt_Twiss(ds):
    fig, (gr_1, gr_2) = plt.subplots(2)

    gr_1.set_title("Linear Optics")
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel(r"$\beta_{x,y}$ [m]")
    gr_1.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), label=r"$\beta_x$")
    gr_1.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), label=r"$\beta_y$")
    gr_1.legend()

    gr_2.set_xlabel("s [m]")
    gr_2.set_ylabel(r"$\eta_x$ [m]")
    gr_2.plot(ds.s, ds.dispersion.sel(phase_coordinate="x"), label=r"$\eta_x$")
    fig.tight_layout()

    
def compute_periodic_solution(lat, model_state, A):
    # Compute the periodic solution for the super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    ds = compute_Twiss_along_lattice(n_dof, lat, model_state, A=A, desc=desc)

    return ds


def prt_Twiss(str, Twiss):
    eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    print(str, end='')
    print(f"  eta   = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
    print(f"  alpha = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
    print(f"  beta  = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")


def set_b_2_fam(lat, fam_name, b_2):
    for q in lat.elements_with_name(fam_name):
        q.get_multipoles().set_multipole(2, b_2)


def set_L_fam(lat, fam_name, L):
    for q in lat.elements_with_name(fam_name):
        q.set_length(L)


# Global variables.
n_iter    = 0
chi_2_min = 1e30


def match_straight(lat, loc, prm_list, bounds, Twiss0, Twiss1, Twiss1_design,
                   weight):

    def get_prm(prm_list):
        prms = []
        for k in range(len(prm_list)):
            q = lat.find(prm_list[k][0], 0)
            if prm_list[k][1] == 'b_2':
                prms.append(q.get_multipoles().get_multipole(2).real)
            elif prm_list[k][1] == 'L':
                prms.append(q.get_length())
        return prms

    def f_match(prms):
        global n_iter
        global chi_2_min

        n_iter += 1
        for k in range(3):
            if prm_list[k][1] == 'b_2':
                set_b_2_fam(lat, prm_list[k][0], prms[k])
            else:
                set_L_fam(lat, prm_list[k][0], prms[k])

        # Local copy.
        A1 = copy.copy(A0)
        lat.propagate(model_state, A1, loc, len(lat)-loc)
        Twiss_k = compute_Twiss_A(A1.jacobian())[0:3]

        chi_2 = 0e0
        for k in range(len(Twiss_k)):
            chi_2 += weight[k]*np.sum((Twiss_k[k]-Twiss1_design[k])**2)

        if chi_2 < chi_2_min:
            chi_2_min = chi_2
            print(f"\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n  prms ="
                  + vec2txt(prms))
            prt_Twiss("\nTwiss_k:\n", Twiss_k)

        return chi_2


    max_iter = 100
    f_tol    = 1e-3
    x_tol    = 1e-3

    print("\nloc = ", loc)
    prt_Twiss("\nTwiss parameters at centre of super period:\n", Twiss0)
    prt_Twiss("\nTwiss parameters at centre of straight:\n",
              Twiss1)
    prt_Twiss("\nDesired Twiss parameters at the centre of straight:\n",
              Twiss1_design)

    A0 = gtpsa.ss_vect_tpsa(desc, 1)
    A0.set_zero()
    A0.set_jacobian(compute_A(Twiss0))

    # Initialise parameters.
    prms1 = prms0 = get_prm(prm_list)
    print("\nprms = ", prms1)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.
    minimum = \
        optimize.minimize(
            f_match, prms1, method="Powell", bounds=bounds,
            options={'xtol':x_tol, 'ftol':f_tol, 'maxiter':max_iter})

    print("\nInitial parameter values: ")
    f_match(prms0)
    print("\nFinal parameter values: ")
    f_match(minimum["x"])
    print("\n Minimum:\n", minimum)


def get_Twiss(loc):
    eta = \
        np.array([data.dispersion.sel(phase_coordinate="x").values[loc],
                  data.dispersion.sel(phase_coordinate="px").values[loc],
                  data.dispersion.sel(phase_coordinate="y").values[loc],
                  data.dispersion.sel(phase_coordinate="py").values[loc]])
    alpha = \
        np.array([data.twiss.sel(plane="x", par="alpha").values[loc],
                  data.twiss.sel(plane="y", par="alpha").values[loc]])
    beta = \
        np.array([data.twiss.sel(plane="x", par="beta").values[loc],
                  data.twiss.sel(plane="y", par="beta").values[loc]])
    return eta, alpha, beta


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

M, A = propagate_and_find_fixed_point(n_dof, lat, model_state, desc=desc)
print("\mM:\n", mat2txt(M))
Twiss = compute_Twiss_A(A)
prt_Twiss("\nTwiss:\n", Twiss)

# Compute Twiss parameters along lattice.
data = compute_periodic_solution(lat, model_state, A)

# Turn interactive mode off.
plt.ioff()
plt_Twiss(data)
# plt.show()

# Get dispersion & Twiss parameters at centre of super period.
loc = lat.find("b_t", 1).index
print(f"\nCentre of super period: {loc:1d} {lat[loc].name:s}")
Twiss0 = get_Twiss(loc)

# Twiss parameters at centre of straight.
Twiss1 = get_Twiss(len(lat)-1)

# Desired Twiss parameters at centre of straight.
eta   = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta  = np.array([2.5, 2.5])
Twiss1_design = eta, alpha, beta

weight = np.array([1e0, 1e3, 1e0])

# Triplet parameter family names & type.
prm_list = [('uq1', 'b_2'), ('uq2', 'b_2'), ('uq3', 'b_2'),
            ('ul1', 'L'), ('ul2', 'L'), ('ul3', 'L')]

# Max parameter range.
b_2_max = 10.0
L_min   = 0.1

bounds = [(-b_2_max, 0.0), (0.0, b_2_max), (-b_2_max, 0.0),
          (L_min, 0.15), (L_min, 0.3), (L_min, 0.35)]

match_straight(lat, loc, prm_list, bounds, Twiss0, Twiss1, Twiss1_design,
               weight)
