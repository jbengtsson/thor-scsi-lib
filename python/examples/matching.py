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


import os

import numpy as np
import matplotlib.pyplot as plt

import gtpsa

import thor_scsi.lib as tslib
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.courant_snyder import compute_A_CS, compute_A
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv
from thor_scsi.utils.linear_optics import propagate_and_find_fixed_point, \
    compute_twiss_along_lattice, compute_dispersion
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

    
def compute_periodic_solution(acc, model_state, A):
    # Compute the periodic solution for the super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    ds = compute_twiss_along_lattice(n_dof, acc, model_state, A=A, desc=desc)

    return ds


def propagate_Twiss_half_cell(loc, A, acc, model_state):
    n_dof = 2
    A0 = A
    acc.propagate(model_state, A0, loc, len(acc)-loc)
    A0.set_jacobian(compute_A_CS(n_dof, A0.jacobian())[0])
    logger.info("\npropagate_Twiss_half_cell\nA:\n" + prt2txt(A0))


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

# Read in & parse lattice file.
acc = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# Degrees of freedom - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

M, A = propagate_and_find_fixed_point(n_dof, acc, model_state, desc=desc)
print("\mM:\n", mat2txt(M))
# Compute Twiss parameters along lattice.
data = compute_periodic_solution(acc, model_state, A)

# Turn interactive mode off.
plt.ioff()
plt_Twiss(data)
# plt.show()

# Get dispersion & Twiss parameters at lattice midpoint.
loc = acc.find("b_t", 1).index
print(f"\nlattice midpoint: {loc:1d} {acc[loc].name:s}")
# loc = 0
eta = [data.dispersion.values[loc, 4], data.dispersion.values[loc, 2]]
alpha = [data.twiss.values[loc, X_, 0], data.twiss.values[loc, Y_, 0]]
beta = [data.twiss.values[loc, X_, 1], data.twiss.values[loc, Y_, 1]]
print(f"  eta   = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
print(f"  alpha = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
print(f"  beta  = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")

A = gtpsa.ss_vect_tpsa(desc, 1)
A.set_zero()
A.set_jacobian(compute_A(eta, alpha, beta))
propagate_Twiss_half_cell(loc, A, acc, model_state)
