"""Use Case:
     Minimise the driving terms.
     See SLS Tech Note 09/97

     https://ados.web.psi.ch/slsnotes/sls0997.pdf
"""

"""Use Case:
     Super period straight section matching with a triplet.
   Constraints:
    [alpha_x,y, beta_x,y, eta_x, eta'_x] at the centre of the straigth.
   Parameters:
     triplet gradients (3) & spacing (3).
"""
import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="INFO")
logger = logging.getLogger("thor_scsi")


from scipy import optimize

import copy as _copy
from dataclasses import dataclass
import math
import os
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as tslib

from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs

# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

tpsa_order = 3

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(6, tpsa_order)

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


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole = 3


def plt_Twiss(ds, file_name, title):
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


def prt_Twiss(str, Twiss):
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


def compute_periodic_solution(lat, model_state):
    """
    Todo:
        model_state: rename to calculation_configuration or calc_config
    """
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    M, A = lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc)
    print("\nM:\n", mat2txt(M.jacobian()))
    res= cs.compute_Twiss_A(A)
    Twiss = res[:3]
    prt_Twiss("\nTwiss:\n", Twiss)
    print("\nend of twiss\n")
    ds = lo.compute_Twiss_along_lattice(n_dof, lat, model_state, A=A, desc=desc)

    return M, A, ds


def prt_twiss_sxt(lat, data, fam_name):
    print("\n  name        s    beta_x beta_y")
    for k in range(len(lat)):
        if (lat[k].name == fam_name):
            print("  {:8s} {:6.3f} {:6.3f} {:6.3f}".
                  format(lat[k].name, data.s.values[k],
                         data.twiss.values[k, 0, 1],
                         data.twiss.values[k, 1, 1]))

    
lat_file = "b3_sfsf_thor_sextsplit.lat"
t_dir = os.path.join("/Volumes", "qst", "thor_scsi", "JB")
t_file = os.path.join(t_dir, lat_file)

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

# Compute Twiss parameters along lattice.
M, A, data = compute_periodic_solution(lat, model_state)

if True:
    plt_Twiss(data, "lin_opt.png", "Linear Optics")
    plt.show()
    print("\nPlots saved as: before.png & after.png")

print("\nR:\n", mat2txt(np.linalg.inv(A) @ M.jacobian() @ A))

prt_twiss_sxt(lat, data, "om_sf")
prt_twiss_sxt(lat, data, "om_sd")
