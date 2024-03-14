"""Use Case:
     Parametric scans of the unit cell.
"""


import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
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

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
     radiate as rad

# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

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


def compute_periodic_solution(lat, model_state, named_index, desc, prt):
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
    if stable:
        if prt:
            print("\nM:\n", mat2txt(M.jacobian()))
        res= cs.compute_Twiss_A(A)
        Twiss = res[:3]
        if prt:
             print_Twiss_params("\nTwiss:\n", Twiss)
        A_map = gtpsa.ss_vect_tpsa(desc, no)
        A_map.set_jacobian(A)
        ds = \
            lo.compute_Twiss_along_lattice(
                n_dof, lat, model_state, A=A_map, desc=desc,
                mapping=named_index)
    else:
        ds = np.nan
        print("\ncompute_periodic_solution: unstable")

    return M, A, ds

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB")
# t_file = os.path.join(t_dir, 'b3_tst.lat')
# t_file = os.path.join(t_dir, "b3_sf_40Grad_JB.lat")
# t_file = os.path.join(t_dir, "b3_sfsf_4Quads_unitcell.lat")
t_file = os.path.join(t_dir, "b3_sfsf4Q_tracy_jb.lat")
# t_file = os.path.join(t_dir, "b3_sfsf4Q_tracy_jb_2.lat")

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

# Compute Twiss parameters along lattice.
M, A, data = \
    compute_periodic_solution(
        lat, model_state, named_index, desc, True)

plot_Twiss(data, "before.png", "Linear Optics - Before")
if False:
    print_Twiss(lat, data)
if False:
    plt.show()

