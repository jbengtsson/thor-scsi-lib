"""Use Case:
     Linear optics gymnastics.
"""

import os

import math
import numpy as np

import gtpsa
import thor_scsi.lib as tslib

from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
     radiate as rad

from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

tpsa_order = 2

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(7, tpsa_order)


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB",
                     "BESSY-III", "ipac_2023")
t_file = os.path.join(t_dir, "b3_cf425cf_thor_scsi.lat")

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

M = lo.compute_map(lat, model_state, desc=desc)
M_mat = M.jacobian()[:6, :6]
stable, A_mat, A_inv_mat, alpha_rad = lo.compute_M_diag(n_dof, M_mat)

print("\nM:\n", mat2txt(M_mat))
print("\nA:\n", mat2txt(A_mat))

stable, U_0, J, tau, eps = \
    rad.compute_radiation(lat, model_state, 2.5e9, 1e-15, desc=desc)

