"""Read lattice file and calculate radiation
"""
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils import linear_optics
from thor_scsi.utils.output import vec2txt, mat2txt
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv, df_header
from thor_scsi.utils.linear_optics import compute_twiss_along_lattice
from thor_scsi.utils.phase_space_vector import ps_jac2ss_vect_tps, map2numpy
from thor_scsi.utils.phase_advance import compute_nus_for_symplectic_matrix, compute_nus
from thor_scsi.utils.courant_snyder import compute_dispersion, compute_twiss_M, compute_twiss_A, compute_twiss_A_A_tp

import thor_scsi.lib as tslib

import os
import sys
import xarray as xr
import pandas as pd
import numpy as np
import copy
from typing import Sequence


# X_, Y_, Z_ = spatial_ind.X_, spatial_ind.Y_, spatial_ind.Z_

# x_ = phase_space_ind.x_
# px_ = phase_space_ind.px_
# y_ = phase_space_ind.y_
# py_ = phase_space_ind.py_
# delta_ = phase_space_ind.delta_
# ct_ = phase_space_ind.ct_
#
#
# def prt_np_vec(name, vec):
#     print(name, end="")
#     print(vec2txt(vec))
#
#
# def prt_np_cmplx_vec(name, vec):
#     prt_np_vec(name, vec)
#
#
# def prt_np_mat(name, mat):
#     print(name, end="")
#     print(mat2txt(mat))
#
#
# def prt_np_cmplx_mat(name, mat):
#     prt_np_mat(name, mat)
#
#
# def get_mat(t_map):
#     mat = tslib.ss_vect_tps_to_mat(t_map)
#     mat = np.array(mat)
#     return mat
#
#
# def get_map(M):
#     [Id, t_map] = [tslib.ss_vect_tps(), tslib.ss_vect_tps()]
#     Id.set_identity()
#     t_map.set_zero()
#     for j in range(len(M)):
#         for k in range(len(M[0])):
#             t_map[j] += M[j][k] * Id[k]
#     return t_map





t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)
calc_config = tslib.ConfigType()


ds = compute_twiss_along_lattice(acc, calc_config)
df = twiss_ds_to_df(ds)
# print(df)

with open("twiss.tsf", "wt") as fp:
    # fp.write(df_header())
    # fp.write("\n")
    fp.write(df_to_tsv(df))

df.to_json("twiss.json")
