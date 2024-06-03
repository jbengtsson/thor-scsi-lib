'''
Author:

  Johan Bengtsson
  03/06/24

Tbt_bpm class to analyse turn-by-turn BPM data to extract the linear optics -
beta functions & phase advance - at the BPMs.
'''


import os
import numpy as np

from thor_scsi.utils.tbt_bpm import tbt_bpm_class

home_dir = os.path.join(
    "/Volumes/Ext. HD (Dropbox)", "Dropbox", "Johan", "MAX V", "MAX IV",
    "T-b-T_BPM_data")

# Allocate the tbt_class.
tbt = tbt_bpm_class()

file_name = home_dir+"/TbTData_20230417.h5"

tbt.rd_tbt_MAX_IV(file_name)

if False:
    bpm_no = 0
    cut    = 322
    n_data = 500

    # Convert from [nm] to [mm].
    tbt._tbt_data = 1e-6*tbt._tbt_data_buf[:, cut:cut+n_data, bpm_no]

    tbt.analyse_tbt_bpm_data(1, True, True, True)

    assert False

if not False:
    cut    = 322
    n_data = 500

    # Convert from [nm] to [mm].
    tbt._tbt_data_buf = 1e-6*tbt._tbt_data_buf[:, cut:cut+n_data, :]

    # If plot is set to True (last parameter):
    #   plot for first BPM and then terminate. 
    tbt.compute_lin_opt_MAX_IV(cut, True, [True, True], False, False)
