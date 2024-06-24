'''
Author:

  Johan Bengtsson
  16/02/24

Tbt_bpm class to analyse turn-by-turn BPM data to extract the linear optics -
beta functions & phase advance - at the BPMs.
'''


import os
import math
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from thor_scsi.utils.tbt_bpm import tbt_bpm_class


home_dir = os.path.join(os.environ["HOME"])

# Allocate the tbt_class.
tbt = tbt_bpm_class()

if False:
    # Test case - analyse one BPM.
    file_name = [
        home_dir+"/Teresia/20240122_TbT_hor.dat",
        home_dir+"/Teresia/20240122_TbT_ver.dat"
    ]

    # Skip initial transient.
    n_data, cut = [2**10, 9999]

    tbt.rd_tbt_txt(file_name, n_data, cut)
    tbt.analyse_tbt_bpm_data(3, True, True, True)

    assert False

if not False:
    # Extract the linear optics.
    file_name = \
        home_dir \
        +"/Teresia/Markus/20240216_testset_injectionkickers_storedbeam.hdf"

    tbt.rd_tbt_df(file_name)

    if False:
        tbt.plt_data("20240216_firstTbTdata.png")
        plt.show()

    if not False:
        cut = 108
        # If plot is set to True (last parameter):
        #   plot for first BPM and then terminate. 
        tbt.compute_lin_opt(cut, True, [True, True], False, False)
