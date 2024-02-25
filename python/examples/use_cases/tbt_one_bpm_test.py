
import os
import numpy as np
import scipy as sp

from thor_scsi.utils.tbt_bpm import tbt_bpm_class


home_dir = os.path.join(os.environ["HOME"])
file_name = [
    home_dir+"/Teresia/20240122_TbT_hor.dat",
    home_dir+"/Teresia/20240122_TbT_ver.dat"
]

tbt = tbt_bpm_class()

# Skip initial transient.
cut, n_data = [9999, 2**10]

tbt.rd_tbt_txt(file_name, n_data, cut)

tbt.analyse_tbt_bpm_data(1, cut, n_data, True)
