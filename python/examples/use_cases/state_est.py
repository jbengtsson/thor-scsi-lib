'''
Author:

  Johan Bengtsson
  17/02/24

State_space_est class to reconstruct the state-space/phase-space from
turn-by-turn BPM data.
'''


import os
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from thor_scsi.utils.fft import fft_class


X_, Y_, Z_ = [0, 1, 2]


#-------------------------------------------------------------------------------

class state_space:
    def __init__(self):
        self.bpm  = ["BPMZ7D8R", "BPMZ1T8R"]
        self.beta = np.array([[5.986, 14.713], [14.626, 7.136]])
        self.dnu  = np.array([16.10234-16.08547, 6.06078-6.03863])
        self.sin  = sin(2e0*np.pi*self.dnu)
        self.cos  = cos(2e0*np.pi*self.dnu)

    def rd_tbt_txt(self, file_name, n_data, cut):
        self.tbt_data = np.zeros((2, n_data), dtype=float)
        for k in range(2):
            # First value is the number of data.
            with open(file_name[k]) as f:
                for line in f:
                    buf = np.array(line.split())
                    # Apply cut & convert from [nm] to [mm].
                    self.tbt_data[k] =  \
                        1e-6*(np.array(buf, dtype='float')[cut:cut+n_data])
 
    def rd_tbt_df(self, file_name):
        self.tbt_data = pd.read_hdf(file_name, key="df")

    def state_space_est(self, xy):
        print("\n A_rel =  [{:7.5f}, {:7.5f}] dnu = [{:7.5f}, {:7.5f}]".
              format(np.sqrt(self.beta[1][X_]/self.beta[0][X_]),
                     np.sqrt(self.beta[1][Y_]/self.beta[0][Y_]),
                     self.dnu[X_], self.dnu[Y_]))
        ps_Fl = np.zeros(4)
        twoJ = np.zeros(2)
        # for k in range(2):
        #     ps_Fl[2*k] = xy[0][k]/self.beta[0][k]
        #     ps_Fl[2*k+1] = \
        #         (np.sqrt(self.beta[0][k]/self.beta[1][k])*xy[1][k] \
        #          -xy[0][k]*self.cos[k]) \
        #         /(np.sqrt(self.beta[0][k])*self.sin[k])
        #     twoJ[k] = ps_Fl[2*k]**2 + ps_Fl[2*k+1]**2
        return twoJ


#-------------------------------------------------------------------------------

# Main program.

home_dir = os.path.join(os.environ["HOME"])

# Allocate the tbt_class.
ss = state_space_est()

file_name = \
    home_dir \
    +"/Teresia/Markus/20240216_testset_injectionkickers_storedbeam.hdf"

ss.rd_tbt_df(file_name)

cut = 108

for name, x, y in zip(
        ss.tbt_data["BPM_name"], self.tbt_data["X"], self.tbt_data["Y"]):
    pass

print(ss.tbt_data[tbt.bpm[0]])
