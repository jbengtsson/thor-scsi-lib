'''
Author:

  Johan Bengtsson
  26/02/24

Script to measure the linear chromaticity by extracting:

  nu( f_RF )

from turn-by-turn BPM data.
'''


import numpy as np

from thor_scsi.utils.bessy_ii_mml import middle_layer
from thor_scsi.utils.tbt_bpm import tbt_bpm_class


X_, Y_, Z_ = [0, 1, 2]

#-------------------------------------------------------------------------------

class epics_class:
    # Private.

    def __init__(self):
        self._on_line = False
        self._rd_only = True

    # Public.

    def init(self, on_line, rd_only):
        self._on_line = on_line
        self._rd_only = rd_only
        if on_line:
            from epics import caget, caput, cainfo
        
    def get_pv_sp(self, pv):
        prt = False
        if self._on_line: 
            val = caget(pv)
        else:
            val = 0.1234
        if prt:
            print("  {:10s} {:8.5f}".format(pv, val))
        return val

    def put_pv_sp(self, pv, val):
        prt = False
        if self._on_line:
            if not rd_only:
                caput(pv, val)
            else:
                print("\nput_pv_sp – rd_only mode")
        elif prt:
            print("  {:10s} {:8.5f}".format(pv, val))
            
    def get_pv_rb(self, pv):
        prt = False
        if self._on_line:
            val = caget(pv)
        else:
            val = 0.1234
        if prt:
            print("  {:10s} {:8.5f}".format(pv, val))
        return val            
            

    def meas_chrom(self, n_point, f_RF_step, pv_sp_RF, alpha_c):
        tbt = tbt_bpm_class()

        if not self._on_line:
            home_dir = "/Users/johan/Teresia"
        else:
            home_dir = \
                "/net/nfs/srv/MachinePhysics/MachineDevelopment/mcateer/" \
                "Jupyter/olsson-test"

        file_name = [
            home_dir+"/20240122_TbT_hor.dat",
            home_dir+"/20240122_TbT_ver.dat"
        ]

        # Skip initial transient.
        cut, n_data = [9999, 2**10]

        if self._on_line:
            f_0_RF = self.get_pv_sp(pv_sp_RF)
        else:
            f_0_RF = 500e6
        print("\nmeas_chrom:\n  f_0_RF [MHz] = {:10.6f}".format(1e-6*f_0_RF))

        nu = np.zeros((2, 2*n_point+1))
        delta = np.zeros(2*n_point+1)
        df_RF = -n_point*f_RF_step
        print("\n        f_RF      delta   nu_x    nu_y")
        print("        [MHz]      [%]")
        for k in range(0, 2*n_point+1):
            df_RF += f_RF_step
            f_RF = f_0_RF + df_RF
            self.put_pv_sp(pv_sp_RF, f_RF)
            delta[k] = df_RF/(f_0_RF*alpha_c)

            tbt.rd_tbt_txt(file_name, n_data, cut)
            tbt.analyse_tbt_bpm_data(1, True, False, False)
            nu[X_][k] = tbt._f[X_, 0] + 1e0*delta[k]
            nu[Y_][k] = tbt._f[Y_, 0] - 1e0*delta[k]

            print("  {:2d} {:10.6f} {:8.5f} {:7.5f} {:7.5f}".
                  format(k-n_point, 1e-6*f_RF, 1e2*delta[k], nu[X_][k],
                         nu[Y_][k]))

        n_pol = 2
        pol = np.zeros((2, n_pol+1))
        for k in range(2):
            pol[k] = np.polyfit(delta, nu[k, :], n_pol)
        print("\n  ksi = [{:9.3e}, {:9.3e}]".format(pol[X_][1], pol[Y_][1]))

#-------------------------------------------------------------------------------

# Main program.

on_line = False
rd_only = not True

mml = epics_class()

mml.init(on_line, rd_only)

if False:
    cav_RF = "cav_RF"
    init_sp = mml.get_pv_sp(cav_RF)
    scaling = 1.01
    mml.put_pv_sp(cav_RF, init_sp*scaling)
    mml.get_pv_rb(cav_RF)
    mml.put_pv_sp(cav_RF, init_sp)

n_point   = 5
f_RF_step = 1e2
pv_sp_RF  = "cav_RF"
alpha_c   = 7.038e-04

mml.meas_chrom(n_point, f_RF_step, pv_sp_RF, alpha_c)
