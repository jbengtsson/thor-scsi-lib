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

from thor_scsi.utils.fft import fft_class


X_, Y_, Z_ = [0, 1, 2]


#-------------------------------------------------------------------------------

class tbt_bpm_class:
    # Private.

    def __init__(self):
        self._tbt_data     = math.nan
        self._tbt_data_fft = math.nan
        self._A_fft        = math.nan
        self._nu           = np.array((math.nan, math.nan))
        self._A            = math.nan
        self._f            = math.nan
        self._phi          = math.nan
        self._bpm          = {}

    # Public.

    def plt_data(self, file_name):
        fig, (gr_1, gr_2) = plt.subplots(2)
        fig.set_size_inches(16, 10)

        for bpm, x, y in zip(
                self._tbt_data["BPM_name"], self._tbt_data["X"],
                self._tbt_data["Y"]):
            gr_1.plot(x*1e-6, label=bpm)
            gr_2.plot(y*1e-6, label=bpm)

        gr_1.set_xlim(50, 150)
        gr_1.grid()
        gr_1.legend(
            fontsize="xx-large", bbox_to_anchor=(1.05, 1.0), loc='upper left')
        gr_1.set_title("Horizontal Displacement", fontsize="xx-large")
        gr_1.set_xlabel("Turn Number", fontsize="xx-large")
        gr_1.set_ylabel("x [mm]", fontsize="xx-large")
        gr_1.tick_params(axis='x', labelsize="xx-large")
        gr_1.tick_params(axis='y', labelsize="xx-large")

        gr_2.set_xlim(50, 150)
        gr_2.grid()
        gr_2.legend(
            fontsize="xx-large", bbox_to_anchor=(1.05, 1.0), loc='upper left')
        gr_2.set_title("Vertical Displacement", fontsize="xx-large")
        gr_2.set_xlabel("Turn Number", fontsize="xx-large")
        gr_2.set_ylabel("y [mm]", fontsize="xx-large")

        fig.tight_layout()

        plt.savefig(file_name)
        print("\nPlot saved as: {:s}".format(file_name))

    def plt_tbt_fft(self, A):
        n = len(A[X_])

        f = np.array(range(n))/n

        fig, ((gr_1, gr_2), (gr_3, gr_4)) = plt.subplots(2, 2)

        file_name = "tbt_fft.png"

        fig.suptitle("FFT of Turn-by-Turn BPM Data")

        gr_1.set_title("Horizontal Position")
        gr_1.set_xlabel("Turn Number")
        gr_1.set_ylabel(r"x [mm]")
        gr_1.stem(
            self._tbt_data[X_], linefmt="b-", markerfmt="", basefmt="k-")

        gr_3.set_title("FFT of Horizontal Position")
        gr_3.set_xlabel("f")
        gr_3.set_ylabel(r"$A_x$")
        gr_3.stem \
            (f[0:n//2+1], A[X_][0:n//2+1], linefmt="b-", markerfmt="",
             basefmt="k-")

        gr_2.set_title("Vertical Position")
        gr_2.set_xlabel("Turn Number")
        gr_2.set_ylabel(r"y [mm]")
        gr_2.stem(
            self._tbt_data[Y_], linefmt="r-", markerfmt="", basefmt="k-")

        gr_4.set_title("FFT of Vertical Position")
        gr_4.set_xlabel("f")
        gr_4.set_ylabel(r"$A_y$")
        gr_4.stem \
            (f[0:n//2+1], A[Y_][0:n//2+1], linefmt="r-", markerfmt="",
             basefmt="k-")

        fig.tight_layout()

        plt.savefig(file_name)
        print("\nPlot saved as: {:s}".format(file_name))

    def rd_tbt_txt(self, file_name, n_data, cut):
        self._tbt_data = np.zeros((2, n_data), dtype=float)
        for k in range(2):
            # First value is the number of data.
            with open(file_name[k]) as f:
                for line in f:
                    buf = np.array(line.split())
                    # Apply cut & convert from [nm] to [mm].
                    self._tbt_data[k] =  \
                        1e-6*(np.array(buf, dtype='float')[cut:cut+n_data])
 
    def rd_tbt_df(self, file_name):
        self._tbt_data = pd.read_hdf(file_name, key="df")

    def prt_f(self, fft, n_max, n_peak):
        print("\nHorizontal Plane")
        print("     f       1-f        A        phi   n_x  n_y    eps")
        for k in range(0, n_peak):
            n_x, n_y, eps = \
                fft.find_harmonic(
                    n_max, self._nu[X_], self._nu[Y_], self._f[X_][k])
            print("  {:7.5f}  {:7.5f}  {:9.3e}  {:6.1f}   {:1d}   {:2d}"
                  "   {:7.1e}".
                  format(self._f[X_][k], 1e0-self._f[X_][k], self._A[X_][k],
                         np.rad2deg(self._phi[X_][k]), n_x, n_y, eps))
        print("\nVertical Plane")
        print("     f       1-f        A        phi   n_x  n_y    eps")
        for k in range(0, n_peak):
            n_x, n_y, eps = \
                fft.find_harmonic(
                    n_max, self._nu[X_], self._nu[Y_], self._f[Y_][k])
            print("  {:7.5f}  {:7.5f}  {:9.3e}  {:6.1f}   {:1d}   {:2d}"
                  "   {:7.1e}".
                  format(self._f[Y_][k], 1e0-self._f[Y_][k], self._A[Y_][k],
                         np.rad2deg(self._phi[Y_][k]), n_x, n_y, eps))

    def analyse_tbt_bpm_data(self, n_peak, rm_avg, prt, plot):
        fft = fft_class()

        n_max  = 4
        n_data = len(self._tbt_data[0])

        if rm_avg:
            for k in range(2):
                # Remove average.
                self._tbt_data[k] -= np.mean(self._tbt_data[k])

        self._tbt_data_fft = np.zeros([2, n_data], dtype="complex")
        self._A_fft = np.zeros([2, n_data], dtype="float")
        sine_window = sp.signal.windows.cosine(n_data)

        self._f = np.zeros((2, n_peak), dtype=float)
        self._A = np.zeros((2, n_peak), dtype=float)
        self._phi = np.zeros((2, n_peak), dtype=float)
        for k in range(2):
            # Use [mm].
            self._tbt_data_fft[k] = \
                sp.fft.fft(self._tbt_data[k]*sine_window)/n_data
            self._A_fft[k] = abs(self._tbt_data_fft[k])
            self._f[k], self._A[k], ind_2 = \
                fft.get_peak_sin(self._A_fft[k], n_peak)
            self._phi[k] = fft.get_phase(ind_2, self._f[k], self._tbt_data[k])

        if prt:
            self._nu = np.array([self._f[X_][0], self._f[Y_][0]])
            self.prt_f(fft, n_max, n_peak)

        if plot:
            self.plt_tbt_fft(self._A_fft)
            plt.show()

    def prt_lin_opt(self):
        # Print the results from compute_lin_opt.
        first = True
        A_k = np.zeros(2, dtype=float)
        beta_rel = np.zeros(2, dtype=float)
        phi_k = np.zeros(2, dtype=float)
        dphi = np.zeros(2, dtype=float)
        print("\n  BPM     Rel. beta_x  Rel. beta_y  dnu_x  dnu_y");
        for name in self._bpm:
            if first:
                A_k = self._bpm[name]["A"]
                phi_k = self._bpm[name]["phi"]
                first = False
            beta_rel = (self._bpm[name]["A"]/A_k)**2
            dphi = self._bpm[name]["phi"] - phi_k
            for k in range(2):
                if dphi[k] < 0e0:
                    dphi[k] += 2e0*np.pi
            A_k = self._bpm[name]["A"]
            phi_k = self._bpm[name]["phi"]
            print("{:10s}   {:5.3f}        {:5.3f}     {:5.3f}  {:5.3f}".
                  format(name, beta_rel[X_], beta_rel[Y_],
                         dphi[X_]/(2e0*np.pi), dphi[Y_]/(2e0*np.pi)))

    def compute_lin_opt(self, cut, rm_avg, alias, prt, plot):
        self._bpm = {}

        # Loop over BPMs.
        for name, x, y in zip(
                self._tbt_data["BPM_name"], self._tbt_data["X"],
                self._tbt_data["Y"]):
            if prt:
                print("\n{:10s}".format(name))
            n_data = len(x)
            # Convert from [nm] to [mm].
            self._tbt_data = \
                np.array([1e-6*x[cut:cut+n_data], 1e-6*y[cut:cut+n_data]],
                         dtype=float)
            self.analyse_tbt_bpm_data(1, rm_avg, prt, plot)
            for k in range(2):
                # Aliasing if the tune is > 0.5.
                if alias[k]:
                    self._phi[k, 0] = -self._phi[k, 0]
            # Collect the results.
            self._bpm.update({
                name : {
                    "nu" :self._f[:, 0], "A" : self._A[:, 0],
                    "phi" : self._phi[:, 0]
                }
            })
            if plot:
                break

        self.prt_lin_opt()

    def compute_lin_opt_MAX_IV(self, cut, rm_avg, alias, prt, plot):
        self._bpm = {}

        # Loop over BPMs.
        for bpm in range(len(self._tbt_data_buf[0, 0, :])):
            self._tbt_data = self._tbt_data_buf[:, :, bpm]
            self.analyse_tbt_bpm_data(1, rm_avg, prt, plot)
            for k in range(2):
                # Aliasing if the tune is > 0.5.
                if alias[k]:
                    self._phi[k, 0] = -self._phi[k, 0]
            # Collect the results.
            self._bpm.update({
                str(bpm) : {
                    "nu" :self._f[:, 0], "A" : self._A[:, 0],
                    "phi" : self._phi[:, 0]
                }
            })
            if plot:
                break

        self.prt_lin_opt()

       
__all__ = ["tbt_bpm_class"]
