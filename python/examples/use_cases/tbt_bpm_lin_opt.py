'''
Author:

  Johan Bengtsson
  16/02/24

Tbt_bpm class to analyse turn-by-turn BPM data to extract the linear optics -
beta functions & phase advance - at the BPMs.
'''


import os
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from thor_scsi.utils.fft import fft_class


X_, Y_, Z_ = [0, 1, 2]


#-------------------------------------------------------------------------------

class tbt_bpm:
    # Private.

    def __init__(self):
        self._tbt_data     = 0
        self._tbt_data_fft = 0
        self._tbt_data     = 0
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

    def plt_tbt_fft(self, nu, A):
        n = len(A[X_])

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
            (nu[0:n//2+1], A[X_][0:n//2+1], linefmt="b-", markerfmt="",
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
            (nu[0:n//2+1], A[Y_][0:n//2+1], linefmt="r-", markerfmt="",
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

    def analyse_tbt_bpm_data(self, rm_avg, prt, plot):
        fft = fft_class()

        n_peak = 3
        n_max  = 5
        n_data = len(self._tbt_data[0])

        if rm_avg:
            for k in range(2):
                # Remove average.
                self._tbt_data[k] -= np.mean(self._tbt_data[k])

        self._tbt_data_fft = np.zeros([2, n_data], dtype="complex")
        A_fft = np.zeros([2, n_data], dtype="float")
        f = np.array(range(n_data))/n_data
        sine_window = sp.signal.windows.cosine(n_data)

        nu = np.zeros((2, n_peak), dtype=float)
        A = np.zeros((2, n_peak), dtype=float)
        j = np.zeros(n_peak, dtype=int)
        phi = np.zeros((2, n_peak), dtype=float)
        for i in range(2):
            # Use [mm].
            self._tbt_data_fft[i] = \
                sp.fft.fft(self._tbt_data[i]*sine_window)/n_data
            A_fft[i] = abs(self._tbt_data_fft[i])
            nu[i], A[i], j = fft.get_peak_sin(A_fft[i], n_peak)
            for k in range(n_peak):
                phi[i][k] = fft.get_phase(j[k], nu[i][k], self._tbt_data[i])

        if prt:
            print("\nnu = [{:7.5f}/{:7.5f}, {:7.5f}/{:7.5f}]".
                  format(nu[X_][0], 1e0-nu[X_][0], nu[Y_][0], 1e0-nu[Y_][0]))
            print("A  = [{:5.3f}, {:5.3f}]".format(A[X_][0], A[Y_][0]))
            print("\nHorizontal Plane")
            print("    nu         A        phi   n_x  n_y    eps")
            for k in range(1, n_peak):
                n_x, n_y, eps = \
                    fft.find_harmonic(n_max, nu[X_][k], nu[Y_][k], nu[X_][k])
                print("  {:7.5f}  {:9.3e}  {:6.1f}   {:1d}    {:1d}   {:7.1e}".
                      format(nu[X_][k], A[X_][k], np.rad2deg(phi[X_][k]), n_x,
                             n_y, eps))
            print("\nVertical Plane")
            print("    nu         A        phi   n_x  n_y    eps")
            for k in range(1, n_peak):
                n_x, n_y, eps = \
                    fft.find_harmonic(n_max, nu[X_][k], nu[Y_][k], nu[Y_][k])
                print("  {:7.5f}  {:9.3e}  {:6.1f}   {:1d}    {:1d}   {:7.1e}".
                      format(nu[Y_][k], A[Y_][k], np.rad2deg(phi[Y_][k]), n_x,
                             n_y, eps))

        if plot:
            self.plt_tbt_fft(f, A_fft)
            plt.show()

        return nu, A, phi

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
            nu, A, phi = tbt.analyse_tbt_bpm_data(rm_avg, prt, plot)
            for k in range(2):
                # Aliasing if the tune is > 0.5.
                if alias[k]:
                    phi[k] = -phi[k]
            # Collect the results.
            self._bpm.update({name : {"nu" : nu, "A" : A, "phi" : phi}})
            if plot:
                break

        self.prt_lin_opt()

#-------------------------------------------------------------------------------

# Main program.

home_dir = os.path.join(os.environ["HOME"])

# Allocate the tbt_class.
tbt = tbt_bpm()

if not False:
    # Test case - analyse one BPM.
    file_name = [
        home_dir+"/Teresia/20240122_TbT_hor.dat",
        home_dir+"/Teresia/20240122_TbT_ver.dat"
    ]

    # Skip initial transient.
    n_data, cut = [2**10, 9999]

    tbt.rd_tbt_txt(file_name, n_data, cut)
    tbt.analyse_tbt_bpm_data(True, True, True)
    
if False:
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
        tbt.compute_lin_opt(cut, True, [True, True], not False, False)
