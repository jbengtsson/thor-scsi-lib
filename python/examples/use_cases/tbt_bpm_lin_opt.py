'''
Author:

  Johan Bengtsson
  17/02/24

Tbt_bpm class to analyse turn-by-turn BPM data to extract the linear optics -
beta functions & phase advance at the BPMs.
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
    def __init__(self):
        self.tbt_data     = 0
        self.tbt_data_fft = 0
        self.tbt_data     = 0
        self.bpm          = {}

    def plt_data(self, file_name):
        fig, (gr_1, gr_2) = plt.subplots(2)
        fig.set_size_inches(16, 10)

        for bpm, x, y in zip(
                self.tbt_data["BPM_name"], self.tbt_data["X"],
                self.tbt_data["Y"]):
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
        print("Plot saved as: {:s}".format(file_name))

    def plt_tbt_fft(self, nu, A, plane):
        n = len(self.tbt_data[plane])

        fig, (gr_1, gr_2) = plt.subplots(2)

        if plane == 0:
            file_name = "tbt_hor.png"

            gr_1.set_title("Horizontal Position")
            gr_2.set_title("FFT of Horizontal Position")
            gr_1.set_xlabel("Turn Number")
            gr_1.set_ylabel(r"x [mm]")
            gr_1.stem(
                self.tbt_data[plane], linefmt="b-", markerfmt="", basefmt="k-")

            gr_2.set_xlabel("f")
            gr_2.set_ylabel(r"$A_x$")
            gr_2.stem \
                (nu[0:n//2+1], A[0:n//2+1], linefmt="b-", markerfmt="",
                 basefmt="k-")
        else:
            file_name = "tbt_ver.png"

            gr_1.set_title("Vertical Position")
            gr_2.set_title("FFT of Vertical Position")
            gr_1.set_xlabel("Turn Number")
            gr_1.set_ylabel(r"y [mm]")
            gr_1.stem(
                self.tbt_data[plane], linefmt="r-", markerfmt="", basefmt="k-")

            gr_2.set_xlabel("f")
            gr_2.set_ylabel(r"$A_y$")
            gr_2.stem \
                (nu[0:n//2+1], A[0:n//2+1], linefmt="r-", markerfmt="",
                 basefmt="k-")

        fig.tight_layout()

        plt.savefig(file_name)
        print("Plot saved as: {:s}".format(file_name))

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


    def analyse_tbt_bpm_data(self, rm_avg, prt, plot):
        fft = fft_class()

        n_data = len(self.tbt_data[0])

        if rm_avg:
            for k in range(2):
                # Remove average.
                self.tbt_data[k] -= np.mean(self.tbt_data[k])

        self.tbt_data_fft = np.zeros([2, n_data], dtype="complex")
        A_fft = np.zeros([2, n_data], dtype="float")
        f = np.array(range(n_data))/n_data
        sine_window = sp.signal.windows.cosine(n_data)

        nu = np.zeros(2, dtype=float)
        A = np.zeros(2, dtype=float)
        phi = np.zeros(2, dtype=float)
        for j in range(2):
            # Use [mm].
            self.tbt_data_fft[j] = \
                sp.fft.fft(self.tbt_data[j]*sine_window)/n_data
            A_fft[j] = abs(self.tbt_data_fft[j])
            nu_buf, A_buf, k = fft.get_peak_sin(A_fft[j], 1)
            phi[j] = fft.get_phase(k, nu_buf[0], self.tbt_data[j])
            nu[j], A[j] = nu_buf[0], A_buf[0]
            if prt:
                print("nu = [{:8.6f}, {:8.6f}] A = {:9.3e} phi = {:5.1f}".
                      format(nu_buf[0], 1e0-nu_buf[0], A_buf[0],
                             np.rad2deg(phi[j])))

        if plot:
            print()
            for j in range(2):
                self.plt_tbt_fft(f, A_fft[j], j)
            plt.show()

        return nu, A, phi

    def prt_lin_opt(self):
        # Print the results from compute_lin_opt.
        first = True
        A_k = np.zeros(2, dtype=float)
        A_rel = np.zeros(2, dtype=float)
        phi_k = np.zeros(2, dtype=float)
        dphi = np.zeros(2, dtype=float)
        print("\n  BPM     Rel. A_x  Rel. A_y  dnu_x  dnu_y");
        for name in self.bpm:
            if first:
                A_k = self.bpm[name]["A"]
                phi_k = self.bpm[name]["phi"]
                first = False
            A_rel = self.bpm[name]["A"]/A_k
            dphi = self.bpm[name]["phi"] - phi_k
            for k in range(2):
                if dphi[k] < 0e0:
                    dphi[k] += 2e0*np.pi
            A_k = self.bpm[name]["A"]
            phi_k = self.bpm[name]["phi"]
            print("{:10s} {:5.3f}     {:5.3f}    {:5.3f}  {:5.3f}".
                  format(name, A_rel[X_], A_rel[Y_], dphi[X_]/(2e0*np.pi),
                         dphi[Y_]/(2e0*np.pi)))


    def compute_lin_opt(self, cut, rm_avg, alias, prt, plot):
        self.bpm = {}

        # Loop over BPMs.
        for name, x, y in zip(
                self.tbt_data["BPM_name"], self.tbt_data["X"],
                self.tbt_data["Y"]):
            if prt:
                print("\n{:10s}".format(name))
            n_data = len(x)
            self.tbt_data = \
                np.array([1e-6*x[cut:cut+n_data], 1e-6*y[cut:cut+n_data]],
                         dtype=float)
            nu, A, phi = tbt.analyse_tbt_bpm_data(rm_avg, prt, plot)
            for k in range(2):
                # Aliasing if the tune is > 0.5.
                if alias[k]:
                    phi[k] = -phi[k]
            # Collect the results.
            self.bpm.update({name : {"nu" : nu, "A" : A, "phi" : phi}})
            if plot:
                break

        self.prt_lin_opt()

#-------------------------------------------------------------------------------

# Main program.

home_dir = os.path.join(os.environ["HOME"])

# Allocate the tbt_class.
tbt = tbt_bpm()

if False:
    # Analyse a single BPM.
    file_name = [
        home_dir+"/Teresia/20240122_TbT_hor.dat",
        home_dir+"/Teresia/20240122_TbT_ver.dat"
    ]

    # Skip initial transient.
    n_data, cut = [2**10, 9999]

    tbt.rd_tbt_txt(file_name, n_data, cut)
    tbt.analyse_tbt_bpm_data(True, True, True)
    
if not False:
    # Compute the linear optics.
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
