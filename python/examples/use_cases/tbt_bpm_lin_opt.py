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

#-------------------------------------------------------------------------------

class tbt_bpm:
    def __init__(self):
        self.tbt_data = 0

    def plt_tbt(self, nu, A, plane):
        # Turn interactive mode off.
        # plt.ioff()

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
        print(f"Plot saved as: {file_name:s}")

    def analyse_tbt_bpm_data(self, plot):
        fft = fft_class()

        n_data = len(self.tbt_data[0])

        tbt_data_fft = np.zeros([2, n_data], dtype="complex")
        A_fft = np.zeros([2, n_data], dtype="float")
        f = np.array(range(n_data))/n_data
        sine_window = sp.signal.windows.cosine(n_data)

        print()
        for j in range(2):
            # Use [mm].
            tbt_data_fft[j] = sp.fft.fft(self.tbt_data[j]*sine_window)/n_data
            A_fft[j] = abs(tbt_data_fft[j])
            nu, A, k = fft.get_peak_sin(A_fft[j], 1)
            phi = fft.get_phase(k, nu[0], self.tbt_data[j])
            print("nu = [{:8.6f}, {:8.6f}] A = {:9.3e} phi = {:5.1f}".
                  format(nu[0], 1e0-nu[0], A[0], phi*180e0/sp.pi))

        if plot:
            print()
            for j in range(2):
                self.plt_tbt(f, A_fft[j], j)
            plt.show()

#-------------------------------------------------------------------------------

def plt_data(data, file_name):
    fig, (gr_1, gr_2) = plt.subplots(2)
    fig.set_size_inches(16, 10)

    for bpm, x in zip(data["BPM_name"], data["X"]):
        gr_1.plot(x*1e-6, label=bpm)

    gr_1.set_xlim(50, 150)
    gr_1.grid()
    gr_1.legend(
        fontsize="xx-large", bbox_to_anchor=(1.05, 1.0), loc='upper left')
    gr_1.set_title("Horizontal Displacement", fontsize="xx-large")
    gr_1.set_xlabel("Turn Number", fontsize="xx-large")
    gr_1.set_ylabel("x [mm]", fontsize="xx-large")
    gr_1.tick_params(axis='x', labelsize="xx-large")
    gr_1.tick_params(axis='y', labelsize="xx-large")

    for bpm, y in zip(data["BPM_name"], data["Y"]):
        gr_2.plot(y*1e-6, label=bpm)

    gr_2.set_xlim(50, 150)
    gr_2.grid()
    gr_2.legend(
        fontsize="xx-large", bbox_to_anchor=(1.05, 1.0), loc='upper left')
    gr_2.set_title("Vertical Displacement", fontsize="xx-large")
    gr_2.set_xlabel("Turn Number", fontsize="xx-large")
    gr_2.set_ylabel("y [mm]", fontsize="xx-large")

    fig.tight_layout()

    plt.savefig(file_name)
    print(f"Plot saved as: {file_name:s}")


def rd_tbt_txt(file_name, n_data, cut):
    # First number is the number of data.
    buf = []
    with open(file_name) as f:
        for line in f:
            buf = np.array(line.split())
    # Apply cut & convert from [nm] to [mm].
    tbt_data = 1e-6*(np.array(buf, dtype='float')[cut:cut+n_data])
    return tbt_data

def rd_tbt_df(file_name):
    data = pd.read_hdf(file_name, key="df")
    return data


tbt = tbt_bpm()

if not False:
    home_dir = os.path.join(os.environ["HOME"])
    file_name = [
        home_dir+"/Teresia/20240122_TbT_hor.dat",
        home_dir+"/Teresia/20240122_TbT_ver.dat"
    ]

    # Skip initial transient.
    n_data, cut = [2**10, 9999]

    tbt.tbt_data = np.zeros((2, n_data))

    for k in range(2):
        tbt.tbt_data[k] = rd_tbt_txt(file_name[k], n_data, cut)
        tbt.tbt_data[k] -= np.mean(tbt.tbt_data[k])
    tbt.analyse_tbt_bpm_data(True)
    
if False:
    file_name = \
        home_dir \
        +"/Teresia/Markus/20240216_testset_injectionkickers_storedbeam.hdf"

    data = rd_tbt_df(file_name)
    for bpm in data["BPM_name"]:
        print(bpm)

        plt_data(data, "20240216_firstTbTdata.png")
        plt.show()
