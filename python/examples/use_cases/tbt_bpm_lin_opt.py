
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from thor_scsi.utils.fft import fft_class

#-------------------------------------------------------------------------------

class tbt_bpm:
    def __init__(self):
        self.n_data   = 0
        self.cut      = 0
        self.tbt_data = np.zeros((2, 1))

    def rd_tbt(self, file_name, k):
        # First number is the number of data.
        buf = []
        with open(file_name[k]) as f:
            for line in f:
                buf = np.array(line.split())
        # Apply cut & convert from [nm] to [mm].
        self.tbt_data[k] = \
            1e-6*(np.array(buf, dtype='float')[self.cut:self.cut+self.n_data])

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
        print(f"File saved as: {file_name:s}")

    def analyse_tbt_bpm_data(self, file_name, cut, n_data, plot):
        fft = fft_class()

        self.n_data = n_data
        self.cut = cut
        self.tbt_data = np.zeros([2, n_data], dtype="float")

        for k in range(2):
            self.rd_tbt(file_name, k)
            self.tbt_data[k] -= np.mean(self.tbt_data[k])

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

home_dir = os.path.join(os.environ["HOME"])
file_name = [
    home_dir+"/Teresia/20240122_TbT_hor.dat",
    home_dir+"/Teresia/20240122_TbT_ver.dat"
]

tbt = tbt_bpm()

# Skip initial transient.
cut, n_data = [9999, 2**10]
tbt.analyse_tbt_bpm_data(file_name, cut, n_data, True)
