'''
Author:

  Johan Bengtsson
  16/02/24

State_est class to reconstruct the state-space/phase-space from turn-by-turn
BPM data.
'''


import os
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from thor_scsi.utils.fft import fft_class


X_, Y_, Z_                      = [0, 1, 2]
[x_, px_, y_, py_, ct_, delta_] = [0, 1, 2, 3, 4, 5]


#-------------------------------------------------------------------------------

class state_space:
    # Private.

    def __init__(self):
        self._bpm      = ["BPMZ6D8R", "BPMZ1T8R"]
        self._beta     = np.array([[4.786, 16.506], [5.98559, 14.713]])
        self._beta_rel = \
            np.array([
                self._beta[0][X_]/self._beta[1][X_],
                self._beta[0][Y_]/self._beta[1][Y_]
            ])
        self._dnu      = np.array([16.08547-15.68144, 6.03863-6.01883])
        self._sin      = np.sin(2e0*np.pi*self._dnu)
        self._cos      = np.cos(2e0*np.pi*self._dnu)

        self._tbt_data = 0
        self._xy_0     = 0
        self._xy_1     = 0
        self._ps_Fl    = 0
        self._twoJ     = 0
        self._twoJ_fft = 0

    # Public.

    def plt_Floquet_space(self):
        fig, (gr_1, gr_2) = plt.subplots(1, 2)

        file_name = "Floquet_space.png"

        fig.suptitle("State Estimator -  Floquet Space")

        gr_1.set_title("Horizontal")
        gr_1.set_xlabel(r"$\tilde{x}$")
        gr_1.set_ylabel(r"$\tilde{p}_x$")
        gr_1.plot(self._ps_Fl[x_], self._ps_Fl[px_], "b+")

        gr_2.set_title("Vertical")
        gr_2.set_xlabel(r"$\tilde{y}$")
        gr_2.set_ylabel(r"$\tilde{p}_y$")
        gr_2.plot(self._ps_Fl[y_], self._ps_Fl[py_], "r+")

        fig.tight_layout()

        plt.savefig(file_name)
        print("\nPlot saved as: {:s}".format(file_name))

    def plt_tbt_fft(self, A):
        n = len(self._twoJ[X_])

        f = np.array(range(n))/n

        fig, ((gr_1, gr_2), (gr_3, gr_4)) = plt.subplots(2, 2)

        file_name = "tbt_fft.png"

        fig.suptitle("FFT of 2J")

        gr_1.set_title(r"$2J_x$")
        gr_1.set_xlabel("Turn Number")
        gr_1.set_ylabel("")
        gr_1.stem(
            self._twoJ[X_], linefmt="b-", markerfmt="", basefmt="k-")

        gr_3.set_title(r"FFT of $2J_x$")
        gr_3.set_xlabel("f")
        gr_3.set_ylabel(r"$A_x$")
        gr_3.stem \
            (f[0:n//2+1], A[X_][0:n//2+1], linefmt="b-", markerfmt="",
             basefmt="k-")

        gr_2.set_title(r"$2J_y$")
        gr_2.set_xlabel("Turn Number")
        gr_2.set_ylabel(r"")
        gr_2.stem(
            self._twoJ[Y_], linefmt="r-", markerfmt="", basefmt="k-")

        gr_4.set_title(r"FFT of $2J_y$")
        gr_4.set_xlabel("f")
        gr_4.set_ylabel(r"$A_y$")
        gr_4.stem \
            (f[0:n//2+1], A[Y_][0:n//2+1], linefmt="r-", markerfmt="",
             basefmt="k-")

        fig.tight_layout()

        plt.savefig(file_name)
        print("\nPlot saved as: {:s}".format(file_name))

    def rd_tbt_df(self, file_name, n_data, cut):
        self._tbt_data = pd.read_hdf(file_name, key="df")

        for name, x, y in zip(
                self._tbt_data["BPM_name"], ss._tbt_data["X"],
                ss._tbt_data["Y"]):
            n = min(len(x), n_data)
            if name == self._bpm[0]:
                # Convert from [nm] to [mm].
                self._xy_0 = \
                    np.array([1e-6*x[cut:cut+n], 1e-6*y[cut:cut+n]])
            elif name == self._bpm[1]:
                # Convert from [nm] to [mm].
                self._xy_1 = \
                    np.array([1e-6*x[cut:cut+n], 1e-6*y[cut:cut+n]])
        for k in range(2):
            # Remove average.
            self._xy_0[k] -= np.mean(self._xy_0[k])
            self._xy_1[k] -= np.mean(self._xy_1[k])

    def optics_est(self):
        beta_rel = np.zeros(2)
        dnu = np.zeros(2)
        for k in range(2):
            cov = np.cov(np.array([self._xy_0[k], self._xy_1[k]]))
            print(cov)
            beta_rel[k] = cov[0][0]/cov[1][1]
            dnu[k] = \
                np.arccos(cov[0][1]/(np.sqrt(cov[0][0]*cov[1][1])))/(2e0*np.pi)
        print("\nDefault vs. Estimated:\n"
              +"  beta_rel = [{:5.3f}, {:5.3f}] dnu = [{:5.3f}, {:5.3f}]".
              format(1e0/self._beta_rel[X_], 1e0/self._beta_rel[Y_],
                     self._dnu[X_], self._dnu[Y_]))
        print("  beta_rel = [{:5.3f}, {:5.3f}] dnu = [{:5.3f}, {:5.3f}]".
              format(1e0/beta_rel[X_], 1e0/beta_rel[Y_], dnu[X_], dnu[Y_]))
        self._beta_rel = beta_rel
        self._dnu = dnu
        self._sin = np.sin(2e0*np.pi*self._dnu)
        self._cos = np.cos(2e0*np.pi*self._dnu)

    def state_est(self):
        print("\nstate_est:\n  beta_rel = [{:5.3f}, {:5.3f}]"
              " dnu = [{:5.3f}, {:5.3f}]".
              format(1e0/self._beta_rel[X_], 1e0/self._beta_rel[Y_],
                     self._dnu[X_], self._dnu[Y_]))
        n = len(self._xy_0[X_])
        self._ps_Fl = np.zeros((4, n))
        self._twoJ = np.zeros((2, n))
        for j in range(n):
            for k in range(2):
                self._ps_Fl[2*k][j] = self._xy_0[k][j]/np.sqrt(self._beta[0][k])
                self._ps_Fl[2*k+1][j] = \
                    (np.sqrt(self._beta_rel[k])*self._xy_1[k][j] \
                     -self._xy_0[k][j]*self._cos[k]) \
                    /(np.sqrt(self._beta[0][k])*self._sin[k])
                self._twoJ[k][j] = \
                    self._ps_Fl[2*k][j]**2 + self._ps_Fl[2*k+1][j]**2

    def get_tune(self):
        '''
        Extract tune from turn-by-turn data.
        '''
        fft = fft_class()

        n_data = len(self._xy_0[X_])

        self._tbt_data_fft = np.zeros([2, n_data], dtype="complex")
        A_fft = np.zeros([2, n_data], dtype="float")
        f = np.array(range(n_data))/n_data
        sine_window = sp.signal.windows.cosine(n_data)

        n_peak = 1
        f = np.zeros((2, n_peak), dtype=float)
        A = np.zeros((2, n_peak), dtype=float)
        phi = np.zeros((2, n_peak), dtype=float)
        for k in range(2):
            # Use [mm].
            self._tbt_data_fft[k] = \
                sp.fft.fft(self._xy_0[k]*sine_window)/n_data
            A_fft[k] = abs(self._tbt_data_fft[k])
            f[k], A[k], ind_2 = fft.get_peak_sin(A_fft[k], n_peak)

        nu = np.array([f[X_][0], f[Y_][0]])
        A = np.array([A[X_][0], A[Y_][0]])
        
        return nu, A

    def prt_f(self, fft, n_max, n_peak, nu, f, A, phi):
        print("\nHorizontal Plane")
        print("     f       1-f        A        phi   n_x  n_y    eps")
        for k in range(0, n_peak):
            n_x, n_y, eps = \
                fft.find_harmonic(n_max, nu[X_], nu[Y_], f[X_][k])
            print("  {:7.5f}  {:7.5f}  {:9.3e}  {:6.1f}   {:1d}   {:2d}"
                  "   {:7.1e}".
                  format(f[X_][k], 1e0-f[X_][k], A[X_][k],
                         np.rad2deg(phi[X_][k]), n_x, n_y, eps))
        print("\nVertical Plane")
        print("     f       1-f        A        phi   n_x  n_y    eps")
        for k in range(0, n_peak):
            n_x, n_y, eps = \
                fft.find_harmonic(n_max, nu[X_], nu[Y_], f[Y_][k])
            print("  {:7.5f}  {:7.5f}  {:9.3e}  {:6.1f}   {:1d}   {:2d}"
                  "   {:7.1e}".
                  format(f[Y_][k], 1e0-f[Y_][k], A[Y_][k],
                         np.rad2deg(phi[Y_][k]), n_x, n_y, eps))

    def analyse_twoJ(self, rm_avg, prt, plot):
        fft = fft_class()

        n_max  = 4
        n_peak = 4
        n_data = len(self._twoJ[0])

        if rm_avg:
            for k in range(2):
                # Remove average.
                self._twoJ[k] -= np.mean(self._twoJ[k])

        self._twoJ_fft = np.zeros([2, n_data], dtype="complex")
        A_fft = np.zeros([2, n_data], dtype="float")
        sine_window = sp.signal.windows.cosine(n_data)

        f = np.zeros((2, n_peak), dtype=float)
        A = np.zeros((2, n_peak), dtype=float)
        phi = np.zeros((2, n_peak), dtype=float)
        for k in range(2):
            # Use [mm].
            self._twoJ_fft[k] = \
                sp.fft.fft(self._twoJ[k]*sine_window)/n_data
            A_fft[k] = abs(self._twoJ_fft[k])
            f[k], A[k], ind_2 = fft.get_peak_sin(A_fft[k], n_peak)
            phi[k] = fft.get_phase(ind_2, f[k], self._twoJ[k])

        if prt:
            nu, A_nu = self.get_tune()
            print("\nnu\n     f       1-f        A")
            print("  {:7.5f}  {:7.5f}  {:9.3e}".
                  format(nu[X_], 1e0-nu[X_], A_nu[X_]))
            print("  {:7.5f}  {:7.5f}  {:9.3e}".
                  format(nu[Y_], 1e0-nu[Y_], A_nu[Y_]))
            self.prt_f(fft, n_max, n_peak, nu, f, A, phi)
            
        if plot:
            self.plt_tbt_fft(A_fft)
            plt.show()

#-------------------------------------------------------------------------------

# Main program.

home_dir = os.path.join(os.environ["HOME"])

# Allocate the tbt_class.
ss = state_space()

file_name = \
    home_dir \
    + "/Teresia/Markus/20240216_testset_injectionkickers_storedbeam.hdf"

n_data = 300
cut    = 108

ss.rd_tbt_df(file_name, n_data, cut)

if False:
    ss.optics_est()

ss.state_est()

if not False:
    ss.plt_Floquet_space()

ss.analyse_twoJ(True, True, True)

plt.show()
