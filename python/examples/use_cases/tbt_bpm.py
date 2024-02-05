
import os
import numpy as np
import scipy as sp
import NAFFlib
import matplotlib.pyplot as plt

from thor_scsi.utils.fft import get_peak_sin, get_phase


def rd_tbt(file_name):
    with open(file_name) as f:
        for line in f:
            tbt_data = line.split()
            # Unit is microns.
            tbt_data = 1e-6*np.array(tbt_data, dtype='float')
    return tbt_data


def rd_tbt_hor_and_ver(file_name_hor, file_name_ver):
    tbt_data_hor = rd_tbt(file_name_hor)
    tbt_data_ver = rd_tbt(file_name_ver)
    return np.array([tbt_data_hor, tbt_data_ver])


def plt_tbt(tbt_data, nu, A, plane):
    # Turn interactive mode off.
    # plt.ioff()

    n = len(tbt_data)

    fig, (gr_1, gr_2) = plt.subplots(2)

    if plane == 0:
        file_name = "tbt_hor.png"
        gr_1.set_title("Horizontal Position")
        gr_2.set_title("FFT of Horizontal Position")
        color = "b-"
    else:
        file_name = "tbt_ver.png"
        gr_1.set_title("Vertical Position")
        gr_2.set_title("FFT of Horizontal Position")
        color = "r-"

    gr_1.set_xlabel("Turn Number")
    gr_1.set_ylabel(r"x [$\mu$m]")
    gr_1.stem(tbt_data, linefmt=color, markerfmt="", basefmt="k-")

    gr_2.set_xlabel("f")
    gr_2.set_ylabel(r"A_x [$\mu$m]")
    gr_2.stem \
        (nu[0:n//2+1], A[0:n//2+1], linefmt=color, markerfmt="", basefmt="k-")

    fig.tight_layout()

    plt.savefig(file_name)
    print(f"File saved as: {file_name:s}")


def get_f_naff(x):
    '''
    Extract f and amplitude with NAFF lib.
    However, somehow, the phase is not being provided (Sigh!).

    (By maximising the Fourier integral numerically:
    J. Bengtsson, Y. Hidaka
    𝑁𝑆𝐿𝑆-𝐼𝐼: 𝑇𝑢𝑟𝑛-𝑏𝑦-𝑇𝑢𝑟𝑛 𝐵𝑃𝑀 𝐷𝑎𝑡𝑎 𝐴𝑛𝑎𝑙𝑦𝑠𝑖𝑠 – 𝐴 𝑈𝑠𝑒 𝐶𝑎𝑠𝑒 𝐴𝑝𝑝𝑟𝑜𝑎𝑐ℎ
    NSLSII-ASD-TN-125 (2014)

    https://doi.org/10.2172/1480956)

    Documentation: https://pypi.org/project/NAFFlib.
    '''
    nu, A_pos, A_neg = NAFFlib.get_tunes(x, 1)
    A_pos, A_neg = np.absolute([A_pos, A_neg])
    return nu, A_pos


home_dir = os.path.join(os.environ["HOME"])
file_name_hor = home_dir + "/Teresia/20240122_TbT_hor.dat"
file_name_ver = home_dir + "/Teresia/20240122_TbT_ver.dat"

tbt_data = rd_tbt_hor_and_ver(file_name_hor, file_name_ver)
# print(tbt_data)
# Skip initial transient.
cut, n_data = 9999, 2**9
tbt_data = np.array([tbt_data[0, cut:cut+n_data], tbt_data[1, cut:cut+n_data]])
for k in range(2):
    tbt_data[k] -= np.mean(tbt_data[k])

tbt_data_fft = np.zeros([2, n_data], dtype="complex")
A_fft = np.zeros([2, n_data], dtype="float")
f = np.array(range(n_data))/n_data
sine_window = sp.signal.windows.cosine(n_data)
print()
for j in range(2):
    tbt_data_fft[j] = sp.fft.fft(tbt_data[j]*sine_window)/n_data
    A_fft[j] = abs(tbt_data_fft[j])
    nu, A, k = get_peak_sin(A_fft[j], 1)
    phi = get_phase(k, nu[0], tbt_data[j])
    print("nu = {:8.6f} A = {:9.3e} phi = {:5.1f}".
          format(nu[0], A[0], phi*180e0/sp.pi))

if True:
    for j in range(2):
        plt_tbt(tbt_data[j], f, A_fft[j], j)

    plt.show()
