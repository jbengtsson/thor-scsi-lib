"""FFT - Module for Elementary Signal Processing
Originally Pascal code for my thesis work for LEAR, CERN - machine translated
to C with P2C.

References:

E. Asséo, J. Bengtsson, M. Chanel 𝐿𝐸𝐴𝑅 𝐵𝑒𝑎𝑚 𝑆𝑡𝑎𝑏𝑖𝑙𝑖𝑡𝑦 𝐼𝑚𝑝𝑟𝑜𝑣𝑒𝑚𝑒𝑛𝑡𝑠 𝑈𝑠𝑖𝑛𝑔 𝐹𝐹𝑇
𝐴𝑛𝑎𝑙𝑦𝑠𝑖𝑠 EPAC 1988.

https://accelconf.web.cern.ch/e88/PDF/EPAC1988_0541.PDF

J. Bengtsson 𝑁𝑜𝑛-𝐿𝑖𝑛𝑒𝑎𝑟 𝑇𝑟𝑎𝑛𝑠𝑣𝑒𝑟𝑠𝑒 𝐷𝑦𝑛𝑎𝑚𝑖𝑐𝑠 𝑓𝑜𝑟 𝑆𝑡𝑜𝑟𝑎𝑔𝑒 𝑅𝑖𝑛𝑔𝑠 𝑤𝑖𝑡ℎ 𝐴𝑝𝑝𝑙𝑖𝑐𝑎𝑡𝑖𝑜𝑛𝑠 𝑡𝑜 𝑡ℎ𝑒
𝐿𝑜𝑤 𝐸𝑛𝑒𝑟𝑔𝑦 𝐴𝑛𝑡𝑖𝑝𝑟𝑜𝑡𝑜𝑛 𝑅𝑖𝑛𝑔 (𝐿𝐸𝐴𝑅) 𝑎𝑡 𝐶𝐸𝑅𝑁 CERN-88-05.

http://dx.doi.org/10.5170/CERN-1988-005
"""

import numpy as np
# NumPy has also an implementation of FFT but the SciPy is more efficient.
import scipy as sp

import NAFFlib


def get_ind(n, k):
    prt = False
    if k == 1:
        ind_1 = ind_3 = 2
    elif k == n//2+1:
        ind_1 = ind_3 = n//2
    else:
        ind_1 = k - 1
        ind_3 = k + 1
    if prt:
        print("get_ind: n = {:d} k = {:d} ind_1 = {:d} ind_3 = {:d}".
              format(n, k, ind_1, ind_3))
    return ind_1, ind_3


def get_ind_cmplx(n, k, ind_1, ind_3):
    prt = False
    if k == 1:
        ind_1 = ind_3 = 2
    elif k == n:
        ind_1 = ind_3 = n - 1
    else:
        ind_1 = k - 1
        ind_3 = k + 1
    if prt:
        print("get_ind_cmplx: n = {:d} k = {:d} ind_1 = {:d} ind_3 = {:d}".
              format(n, k, ind_1, ind_3))
    return ind_1, ind_3


def get_peak(x):
    '''
    Locate peak in FFT spectrum.
    '''
    peak = 0e0
    k = 1;
    n = len(x)
    for ind_2 in range(1, n//2+2):
        ind_1, ind_3 = get_ind(n, ind_2)
        if (x[ind_2-1] > peak) and (x[ind_1-1] < x[ind_2-1]) \
           and (x[ind_3-1] < x[ind_2-1]):
            peak = x[ind_2-1]
            k = ind_2
    return k


def get_peak_cmplx(x):
    '''
    Locate peak in FFT spectrum.
    '''
    peak = 0e0
    k = 1
    n = len(x)
    for ind_2 in range(1, n+1):
        ind_1, ind_3 = get_ind_cmplx(n, ind_2)
        if (x[ind_2-1] > peak) and (x[ind_1-1] < x[ind_2-1]) \
           and (x[ind_3-1] < x[ind_2-1]):
            peak = x[ind_2-1]
            k = ind_2
    return k


def interpol_sin_nu(x, k):
    '''
    Extract frequency by 2-point nonlinear interpolation with sine window:

            1              2 A(k)       1
       nu = - [ k - 1 + ------------- - - ] ,    k-1 <= N*nu <= k
            N           A(k-1) + A(k)   2
    '''
    n = len(x)
    ind_1, ind_3 = get_ind(n, k)
    if x[ind_3 - 1] > x[ind_1 - 1]:
        ampl_1 = x[k - 1]
        ampl_2 = x[ind_3 - 1]
        ind = k
    else:
        ampl_1 = x[ind_1 - 1]
        ampl_2 = x[k - 1]
        # Interpolate in right direction for 0 frequency.
        if k != 1:
            ind = ind_1
        else:
            ind = 0
    # Avoid division by zero.
    if ampl_1+ampl_2 != 0e0:
        return (ind-1+2e0*ampl_2/(ampl_1+ampl_2)-0.5e0)/n
    else:
        return 0e0


def interpol_sin_nu_cmplx(x, k):
    '''
    Extract frequency by 2-point nonlinear interpolation with sine window:

           1              2 A(k)       1
      nu = - [ k - 1 + ------------- - - ] ,    k-1 <= N*nu <= k
           N           A(k-1) + A(k)   2
    '''
    n = len(x)
    ind_1, ind_3 = get_ind_cmplx(n, k)
    if x[ind_3 - 1] > x[ind_1 - 1]:
        ampl_1 = x[k - 1]
        ampl_2 = x[ind_3 - 1]
        ind = k
    else:
        ampl_1 = x[ind_1 - 1]
        ampl_2 = x[k - 1]
        # Interpolate in right direction for 0 frequency.
        if k != 1:
            ind = ind_1
        else:
            ind = 0
    # Avoid division by zero.
    if ampl_1+ampl_2 != 0e0:
        return (ind-1+2*ampl_2/(ampl_1+ampl_2)-0.5e0)/n
    else:
        return 0e0


def interpol_sin_ampl(x, nu, k):
    '''
    Extract amplitude by nonlinear interpolation for sine window.
    The distribution is given by:

              1    sin pi ( k + 1/2 )     sin pi ( k - 1/2 )
      F(k) =  - ( -------------------- + -------------------- )
              2      pi ( k + 1/2 )          pi ( k - 1/2 )
    '''
    n = len(x)
    corr = \
        (np.sinc((k-1e0+0.5e0-nu*n))+np.sinc((k-1e0-0.5e0-nu*n)))/2e0
    return x[k-1]/corr


def interpol_sin_ampl_cmplx(x, nu, k):
    '''
    Extract amplitude by nonlinear interpolation for sine window.
    The distribution is given by:

              1    sin pi ( k + 1/2 )     sin pi ( k - 1/2 )
      F(k) =  - ( -------------------- + -------------------- )
              2      pi ( k + 1/2 )          pi ( k - 1/2 )
    '''
    n = len(x)
    corr = \
        (sinc(sp.pi*(k-1e0+0.5e0-nu*n))+sinc(sp.pi*(k-1e0-0.5e0-nu*n)))/2e0
    return x[k-1]/corr


def get_phase(k, nu, x):
    '''
    Extract phase by linear interpolation for rectangular window with:
      -pi <= phi <= pi
    '''
    n = len(x)
    x_fft = sp.fft.fft(x)
    phi = np.arctan2(x_fft[k-1].imag, x_fft[k-1].real) - (n*nu-k+1e0)*sp.pi
    if phi > sp.pi:
        phi -= 2e0*sp.pi
    elif phi < -sp.pi:
        phi += 2e0*sp.pi
    return phi


def find_harmonic_eps(nu_x, nu_y, f, eps):
    n = len(x)
    for j in range(n):
        for k in range(-n, n+1):
            delta = abs(j*nu_x+k*nu_y)
            delta -= math.trunc(delta)
            if delta > 0.5e0:
                delta = 1 - delta
            delta = abs(delta-f)
            delta -= math.trunc(delta)
            if delta > 0.5e0:
                delta = 1 - delta
            if delta < eps:
                if (abs(j) + abs(k) < n) and (j != 0 or k >= 0):
                    found = True
                    n_x = j
                    n_y = k
    return n_x, n_y


def find_harmonic(nu_x, nu_y, f):
    '''
    Match f by a linear combination of nu_x and nu_y.
    '''
    found = False
    eps = 0.5e-6
    while True:
        eps *= 1e1
        find_harmonic_eps(nu_x, nu_y, f)
        if found:
            break
    return n_x, n_y


def get_peak_sin(x, n_peaks):
    n = len(x)
    x1 = x
    nu = np.zeros(n_peaks, dtype="float")
    A = np.zeros(n_peaks, dtype="float")
    for j in range(n_peaks):
        k = get_peak(x1)
        nu[j] = interpol_sin_nu(x1, k)
        A[j] = interpol_sin_ampl(x1, nu[j], k)
        # Flatten peak to enable new call.
        ind_1, ind_3 = get_ind(n, k)
        if x1[ind_1-1] > x1[ind_3-1]:
            x1[k-1] = x1[ind_1-1]
        else:
            x1[k-1] = x1[ind_3-1]
    return nu, A, k


def get_peak_sin_cmplx(x, n_peaks):
    n = len(x)
    x1 = x
    nu = np.zeros(n_peaks, dtype="float")
    A = np.zeros(n_peaks, dtype="float")
    for j in range(n_peaks):
        k = get_peak_cmplx(n, x1)
        nu[j] = interpol_sin_nu_cmplx(x1, k)
        A[j] = interpol_sin_ampl_cmplx(x1, nu[j], k)
        # Flatten peak to enable new call.
        ind_1, ind_3 = get_ind_cmplx(n, k)
        if x1[ind_1-1] > x1[ind_3-1]:
            x1[k-1] = x1[ind_1-1]
        else:
            x1[k-1] = x1[ind_3-1]
    return nu, A, k


def get_f_naff(x):
    '''
    Extract f & amplitude from turn-by-turn BPM data with NAFF-lib.
    However, somehow, the phase is not being provided (Sigh!).

    (By maximising the Fourier integral by a numerical search:
    J. Bengtsson, Y. Hidaka
    𝑁𝑆𝐿𝑆-𝐼𝐼: 𝑇𝑢𝑟𝑛-𝑏𝑦-𝑇𝑢𝑟𝑛 𝐵𝑃𝑀 𝐷𝑎𝑡𝑎 𝐴𝑛𝑎𝑙𝑦𝑠𝑖𝑠 – 𝐴 𝑈𝑠𝑒 𝐶𝑎𝑠𝑒 𝐴𝑝𝑝𝑟𝑜𝑎𝑐ℎ
    NSLSII-ASD-TN-125 (2014)

    https://doi.org/10.2172/1480956)

    Documentation: https://pypi.org/project/NAFFlib.
    '''
    f, A_pos, A_neg = NAFFlib.get_tunes(x, 1)
    A_pos, A_neg = np.absolute([A_pos, A_neg])
    # f & A_pos are arrays.
    return f, A_pos


__all__ = ["get_peak_sin", "get_peak_sin_cmplx", "get_phase", "get_f_naff"]
