"""
Author:

  Johan Bengtsson
  02/02/24

FFT - Class for Elementary Signal Processing
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

from math import trunc
import numpy as np
# NumPy has also an implementation of FFT but the SciPy is more efficient.
import scipy as sp

naff_lib = False

if naff_lib:
    import NAFFlib


class fft_class:
    # Private

    def __init__(self):
        pass

    # Public.

    def get_phase(self, ind_2, nu, x):
        '''
        Extract phase by linear interpolation for rectangular window with:
          -pi <= phi <= pi
        '''
        n = len(x)
        n_peak = len(ind_2)
        phi = np.zeros(n_peak)
        x_fft = sp.fft.fft(x)
        for k in range(n_peak):
            phi[k] = \
                np.arctan2(x_fft[ind_2[k]-1].imag, x_fft[ind_2[k]-1].real) \
                - (n*nu[k]-ind_2[k]+1e0)*sp.pi
            if phi[k] > sp.pi:
                phi[k] -= 2e0*sp.pi
            elif phi[k] < -sp.pi:
                phi[k] += 2e0*sp.pi
        return phi


    def find_harmonic(self, n, nu_x, nu_y, f):
        '''
        Match f by a linear combination of nu_x and nu_y.
        '''
        found = False
        eps = 0.5e-6
        while True:
            eps *= 1e1
            found, n_x, n_y, delta = find_harmonic_eps(n, nu_x, nu_y, f, eps)
            if found:
                break
        return n_x, n_y, delta


    def get_peak_sin(self, x, n_peaks):
        n = len(x)
        x1 = x
        ind_2 = np.zeros(n_peaks, dtype=int)
        nu = np.zeros(n_peaks, dtype=float)
        A = np.zeros(n_peaks, dtype=float)
        for k in range(n_peaks):
            ind_2[k] = get_peak(x1)
            nu[k] = interpol_sin_nu(x1, ind_2[k])
            A[k] = interpol_sin_ampl(x1, nu[k], ind_2[k])
            # Flatten peak for a recursive approach - i.e., the function can be
            # called again to obtain the next, etc.
            ind_1, ind_3 = get_ind(n, ind_2[k])
            if x1[ind_1-1] > x1[ind_3-1]:
                x1[ind_2[k]-1] = x1[ind_1-1]
            else:
                x1[ind_2[k]-1] = x1[ind_3-1]
        return nu, A, ind_2


    def get_peak_sin_cmplx(self, x, n_peaks):
        n = len(x)
        x1 = x
        ind_2 = np.zeros(n_peaks, dtype=int)
        nu = np.zeros(n_peaks, dtype="float")
        A = np.zeros(n_peaks, dtype="float")
        for k in range(n_peaks):
            ind_2[k] = get_peak_cmplx(n, x1)
            nu[k] = interpol_sin_nu_cmplx(x1, ind_2)
            A[k] = interpol_sin_ampl_cmplx(x1, nu[k], ind_2[k])
            # Flatten peak to enable new call.
            ind_1, ind_3 = get_ind_cmplx(n, ind_2)
            if x1[ind_1-1] > x1[ind_3-1]:
                x1[ind_2-1] = x1[ind_1-1]
            else:
                x1[ind_2-1] = x1[ind_3-1]
        return nu, A, ind_2

# ------------------------------------------------------------------------------

def get_ind(n, ind_2):
    if ind_2 == 1:
        ind_1 = ind_3 = 2
    elif ind_2 == n//2+1:
        ind_1 = ind_3 = n//2
    else:
        ind_1 = ind_2 - 1
        ind_3 = ind_2 + 1
    return ind_1, ind_3


def get_ind_cmplx(n, ind_2):
    prt = False
    if ind_2 == 1:
        ind_1 = ind_3 = 2
    elif ind_2 == n:
        ind_1 = ind_3 = n - 1
    else:
        ind_1 = ind_2 - 1
        ind_3 = ind_2 + 1
    if prt:
        print("get_ind_cmplx: n = {:d} ind_2 = {:d} ind_1 = {:d} ind_3 = {:d}".
              format(n, ind_2, ind_1, ind_3))
    return ind_1, ind_3


def get_peak(x):
    '''
    Locate peak in FFT spectrum.
    '''
    peak = 0e0
    ind_2 = 1;
    n = len(x)
    for k in range(1, n//2+2):
        ind_1, ind_3 = get_ind(n, k)
        if (x[k-1] > peak) and (x[ind_1-1] < x[k-1]) \
           and (x[ind_3-1] < x[k-1]):
            peak = x[k-1]
            ind_2 = k
    return ind_2


def get_peak_cmplx(x):
    '''
    Locate peak in FFT spectrum.
    '''
    peak = 0e0
    ind_2 = 1
    n = len(x)
    for k in range(1, n+1):
        ind_1, ind_3 = get_ind_cmplx(n, k)
        if (x[k-1] > peak) and (x[ind_1-1] < x[k-1]) \
           and (x[ind_3-1] < x[k-1]):
            peak = x[k-1]
            ind_2 = k
    return ind_2


def interpol_sin_nu(x, ind_2):
    '''
    Extract frequency by 2-point nonlinear interpolation with sine window:

            1                     2 A(ind_2)        1
       nu = - [ ind_2 - 1 + --------------------- - - ] ,
            N               A(ind_2-1) + A(ind_2)   2

            ind_2-1 <= N*nu <= ind_2
    '''
    n = len(x)
    ind_1, ind_3 = get_ind(n, ind_2)
    if x[ind_3 - 1] > x[ind_1 - 1]:
        ampl_1 = x[ind_2 - 1]
        ampl_2 = x[ind_3 - 1]
        ind = ind_2
    else:
        ampl_1 = x[ind_1 - 1]
        ampl_2 = x[ind_2 - 1]
        # Interpolate in right direction for 0 frequency.
        if ind_2 != 1:
            ind = ind_1
        else:
            ind = 0
    # Avoid division by zero.
    if ampl_1+ampl_2 != 0e0:
        return (ind-1+2e0*ampl_2/(ampl_1+ampl_2)-0.5e0)/n
    else:
        return 0e0


def interpol_sin_nu_cmplx(x, ind_2):
    '''
    Extract frequency by 2-point nonlinear interpolation with sine window:

           1                    2 A(ind_2)         1
      nu = - [ ind_2 - 1 + --------------------- - - ] ,
           N               A(ind_2-1) + A(ind_2)   2

           ind_2-1 <= N*nu <= ind_2
    '''
    n = len(x)
    ind_1, ind_3 = get_ind_cmplx(n, ind_2)
    if x[ind_3 - 1] > x[ind_1 - 1]:
        ampl_1 = x[ind_2 - 1]
        ampl_2 = x[ind_3 - 1]
        ind = ind_2
    else:
        ampl_1 = x[ind_1 - 1]
        ampl_2 = x[ind_2 - 1]
        # Interpolate in right direction for 0 frequency.
        if ind_2 != 1:
            ind = ind_1
        else:
            ind = 0
    # Avoid division by zero.
    if ampl_1+ampl_2 != 0e0:
        return (ind-1+2*ampl_2/(ampl_1+ampl_2)-0.5e0)/n
    else:
        return 0e0


def interpol_sin_ampl(x, nu, ind_2):
    '''
    Extract amplitude by nonlinear interpolation for sine window.
    The distribution is given by:

                 1       sin pi ( ind_2 + 1/2 )   sin pi ( ind_2 - 1/2 )
      F(ind_2) = - ( ---------------------- + ---------------------- )
                 2         pi ( ind_2 + 1/2 )       pi ( ind_2 - 1/2 )
    '''
    n = len(x)
    corr = \
        (np.sinc((ind_2-1e0+0.5e0-nu*n))+np.sinc((ind_2-1e0-0.5e0-nu*n)))/2e0
    return x[ind_2-1]/corr


def interpol_sin_ampl_cmplx(x, nu, ind_2):
    '''
    Extract amplitude by nonlinear interpolation for sine window.
    The distribution is given by:

                 1   sin pi ( ind_2 + 1/2 )   sin pi ( ind_2 - 1/2 )
      F(ind_2) = - ( ---------------------- + ---------------------- )
                 2     pi ( ind_2 + 1/2 )       pi ( ind_2 - 1/2 )
    '''
    n = len(x)
    corr = \
        (sinc(sp.pi*(ind_2-1e0+0.5e0-nu*n))
         +sinc(sp.pi*(ind_2-1e0-0.5e0-nu*n)))/2e0
    return x[ind_2-1]/corr


def find_harmonic_eps(n, nu_x, nu_y, f, eps):
    prt = False
    found = False
    n_x = n_y = 0
    delta_min = eps
    if prt:
        print()
    for j in range(n+1):
        for k in range(-n, n+1):
            delta = abs(j*nu_x+k*nu_y)
            delta -= int(delta)
            if delta > 0.5e0:
                delta = 1 - delta
            delta = abs(delta-f)
            delta -= int(delta)
            if delta > 0.5e0:
                delta = 1 - delta
            if prt:
                print("find_harmonic_eps:"
                      +" {:9.3e} {:1d} {:2d} {:7.5f} {:7.5f} {:8.1e} {:8.1e}".
                      format(f, j, k, nu_x, nu_y, delta, eps))
            if (delta < delta_min) and (abs(j)+abs(k) <= n) \
               and (j != 0 or k >= 0):
                found = True
                delta_min = min(delta, delta_min)
                n_x = j
                n_y = k
    return found, n_x, n_y, delta_min

# ------------------------------------------------------------------------------

def get_f_naff(x):
    '''
    Extract f & amplitude from turn-by-turn BPM data with NAFF-lib.
    However, somehow, the phase is not being provided (Sigh!).

    (By maximising the Fourier integral by a numerical search see:
    J. Bengtsson, Y. Hidaka
    𝑁𝑆𝐿𝑆-𝐼𝐼: 𝑇𝑢𝑟𝑛-𝑏𝑦-𝑇𝑢𝑟𝑛 𝐵𝑃𝑀 𝐷𝑎𝑡𝑎 𝐴𝑛𝑎𝑙𝑦𝑠𝑖𝑠 – 𝐴 𝑈𝑠𝑒 𝐶𝑎𝑠𝑒 𝐴𝑝𝑝𝑟𝑜𝑎𝑐ℎ
    NSLSII-ASD-TN-125 (2014)

    https://doi.org/10.2172/1480956
    )

    Documentation:

    https://pypi.org/project/NAFFlib
    '''
    if not naff_lib:
        print("\nget_f_naff: enable import NAFF-lib")
    f, A_pos, A_neg = NAFFlib.get_tunes(x, 1)
    A_pos, A_neg = np.absolute([A_pos, A_neg])
    # f & A_pos are arrays.
    return f, A_pos

__all__ = ["fft_class", "get_f_naff"]
