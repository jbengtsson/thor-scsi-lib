'''Use Case:
     On & off-momentum tune footprintg.
'''
import enum
import logging
# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level='INFO')
logger = logging.getLogger('thor_scsi')


from scipy import optimize

import os
import copy

import math
import numpy as np
import matplotlib.pyplot as plt
import NAFFlib

from scipy.constants import c as c0

import gtpsa

import thor_scsi.lib as tslib
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv
from thor_scsi.utils.linear_optics import compute_map, compute_nus, \
    compute_nu_symp, check_if_stable_2D, compute_nu_xi, compute_map_and_diag, \
    compute_Twiss_along_lattice, compute_dispersion
from thor_scsi.utils.courant_snyder import compute_A_CS, compute_A, \
    compute_Twiss_A
from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

tpsa_order = 2

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(6, tpsa_order)

# Configuration space coordinates.
X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta,
]

[quadrupole, sextupole] = [2, 3]


def plt_Twiss(ds, file_name, title):
    # Turn interactive mode off.
    plt.ioff()

    fig, (gr_1, gr_2) = plt.subplots(2)

    gr_1.set_title(title)
    gr_1.set_xlabel('s [m]')
    gr_1.set_ylabel(r'$\beta_{x,y}$ [m]')
    gr_1.plot(ds.s, ds.twiss.sel(plane='x', par='beta'), 'b-',
              label=r'$\beta_x$')
    gr_1.plot(ds.s, ds.twiss.sel(plane='y', par='beta'), 'r-',
              label=r'$\beta_y$')
    gr_1.legend()

    gr_2.set_xlabel('s [m]')
    gr_2.set_ylabel(r'$\eta_x$ [m]')
    gr_2.plot(ds.s, ds.dispersion.sel(phase_coordinate='x'), label=r'$\eta_x$')

    fig.tight_layout()
    plt.savefig(file_name)
    print(f'File saved as: {file_name:s}')


def plt_nu_A(nu, file_name, title, plane):
    # Turn interactive mode off.
    # plt.ioff()

    fig, ax_1 = plt.subplots()

    ax_1.set_title(title)
    if plane == 0:
        ax_1.set_xlabel(r'$A_{x,y}$ [mm]')
        ax_1.set_ylabel(r'$\nu_x$')
    else:
        ax_1.set_xlabel(r'$A_{x,y}$ [mm]')
        ax_1.set_ylabel(r'$\nu_y$')
    ax_1.plot(1e3*nu[0][0, :], nu[0][1, :], 'b-', label=r'$\nu_x$')
    ax_1.plot(1e3*nu[1][0, :], nu[1][1, :], 'r-', label=r'$\nu_y$')
    ax_1.legend(loc='upper left')

    fig.tight_layout()
    plt.savefig(file_name)
    print(f'File saved as: {file_name:s}')


def compute_eng_units(f0, f_RF, alpha_c, delta, nu):
    '''Convert from physics to engineering units.
    '''
    # delta = np.asarray(delta)
    # f_beta = np.zeros((len(delta), 2), dtype='float')
    f_beta = np.zeros((2,), dtype='float')
    df_RF = -alpha_c * f_RF * delta
    f_beta = nu * np.array(f0)[np.newaxis]
    return df_RF, f_beta


def plt_nu_delta(delta, nu, file_name, title, phys_units):
    # Turn interactive mode off.
    plt.ioff()

    fig, ax_1 = plt.subplots()

    ax_1.set_title(title)
    if phys_units:
        ax_1.set_xlabel(r'$\delta$ [%]')
        ax_1.set_ylabel(r'$\nu_x$')
        ax_1.plot(1e2*delta, nu[:, X_], 'b-', label=r'$\nu_x$')
    else:
        ax_1.set_xlabel(r'$\Delta \rm f_{RF}$ [kHz]')
        ax_1.set_ylabel(r'$\Delta \rm f_{\beta,x}$ [kHz]')
        ax_1.plot(1e-3*delta, 1e-3*nu[:, X_], 'b-', label=r'$\nu_x$')
    ax_1.legend(loc='upper right')
    ax_2 = ax_1.twinx()
    if phys_units:
        ax_2.set_ylabel(r'$\nu_y$')
        ax_2.plot(1e2*delta, nu[:, Y_], 'r-', label=r'$\nu_y$')
    else:
        ax_2.set_ylabel(r'$\Delta \rm f_{\beta,x}$ [kHz]')
        ax_2.plot(1e-3*delta, 1e-3*nu[:, Y_], 'r-', label=r'$\nu_y$')
    ax_2.legend(loc='upper right')
    fig.tight_layout()
    plt.savefig(file_name)
    print(f'File saved as: {file_name:s}')


def get_b_n_elem(lat, fam_name, kid_num, n):
    mp = lat.find(fam_name, kid_num)
    return mp.get_multipoles().get_multipole(n).real

def set_b_n_elem(lat, fam_name, kid_num, n, b_n):
    mp = lat.find(fam_name, kid_num)
    mp.get_multipoles().set_multipole(n, b_n)

def set_b_n_fam(lat, fam_name, n, b_n):
    for mp in lat.elements_with_name(fam_name):
        mp.get_multipoles().set_multipole(n, b_n)


def b_n_zero(lat, n):
    for k in range(len(lat)):
        if n == 2:
            if type(lat[k]) == tslib.Quadrupole:
                lat[k].get_multipoles().set_multipole(quadrupole, 0e0)
                print(lat[k].name)
        elif n == 3:
            if type(lat[k]) == tslib.Sextupole:
                lat[k].get_multipoles().set_multipole(sextupole, 0e0)
        else:
            print(f'\nb_3_zero: undef. multipole order {n:1d}')


def compute_periodic_solution(lat, model_state):
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    M, A = \
        compute_map_and_diag(n_dof, lat, model_state, desc=desc,
                             tpsa_order=tpsa_order)

    stable, nu, xi = compute_nu_xi(desc, tpsa_order, M)
    Twiss = compute_Twiss_A(A)

    print('\nM:\n', mat2txt(M.jacobian()))
    prt_tune_chrom(n_dof, nu, xi)
    prt_Twiss('\nTwiss:\n', Twiss)

    ds = compute_Twiss_along_lattice(n_dof, lat, model_state, A=A, desc=desc)

    return M, A, ds


def prt_tune_chrom(n_dof, nu, xi):
    if n_dof == 2:
        print(f'\nnu = [{nu[X_]:7.5f}, {nu[Y_]:7.5f}]')
        print(f'xi = [{xi[X_]:7.5f}, {xi[Y_]:7.5f}]')
    else:
        print(f'\nnu = [{nu[X_]:7.5f}, {nu[Y_]:7.5f}, {nu[Z_]:7.5f}]')


def prt_Twiss(str, Twiss):
    eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    print(str, end='')
    print(f'  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]')
    print(f'  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]')
    print(f'  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]')


def compute_nu_delta(lat, deltas_to_probe) -> float:
    '''
    Args:
        lat: lattice
        deltas_to_probe: off momenta deltas to check

    Returns:
       :math:`nu` s correspoinding to the deltas to probe
    '''
    n_dof = 2
    nu = np.zeros(deltas_to_probe.shape + (2,), dtype=float)
    prt   = not False

    if prt:
        print('\n#    delta      nu_x     nu_y')

    for k, delta in enumerate(deltas_to_probe):
        map = compute_map(lat, model_state, delta=delta, desc=desc)
        M = np.array(map.jacobian())
        if check_if_stable_2D(M):
            nu[k] = compute_nu_symp(n_dof, M)
        if prt:
            print('  {:10.3e} {:8.5f} {:8.5f}'.format(
                delta, nu[k, X_], nu[k, Y_]))
        else:
            nu[k] = np.array([np.nan, np.nan])

    return nu


def np2tps(a):
    '''Numpy array to ss_vect<double>.
    Todo:
        usse / implement the funcionality in gtspa ss_vbect wrapper
    '''
    dof = 3
    a_tps = tslib.ss_vect_double()
    for k in range(2*dof):
        a_tps[k] = a[k]
    return a_tps


def tps2np(a_tps):
    '''ss_vect<double> to numpy array.
    '''
    dof = 3
    a = np.zeros(2*dof)
    for k in range(2*dof):
        a[k] = a_tps[k]
    return a


def chk_if_lost(ps):
    '''Check if particle is lost during tracking by checking for NAN result in
    phase-space vector.

    Todo:
       rename to check_if_lost
       Identifz used type and see if it uses cst of a gtpsa ss_vect
    '''
    dof = 3
    lost = False

    for k in range(2*dof):
        if math.isnan(ps[k]):
            lost = True
            break

    return lost


def rd_track(file_name):
    '''Read tracking data from file <file_name>.

    todo:
        see if pandas.read_csv does the job?
    '''
    dof = 3
    ps_list = []

    with open(file_name) as f:
        for line in f:
            if line[0] != '#':
                values = line.split()
                n = int(values[0])
                ps_k = map(float, values[1:2*dof+1])
                ps_list.append(np.array(list(ps_k)))
    ps_list = np.array(ps_list)

    return ps_list


def track(lat, model_state, n_turn, ps):
    '''Track particle with initial conditions ps for n turns.

    Returns: (lost, track)

    Todo:
        check output of twiss e.g.
        return an xarraz object for the result?


    '''

    def prt_track(k, ps):
        print(f' {k:4d} {ps[x_]:12.5e} {ps[px_]:12.5e}',
              f'{ps[y_]:12.5e} {ps[py_]:12.5e}',
              f'{ps[delta_]:12.5e} {ps[ct_]:12.5e}', file=outf)

    prt = False

    if prt:
        outf = open('track.out', 'w')
        print('#   n      x            p_x          y            p_y'
              '          delta_       ct_',
              file=outf)
        prt_track(0, ps.cst())
    ps_list = []
    ps_list.append(ps)
    for k in range(n_turn):
        lat.propagate(model_state, ps)
        lost = chk_if_lost(ps)
        if lost:
            break
        ps_list.append(ps.cst())
        if prt:
            prt_track(k+1, ps.cst())
    if prt:
        outf.close()

    return lost, ps_list


def compute_cmplx_ps(ps_list):
    ps_list = np.asarray(ps_list)
    dof = 3
    n = len(ps_list)
    ps_list_cmplx = np.zeros((n, dof), dtype = complex)

    for j in range(n):
        for k in range(dof):
            if k < 2:
                ps_list_cmplx[j, k] = \
                    complex(ps_list[j, 2*k], -ps_list[j,2*k+1])
            else:
                ps_list_cmplx[j, k] = \
                    complex(ps_list[j, 2*k+1], -ps_list[j, 2*k])

    return ps_list_cmplx


def compute_nu(x):
    '''Compute nu from NAFF with complex turn-by-turn data: [x-i*p_x, y-i*p_y].

       Documentation: https://pypi.org/project/NAFFlib
    '''
    # Amplitudes for positive & negative nu.
    nu, A_pos, A_neg = NAFFlib.get_tunes(x, 1)
    A_pos, A_neg = np.absolute([A_pos, A_neg])
    if A_pos < A_neg:
        nu = 1e0 - nu

    return nu


def get_nu(lat, model_state, n_turn, A, delta, eps):
    '''

    '''
    ps = gtpsa.ss_vect_double([A[X_], 0e0, A[Y_], 0e0, delta, 0e0])
    nu = np.zeros(2, dtype=float)

    lost, ps_list = track(lat, model_state, n_turn, ps)
    # ps_list = rd_track('track.out')
    ps_list_cmplx = compute_cmplx_ps(ps_list)
    if not lost:
        for k in range(2):
            nu[k] = compute_nu(ps_list_cmplx[:, k])

    return lost, nu


def validate_amplitude_above_minimum(
        amplitudes: np.ndarray,  minimum_amplitude: float=1e-5,
        strict: bool=False) -> np.ndarray:
    '''Check that amplitudes are not zero

    If not in strict mode, the ones whose absolute value is below the minimum
    amplitude, it is replaced by miminum value
    '''
    chk = np.less(np.absolute(amplitudes), minimum_amplitude)
    if np.sum(chk):
        txt = f'compute amplitude dependent tune shift number of values below',
        f'{minimum_amplitude}:{chk}'
        logger.info(txt)
        logger.debug(f'amplitude values changed {amplitudes[chk]}')
        if strict:
            raise AssertionError(txt)
        amplitudes[chk] = minimum_amplitude
    return amplitudes

def compute_nu_A(lat, model_state, amplitudes: np.ndarray, delta: float, plane,
                 minimum_amplitude: float=0.01e-3, strict: bool=False):
    '''Compute nu_x,y(A_x, A_y, delta).

    Compute amplitude dependent tune shift?
       plane:
         horizontal - 0
         vertical   - 1

       Todo:
           use a sequence of amplitudes instead of A_max and A
           return an xarray object?

           Review how to handle A_min
    '''

    # make to
    validate_amplitude_above_minimum(
        amplitudes, minimum_amplitude=minimum_amplitude, strict=strict)
    A       = np.zeros(2, dtype='float')
    nu_list = []

    prt = not False

    A_min = 0.01e-3
    eps    = 0.01
    n_turn = 511 # n_turn + initial conditions = 512 = 2^9 for FFT.

    if prt:
        print('\n#     A_x        A_y       nu_x     nu_y')
    A[(plane+1) % 2] = A_min
    for k in range(-n_points, n_points+1):
        if k != 0:
            A[plane] = k*A_max[plane]/n_points
        else:
            A[plane] = A_min
        lost, nu = get_nu(lat, model_state, n_turn, A, delta, eps)
        if lost:
            nu[X_], nu[Y_] = [math.nan, math.nan]
            print('#  {:10.3e} {:10.3e}: particle lost', A[X_], A[Y_])
        nu_k = np.array([A[X_], A[Y_], nu[X_], nu[Y_]])
        nu_list.append(nu_k)
        if prt:
            print('  {:10.3e} {:10.3e} {:8.5f} {:8.5f}'.format(
                A[X_], A[Y_], nu[X_], nu[Y_]))
    nu_list = np.array(nu_list)

    # res = xr.DataArray(
    #    data = []
    #    dims = ['amplitude', 'plane'],
    #    coords=[amplitudes, 'x', 'y']
    #
    #)

    return nu_list


t_dir = os.path.join(os.environ['HOME'], 'Nextcloud', 'thor_scsi')
t_file = os.path.join(t_dir, 'b3_sf4375sf_tracy.lat')
# t_dir = os.path.join(os.environ['HOME'], 'git', 'dt4acc', 'lattices')
# t_file = os.path.join(t_dir, 'b2_stduser_beamports_blm_tracy_corr.lat')
# t_file = os.path.join(t_dir, 'BII_JB.lat')

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

# Zero sextupoles.
if False:
    b_n_zero(lat, sextupole)

if not True:
    # Compute Twiss parameters along lattice.
    M, A, data = compute_periodic_solution(lat, model_state)
    plt_Twiss(data, 'bessy-ii_1.png', 'BESSY-II - Linear Optics')

n_points = 15

if True:
    # BESSY III.
    n_x   = 2
    n_y   = 0
    n_sym = 16
    A_max     = np.array([3e-3, 3e-3])
    delta_max = 5e-2
else:
    # BESSY II.
    n_sym = 1
    n_x   = 0
    n_y   = 0
    A_max = np.array([15e-3, 10e-3])
    delta_max = 5e-2

if True:
    amplitudes_x = np.linspace(-A_max[0], A_max[0], n_points)
    amplitudes_y = np.linspace(-A_max[1], A_max[1], n_points)
    nu_A_x = compute_nu_A(lat, model_state, amplitudes_x, 0e0, 0)
    nu_A_y = compute_nu_A(lat, model_state, amplitudes_y, 0e0, 1)
    nu_A_x[:, 2] = n_sym*(nu_A_x[:, 2]+n_x)
    nu_A_y[:, 2] = n_sym*(nu_A_y[:, 2]+n_x)
    nu_A_x[:, 3] = n_sym*(nu_A_x[:, 3]+n_y)
    nu_A_y[:, 3] = n_sym*(nu_A_y[:, 3]+n_y)
    nu_x_A_xy = \
        np.array([nu_A_x[:, 0], nu_A_x[:, 2]]), \
        np.array([nu_A_y[:, 1], nu_A_y[:, 2]])
    nu_y_A_xy = \
        np.array([nu_A_x[:, 0], nu_A_x[:, 3]]), \
        np.array([nu_A_y[:, 1], nu_A_y[:, 3]])
    plt_nu_A(nu_x_A_xy, 'bessy-ii_2.png', 'BESSY-III - $\\nu_x ( A_{x,y} )$', 0)
    plt_nu_A(nu_y_A_xy, 'bessy-ii_3.png', 'BESSY-III - $\\nu_y ( A_{x,y} )$', 1)

if True:
    delta = 5e-2
    delta_to_probe = np.linspace(-delta_max, delta_max, n_points)
    nu = compute_nu_delta(lat, delta_to_probe)
    nu[:, X_] = n_sym*(nu[:, X_]+n_x)
    nu[:, Y_] = n_sym*(nu[:, Y_]+n_y)

    phys_units = True
    if phys_units:
        plt_nu_delta(
            delta_to_probe, nu,
            'bessy-ii_4.png', 'BESSY-III - $\\nu_{x,y} ( \delta )$', True)
    else:
        # BESSY II.
        circ = 240.0
        f_RF = 499.6366302e6
        alpha_c = 7.038e-4
        f0 = c0/circ
        df_RF, f_beta = \
            compute_eng_units(f0, f_RF, alpha_c, delta_to_probe, nu)

        plt_nu_delta(
            df_RF, f_beta,
            'bessy-ii_4.png', 'BESSY-II - $\\nu_{x,y} ( \\Delta \\rm f_{RF} )$',
            False)

if not False:
    plt.show()
