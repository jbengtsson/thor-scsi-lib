"""Use Case:
     Sextupole strength calibration.
"""

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


# BESSY-II Sextupole Families.
b_3_s1 = [
    's1mt1rl', 's1mt1rr', 's1md2rl', 's1md2rr', 's1mt2rl', 's1mt2rr',
    's1md3rl', 's1md3rr', 's1mt3rl', 's1mt3rr', 's1md4rl', 's1md4rr',
    's1mt4rl', 's1mt4rr', 's1md5rl', 's1md5rr', 's1mt5rl', 's1mt5rr',
    's1md6rl', 's1md6rr', 's1mt6rl', 's1mt6rr', 's1md7rl', 's1md7rr',
    's1mt7rl', 's1mt7rr', 's1md8rl', 's1md8rr', 's1mt8rl', 's1mt8rr',
    's1md1rl', 's1md1rr'
]

b_3_s2 = [
    's2m2d1rl', 's2m2d1rr', 's2m1t1rl', 's2m1t1rr', 's2m2t1rl', 's2m2t1rr',
    's2m1d2rl', 's2m1d2rr', 's2m2d2rl', 's2m2d2rr', 's2m1t2rl', 's2m1t2rr',
    's2m2t2rl', 's2m2t2rr', 's2m1d3rl', 's2m1d3rr', 's2m2d3rl', 's2m2d3rr',
    's2m1t3rl', 's2m1t3rr', 's2m2t3rl', 's2m2t3rr', 's2m1d4rl', 's2m1d4rr',
    's2m2d4rl', 's2m2d4rr', 's2m1t4rl', 's2m1t4rr', 's2m2t4rl', 's2m2t4rr',
    's2m1d5rl', 's2m1d5rr', 's2m2d5rl', 's2m2d5rr', 's2m1t5rl', 's2m1t5rr',
    's2m2t5rl', 's2m2t5rr', 's2m1d6rl', 's2m1d6rr', 's2m2d6rl', 's2m2d6rr',
    's2m1t6rl', 's2m1t6rr', 's2m2t6rl', 's2m2t6rr', 's2m1d7rl', 's2m1d7rr',
    's2m2d7rl', 's2m2d7rr', 's2m1t7rl', 's2m1t7rr', 's2m2t7rl', 's2m2t7rr',
    's2m1d8rl', 's2m1d8rr', 's2m2d8rl', 's2m2d8rr', 's2m1t8rl', 's2m1t8rr',
    's2m2t8rl', 's2m2t8rr', 's2m1d1rl', 's2m1d1rr'
]

b_3_s3_1 = [
    's3m1d2rl', 's3m1d2rr', 's3m2d2rl', 's3m2d2rr', 's3m1d3rl', 's3m1d3rr',
    's3m2d3rl', 's3m2d3rr', 's3m1d4rl', 's3m1d4rr', 's3m2d4rl', 's3m2d4rr',
    's3m1d5rl', 's3m1d5rr', 's3m2d5rl', 's3m2d5rr', 's3m1d6rl', 's3m1d6rr',
    's3m2d6rl', 's3m2d6rr', 's3m1d7rl', 's3m1d7rr', 's3m2d7rl', 's3m2d7rr',
    's3m1d8rl', 's3m1d8rr', 's3m2d8rl', 's3m2d8rr'
]

b_3_s3_2 = [
    's3m2d1rl', 's3m2d1rr', 's3m1d1rl', 's3m1d1rr'
]

b_3_s3_3 = [
    's3m1t1rl', 's3m1t1rr', 's3m2t1rl', 's3m2t1rr', 's3m1t2rl', 's3m1t2rr',
    's3m2t2rl', 's3m2t2rr', 's3m1t3rl', 's3m1t3rr', 's3m2t3rl', 's3m2t3rr',
    's3m1t4rl', 's3m1t4rr', 's3m2t4rl', 's3m2t4rr', 's3m1t5rl', 's3m1t5rr',
    's3m2t5rl', 's3m2t5rr', 's3m1t6rl', 's3m1t6rr', 's3m2t6rl', 's3m2t6rr',
    's3m1t7rl', 's3m1t7rr', 's3m2t7rl', 's3m2t7rr', 's3m1t8rl', 's3m1t8rr',
    's3m2t8rl', 's3m2t8rr'
]

b_3_s4_1 = [
    's4m1d2rl', 's4m1d2rr', 's4m2d2rl', 's4m2d2rr', 's4m1d3rl', 's4m1d3rr',
    's4m2d3rl', 's4m2d3rr', 's4m1d4rl', 's4m1d4rr', 's4m2d4rl', 's4m2d4rr',
    's4m1d5rl', 's4m1d5rr', 's4m2d5rl', 's4m2d5rr', 's4m1d6rl', 's4m1d6rr',
    's4m2d6rl', 's4m2d6rr', 's4m1d7rl', 's4m1d7rr', 's4m2d7rl', 's4m2d7rr',
    's4m1d8rl', 's4m1d8rr', 's4m2d8rl', 's4m2d8rr'
]

b_3_s4_2 = [
    's4m2d1rl', 's4m2d1rr', 's4m1d1rl', 's4m1d1rr'
]

b_3_s4_3 = [
    's4m1t1rl', 's4m1t1rr', 's4m2t1rl', 's4m2t1rr', 's4m1t2rl', 's4m1t2rr',
    's4m2t2rl', 's4m2t2rr', 's4m1t3rl', 's4m1t3rr', 's4m2t3rl', 's4m2t3rr',
    's4m1t4rl', 's4m1t4rr', 's4m2t4rl', 's4m2t4rr', 's4m1t5rl', 's4m1t5rr',
    's4m2t5rl', 's4m2t5rr', 's4m1t7rl', 's4m1t7rr', 's4m2t7rl', 's4m2t7rr',
    's4m1t8rl', 's4m1t8rr', 's4m2t8rl', 's4m2t8rr'
]

b_3_s4_4 = [
    's4m1t6rl', 's4m1t6rr', 's4m2t6rl', 's4m2t6rr'
]


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


def compute_eng_units(f0, f_RF, alpha_c, delta, nu):
    """Convert from physics to engineering units.
    """
    f_beta = np.zeros((len(delta), 2), dtype='float')
    df_RF = -alpha_c*f_RF*delta
    for k in range(2):
        f_beta[:, k] = nu[:, k]*f0
    return df_RF, f_beta


def plt_nu_vs_delta(delta, nu, file_name, title, phys_units):
    # Turn interactive mode off.
    plt.ioff()

    fig, ax_1 = plt.subplots()

    ax_1.set_title(title)
    if phys_units:
        ax_1.set_xlabel(r'$\delta$ [%]')
        ax_1.set_ylabel(r'$\nu_x$ [m]')
        ax_1.plot(1e2*delta, nu[:, X_], 'b-', label=r'$\nu_x$')
    else:
        ax_1.set_xlabel(r'$\Delta \rm f_{RF}$ [kHz]')
        ax_1.set_ylabel(r'$\Delta \rm f_{\beta,x}$ [kHz]')
        ax_1.plot(1e-3*delta, 1e-3*nu[:, X_], 'b-', label=r'$\nu_x$')
    ax_1.legend(loc='upper left')
    ax_2 = ax_1.twinx()
    if phys_units:
        ax_2.set_ylabel(r'$\nu_y$ [m]')
        ax_2.plot(1e2*delta, nu[:, Y_], 'r-', label=r'$\nu_y$')
    else:
        ax_2.set_ylabel(r'$\Delta \rm f_{\beta,x}$ [kHz]')
        ax_2.plot(1e-3*delta, 1e-3*nu[:, Y_], 'r-', label=r'$\nu_y$')
    ax_2.legend(loc='upper right')
    fig.tight_layout()
    plt.savefig(file_name)


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


def compute_dnu_ddelta(lat, n, delta_max):
    n_dof = 2
    delta = np.zeros(2*n+1, dtype='float')
    nu = np.zeros((2*n+1, 2), dtype='float')
    for k in range(-n, n+1):
        delta[k+n] = k*delta_max/n
        map = compute_map(lat, model_state, delta=delta[k+n], desc=desc)
        M = np.array(map.jacobian())
        if check_if_stable_2D(M):
            nu[k+n] = compute_nu_symp(n_dof, M)
        else:
            nu[k+n] = np.array([np.nan, np.nan])
    return delta, nu


t_dir = os.path.join(os.environ['HOME'], 'git', 'dt4acc', 'lattices')
t_file = os.path.join(t_dir, 'b2_stduser_beamports_blm_tracy_corr.lat')

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

# Compute Twiss parameters along lattice.
M, A, data = compute_periodic_solution(lat, model_state)

plt_Twiss(data, 'bessy-ii_1.png', 'BESSY-II - Linear Optics')

delta, nu = compute_dnu_ddelta(lat, 10, 5e-2)

phys_units = not True
if phys_units:
    plt_nu_vs_delta(
        delta, nu, 'bessy-ii_2.png', 'BESSY-II - $\\nu_{x,y} ( \delta )$', True)
else:
    circ = 240.0
    f_RF = 499.6366302e6
    alpha_c = 7.038e-4
    f0 = c0/circ
    df_RF, f_beta = compute_eng_units(f0, f_RF, alpha_c, delta, nu)

    plt_nu_vs_delta(
        df_RF, f_beta, 'bessy-ii_2.png', 'BESSY-II - $\\nu_{x,y} ( \delta )$',
        not True)
if not False:
    plt.show()
