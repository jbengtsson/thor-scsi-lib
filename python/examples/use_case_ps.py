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

import gtpsa

import thor_scsi.lib as tslib
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv
from thor_scsi.utils.linear_optics import compute_map, compute_nu_symp, \
    compute_map_and_diag, compute_Twiss_along_lattice, compute_dispersion, \
    compute_nus
from thor_scsi.utils.courant_snyder import compute_A_CS, compute_A, \
    compute_Twiss_A
from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(6, 1)

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
    gr_1.plot(ds.s, ds.twiss.sel(plane='x', par='beta'), label=r'$\beta_x$')
    gr_1.plot(ds.s, ds.twiss.sel(plane='y', par='beta'), label=r'$\beta_y$')
    gr_1.legend()

    gr_2.set_xlabel('s [m]')
    gr_2.set_ylabel(r'$\eta_x$ [m]')
    gr_2.plot(ds.s, ds.dispersion.sel(phase_coordinate='x'), label=r'$\eta_x$')
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


def check_if_stable_one_dim(dim, M):
    # Dim is [0, 1].
    return math.fabs(M[2*dim:2*dim+2, 2*dim:2*dim+2].trace()) < 2e0

def check_if_stable_two_dim(M):
 return check_if_stable_one_dim(0, M) and check_if_stable_one_dim(1, M)


def compute_chromaticity(lat, model_state, eps):
    # If not stable return a large value.
    ksi_1_max = 1e30

    map = compute_map(lat, model_state, delta=-eps, desc=desc)
    M = np.array(map.jacobian())
    if check_if_stable_two_dim(M):
        nu_m = compute_nu_symp(n_dof, M)
    else:
        nu_m = np.array([ksi_1_max, ksi_1_max])

    map = compute_map(lat, model_state, delta=eps, desc=desc)
    M = np.array(map.jacobian())
    if check_if_stable_two_dim(M):
        nu_p = compute_nu_symp(n_dof, M)
    else:
        nu_p = np.array([ksi_1_max, ksi_1_max])

    ksi = (nu_p-nu_m)/(2e0*eps)
    return ksi


def compute_periodic_solution(lat, model_state):
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    M, A = compute_map_and_diag(n_dof, lat, model_state, desc=desc)
    nus = compute_nus(n_dof, M)
    ksi_1 = compute_chromaticity(lat, model_state, 1e-6)
    Twiss = compute_Twiss_A(A)

    print('\nM:\n', mat2txt(M))
    prt_tune_chrom(n_dof, nus, ksi_1)
    prt_Twiss('\nTwiss:\n', Twiss)

    ds = compute_Twiss_along_lattice(n_dof, lat, model_state, A=A, desc=desc)

    return M, A, ds


def prt_tune_chrom(n_dof, nus, ksi):
    if n_dof == 2:
        print(f'\nnu  = [{nus[X_]:7.5f}, {nus[Y_]:7.5f}]')
        print(f'ksi = [{ksi[X_]:7.5f}, {ksi[Y_]:7.5f}]')
    else:
        print(f'\nnu = [{nus[X_]:7.5f}, {nus[Y_]:7.5f}, {nus[Z_]:7.5f}]')


def prt_Twiss(str, Twiss):
    eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    print(str, end='')
    print(f'  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]')
    print(f'  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]')
    print(f'  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]')


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
b_n_zero(lat, sextupole)

# Compute Twiss parameters along lattice.
M, A, data = compute_periodic_solution(lat, model_state)

plt_Twiss(data, 'bessy-ii.png', 'BESSY-II - Linear Optics')

