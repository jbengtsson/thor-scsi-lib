"""Use Case:
     Super period straight section matching with a triplet.
   Constraints:
    [alpha_x,y, beta_x,y, eta_x, eta'_x] at the centre of the straigth.
   Parameters:
     triplet gradients (3) & spacing (3).
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
    compute_map_and_diag, compute_Twiss_along_lattice, compute_dispersion
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
    print('\nM:\n', mat2txt(M))
    Twiss = compute_Twiss_A(A)
    prt_Twiss('\nTwiss:\n', Twiss)

    ds = compute_Twiss_along_lattice(n_dof, lat, model_state, A=A, desc=desc)

    return M, A, ds


def prt_Twiss(str, Twiss):
    eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    print(str, end='')
    print(f'  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]')
    print(f'  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]')
    print(f'  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]')


def get_b_n_elem(lat, fam_name, kid_num, n):
    mp = lat.find(fam_name, kid_num)
    return mp.get_multipoles().get_multipole(n).real

def set_b_n_elem(lat, fam_name, kid_num, n, b_n):
    mp = lat.find(fam_name, kid_num)
    mp.get_multipoles().set_multipole(n, b_n)

def set_b_n_fam(lat, fam_name, n, b_n):
    for mp in lat.elements_with_name(fam_name):
        mp.get_multipoles().set_multipole(n, b_n)

def get_L_elem(lat, fam_name, n_kid):
    elem = lat.find(fam_name, n_kid)
    return elem.get_length()

def set_L_fam(lat, fam_name, L):
    for elem in lat.elements_with_name(fam_name):
        elem.set_length(L)


def get_b_2_elem(lat, fam_name, kid_num):
    return get_b_n_elem(lat, fam_name, kid_num, quadrupole)

def set_b_2_fam(lat, fam_name, b_2):
    set_b_n_fam(lat, fam_name, quadrupole, b_2)


# Global variables.
n_iter    = 0
chi_2_min = 1e30


def get_prm(lat, prm_list):
    # Dictionary of parameter types and corresponding get functions.
    how_to_get_prm = {
        'b_2' : get_b_2_elem,
        'L'   : get_L_elem
    }

    prms = []
    for k in range(len(prm_list)):
        prms.append(how_to_get_prm[prm_list[k][1]](lat, prm_list[k][0], 0))
    return prms


def set_prm(lat, prm_list, prms):
    # Dictionary of parameter types and corresponding set functions.
    how_to_set_prm = {
        'b_2' : set_b_2_fam,
        'L'   : set_L_fam
    }

    for k in range(len(prm_list)):
        how_to_set_prm[prm_list[k][1]](lat, prm_list[k][0], prms[k])


def compute_chi_2(lat, model_state, A0, chi_2_min, prms):
    # Local copy.
    A1 = copy.copy(A0)
    lat.propagate(model_state, A1, loc, len(lat)-loc)

    chi_2 = 0e0
    Twiss_k = compute_Twiss_A(A1.jacobian())[0:3]
    for k in range(3):
        chi_2 += weight[k]*np.sum((Twiss_k[k]-Twiss1_design[k])**2)

    if weight[3] != 0e0:
        ksi_1 = compute_chromaticity(lat, model_state, 1e-6)
        dchi_2 = weight[3]*np.sum(ksi_1**2)
        chi_2 += dchi_2

    if chi_2 < chi_2_min:
        chi_2_min = chi_2
        print(f'\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n\n  prms   ='
              + vec2txt(prms))
        prt_Twiss('\n', Twiss_k)
        if weight[3] != 0e0:
            print(f'\n  ksi_1  = [{ksi_1[X_]:5.3f}, {ksi_1[Y_]:5.3f}]')
            print(f'  dchi_2 = {dchi_2:9.3e}')
    return chi_2, chi_2_min


def match_straight(lat, loc, prm_list, bounds, Twiss0, Twiss1, Twiss1_design,
                   weight):

    def f_match(prms):
        global n_iter
        global chi_2_min

        n_iter += 1
        set_prm(lat, prm_list, prms)
        chi_2, chi_2_min = compute_chi_2(lat, model_state, A0, chi_2_min, prms)

        return chi_2

    def prt_result(f_match, prms0, minimum):
        global n_iter
        global chi_2_min

        # Compute new Twiss parameters along lattice.
        M, A, data = compute_periodic_solution(lat, model_state)
        plt_Twiss(data, 'after.png', 'Linear Optics - After')

        print('\nInitial parameter values:')
        n_iter = 0
        chi_2_min = 1e30
        f_match(prms0)
        print('\nFinal parameter values:')
        n_iter = 0
        chi_2_min = 1e30
        f_match(minimum['x'])
        print('\n Minimum:\n', minimum)


    global n_iter
    global chi_2_min

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    print('\nmatch_straight:\n\nloc = ', loc)
    prt_Twiss('\nTwiss parameters at centre of super period:\n', Twiss0)
    prt_Twiss('\nTwiss parameters at centre of straight:\n',
              Twiss1)
    prt_Twiss('\nDesired Twiss parameters at the centre of straight:\n',
              Twiss1_design)

    A0 = gtpsa.ss_vect_tpsa(desc, 1)
    A0.set_zero()
    A0.set_jacobian(compute_A(*Twiss0))

    # Initialise parameters.
    prms1 = prms0 = get_prm(lat, prm_list)
    print('\nprms = ', prms1)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.
    minimum = \
        optimize.minimize(
            f_match, prms1, method='Powell', bounds=bounds,
            options={'xtol':x_tol, 'ftol':f_tol, 'maxiter':max_iter})

    prt_result(f_match, prms0, minimum)


def get_Twiss(loc):
    eta = \
        np.array([data.dispersion.sel(phase_coordinate='x').values[loc],
                  data.dispersion.sel(phase_coordinate='px').values[loc],
                  data.dispersion.sel(phase_coordinate='y').values[loc],
                  data.dispersion.sel(phase_coordinate='py').values[loc]])
    alpha = \
        np.array([data.twiss.sel(plane='x', par='alpha').values[loc],
                  data.twiss.sel(plane='y', par='alpha').values[loc]])
    beta = \
        np.array([data.twiss.sel(plane='x', par='beta').values[loc],
                  data.twiss.sel(plane='y', par='beta').values[loc]])
    return eta, alpha, beta


t_dir = os.path.join(os.environ['HOME'], 'Nextcloud', 'thor_scsi')
# t_file = os.path.join(t_dir, 'b3_tst.lat')
t_file = os.path.join(t_dir, 'b3_sf_40Grad_JB.lat')

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

# Compute Twiss parameters along lattice.
M, A, data = compute_periodic_solution(lat, model_state)
plt_Twiss(data, 'before.png', 'Linear Optics - Before')

# Get dispersion & Twiss parameters at centre of super period.
loc = lat.find('sf_m', 3).index
print(f'\nCentre of super period: {loc:1d} {lat[loc].name:s}')
Twiss0 = get_Twiss(loc)

# Compute linear chromaticity.
ksi_1 = compute_chromaticity(lat, model_state, 1e-6)
print(ksi_1)
print(f'\nksi = [{ksi_1[X_]:5.3f}, {ksi_1[Y_]:5.3f}]')

# Twiss parameters at centre of straight.
Twiss1 = get_Twiss(len(lat)-1)

# Desired Twiss parameters at centre of straight.
eta   = np.array([0e0, 0e0, 0e0, 0e0])
alpha = np.array([0e0, 0e0])
beta  = np.array([3.0, 3.0])
Twiss1_design = eta, alpha, beta

weight = np.array([
    1e0,   # Eta.
    1e3,   # Alpha.
    1e0,   # Beta.
    0*1e-4   # ksi_1.
])

# Triplet parameter family names & type.
prm_list = [
    ('uq1', 'b_2'),
    ('uq2', 'b_2'),
    ('uq3', 'b_2'),
    ('ul1', 'L'),
    ('ul2', 'L'),
    ('ul3', 'L')
]

# Max parameter range.
b_2_max = 10.0
L_min   = 0.1

bounds = [
    (-b_2_max, 0.0),
    (0.0,      b_2_max),
    (-b_2_max, 0.0),
    (L_min,    0.25),
    (L_min,    0.3),
    (L_min,    0.25)
]

# Zero sextopoles.
# Compute linear chromaticity.
set_b_n_fam(lat, 'sf_h', sextupole, 0e0)
set_b_n_fam(lat, 'sd_h', sextupole, 0e0)
ksi_1 = compute_chromaticity(lat, model_state, 1e-6)
print(f'\nksi = [{ksi_1[X_]:5.3f}, {ksi_1[Y_]:5.3f}]')

match_straight(lat, loc, prm_list, bounds, Twiss0, Twiss1, Twiss1_design,
               weight)

if not False:
    plt.show()
print('\nPlots saved as: before.png & after.png')
