
import os.path

import numpy as np

import gtpsa
from thor_scsi.factory import accelerator_from_config
import thor_scsi.lib as tslib
from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils.phase_space_vector import ss_vect_tps2ps_jac
from thor_scsi.utils.output import mat2txt

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


def ind_0():
    ind = (np.array([0, 0, 0, 0, 0, 0]))
    return ind


def ind_1(k):
    ind = (np.array([0, 0, 0, 0, 0, 0]))
    ind[k] = 1
    return ind


def ind_2(j, k):
    ind = (np.array([0, 0, 0, 0, 0, 0]))
    ind[j] = 1
    ind[k] = 1
    return ind


def compute_nu_ksi(M):
    nu   = np.zeros(2)
    ksi  = np.zeros(2)
    m_11 = gtpsa.tpsa(desc, tpsa_order)
    m_22 = gtpsa.tpsa(desc, tpsa_order)
    tr   = gtpsa.tpsa(desc, tpsa_order)

    for k in range(2):
        m_11.set(ind_0(), 0e0, M[2*k].get(ind_1(2*k)))
        m_11.set(ind_1(delta_), 0e0, M[2*k].get(ind_2(2*k, delta_)))
        m_22.set(ind_0(), 0e0, M[2*k+1].get(ind_1(2*k+1)))
        m_22.set(ind_1(delta_), 0e0, M[2*k+1].get(ind_2(2*k+1, delta_)))

        tr = m_11 + m_22
        print('\ntr:\n', tr)

        if True:
            nu[k] = np.arccos(tr.get(ind_0())/2e0)/(2e0*np.pi)
            if M.jacobian()[2*k][2*k+1] < 0e0:
                nu[k] = 1e0 - nu[k]
            ksi[k] = \
                -tr.get(ind_1(delta_))/(4e0*np.pi*np.sin(2e0*np.pi*nu[k]))
        else:
            nu = acos(tr/2e0)/(2e0*np.pi)
            print('\nnu:\n', nu)

    return nu, ksi


def prt_nu_ksi(M):
    print('\nM:\n', mat2txt(M.jacobian()[0:2, 0:2]))
    nu, ksi = compute_nu_ksi(M)
    print('\n  nu  = [{:7.5f}, {:7.5f}]'.format(nu[X_], nu[Y_]))
    print('  ksi = [{:7.5}, {:7.5}]'.format(ksi[X_], ksi[Y_]))
    return ksi


t_dir = os.path.join(os.environ['HOME'], 'git', 'dt4acc', 'lattices')
t_file = os.path.join(t_dir, 'b2_stduser_beamports_blm_tracy_corr.lat')

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

# Zero sextupoles.
if not False:
    b_n_zero(lat, sextupole)

M = gtpsa.ss_vect_tpsa(desc, tpsa_order)
M.set_identity()
lat.propagate(model_state, M)

prt_nu_ksi(M)

exit()

cst = M.cst()
J = M.jacobian()
H = M.hessian()

# print('\ncst:\n', cst)
# print('\nJ:\n', J)
# print('\nH:\n', H)

x      = M[0]
order  = 2
dof    = 3
ps_dim = 2*dof

indices = np.zeros(ps_dim, int)
indices = [1, 0, 0, 0, 1, 0]
r = x.index(indices)
assert r >= 0
print('\n', x.get(indices))

scl = 0.0
x.set(indices, scl, 355.0/113.0)
print('\n', x.get(indices))
