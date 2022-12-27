
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
    """Zero all the n-multipoles in the lattice.
    """
    for k in range(len(lat)):
        if n == quadrupole:
            if type(lat[k]) == tslib.Quadrupole:
                lat[k].get_multipoles().set_multipole(quadrupole, 0e0)
                print(lat[k].name)
        elif n == sextupole:
            if type(lat[k]) == tslib.Sextupole:
                lat[k].get_multipoles().set_multipole(sextupole, 0e0)
        else:
            print(f'\nb_3_zero: undef. multipole order {n:1d}')


def ind_0():
    """Index for constant TPSA term.
    """
    return np.array([0, 0, 0, 0, 0, 0])


def ind_1(k):
    """Index for linear TPSA term.
    """
    ind = (np.array([0, 0, 0, 0, 0, 0]))
    ind[k] = 1
    return ind


def ind_2(j, k):
    """Index for quadratic TPSA term.
    """
    ind = (np.array([0, 0, 0, 0, 0, 0]))
    ind[j] = 1
    ind[k] = 1
    return ind


def compute_nu_xi(M):
    """Compute the tune & linear chromaticity from the trace of the Poincaré
       map:
          nu + xi * delta = arccos( Trace{M} / 2 ) / = ( 2 * pi )
    """
    nu  = np.zeros(2)
    xi = np.zeros(2)
    for k in range(2):
        tr = gtpsa.tpsa(desc, tpsa_order)
        # m_11 + delta * m_16.
        tr.set(ind_0(), 1e0, M[2*k].get(ind_1(2*k)))
        tr.set(ind_1(delta_), 1e0, M[2*k].get(ind_2(2*k, delta_)))
        # m_22 + delta * m_26.
        tr.set(ind_0(), 1e0, M[2*k+1].get(ind_1(2*k+1)))
        tr.set(ind_1(delta_), 1e0, M[2*k+1].get(ind_2(2*k+1, delta_)))
        nu_tpsa = gtpsa.acos(tr/2e0)/(2e0*np.pi)
        if M.jacobian()[2*k][2*k+1] < 0e0:
            nu_tpsa = 1e0 - nu_tpsa
        nu[k] = nu_tpsa.get(ind_0())
        xi[k] = nu_tpsa.get(ind_1(delta_))
    return nu, xi


def prt_nu_xi(M):
    nu, xi = compute_nu_xi(M)
    print('\nM:\n', mat2txt(M.jacobian()))
    print('\n  nu = [{:7.5f}, {:7.5f}]'.format(nu[X_], nu[Y_]))
    print('  xi = [{:7.5}, {:7.5}]'.format(xi[X_], xi[Y_]))
    return nu, xi


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

prt_nu_xi(M)

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
