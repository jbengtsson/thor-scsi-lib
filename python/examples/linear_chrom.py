
import os.path

import numpy as np

import gtpsa
from thor_scsi.factory import accelerator_from_config
import thor_scsi.lib as tslib
from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils.linear_optics import compute_map, compute_nu_xi
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


def prt_nu_xi(desc, tpsa_order, M):
    stable, nu, xi = compute_nu_xi(desc, tpsa_order, M)
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

M = compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
prt_nu_xi(desc, tpsa_order, M)

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
