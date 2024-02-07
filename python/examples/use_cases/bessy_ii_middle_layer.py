'''
Module to convert from Engineering to Physics units & vice versa,
and I/O - get/put - get/set to the EPICS control system.
'''


import os
import numpy as np

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
    radiate as rad, closed_orbit as co

from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt


# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

# Configuration space coordinates.
X_, Y_, Z_ = [
    ts.spatial_index.X,
    ts.spatial_index.Y,
    ts.spatial_index.Z
]
# Phase-space coordinates.
[x_, px_, y_, py_, ct_, delta_] = [
    ts.phase_space_index_internal.x,
    ts.phase_space_index_internal.px,
    ts.phase_space_index_internal.y,
    ts.phase_space_index_internal.py,
    ts.phase_space_index_internal.ct,
    ts.phase_space_index_internal.delta,
]


def rd_sext_conv_fact(file_name):
    # Unit is [1/A] for thick sextupoles.
    b3L = {}
    with open(file_name) as f:
        for line in f:
            if line[0] != "#":
                token = line.split()
                b3L.update({token[0] : float(token[1])})
    return b3L


def get_pv(param_name):
    return pv


def set_pv(param_name, pv):
    return


def read_lattice(t_file):
    # Read in & parse lattice file.
    lat = accelerator_from_config(t_file)

    # Set lattice state (Rf cavity on/off, etc.)
    model_state = ts.ConfigType()

    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    return n_dof, lat, model_state


def print_Twiss_param(str, Twiss):
    # eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    # that way I also check that Twiss has exactly three parameters
    eta, alpha, beta = Twiss
    print(str, end="")
    print(f"  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
    print(f"  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
    print(f"  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")


def compute_periodic_solution(lat, model_state, named_index, desc):
    """
    Todo:
        model_state: rename to calculation_configuration or calc_config
    """
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    stable, M, A = lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc)
    print("\nM:\n" + mat2txt(M.jacobian()[:6, :6]))
    res = cs.compute_Twiss_A(A)
    Twiss = res[:3]
    print_Twiss_param("\nTwiss:\n", Twiss)
    A_map = gtpsa.ss_vect_tpsa(desc, no)
    A_map.set_jacobian(A)
    ds = \
        lo.compute_Twiss_along_lattice(
            n_dof, lat, model_state, A=A_map, desc=desc, mapping=named_index)

    return M, A, ds


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 1
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

# Beam Energy [eV].
E_0 = 1.7e9

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
home_dir = os.path.join(os.environ["HOME"])
file_name_lat = \
    os.path.join \
    (home_dir, "Nextcloud", "thor_scsi", "JB", "BESSY-II",
     "b2_stduser_beamports_blm_tracy_corr.lat")
file_name_conv_coeff = \
    os.path.join(home_dir, "Teresia/conversion-factors-sextupoles.txt")

b3L = rd_sext_conv_fact(file_name_conv_coeff)

n_dof, lat, model_state = read_lattice(file_name_lat)

M, A, data = \
    compute_periodic_solution(lat, model_state, named_index, desc)
