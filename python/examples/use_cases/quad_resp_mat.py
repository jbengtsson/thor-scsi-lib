import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

from scipy import optimize

import copy as _copy
from dataclasses import dataclass
import math
import os
from typing import Sequence
from typing import Tuple

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from thor_scsi.pyflame import Config

import gtpsa
import thor_scsi.lib as tslib

from thor_scsi.utils import knobs
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.accelerator import instrument_with_standard_observers
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv

from thor_scsi.utils.extract_info import accelerator_info
from thor_scsi.utils import linear_optics as lo, courant_snyder as cs

# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt


# Configuration space coordinates.
X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]


def prt_Twiss(str, Twiss):
    """

    todo:
        rename str e.g. header? prefix?

    """
    # eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    # that way I also check that Twiss has exactly three parameters
    eta, alpha, beta = Twiss
    print(str, end="")
    print(f"  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
    print(f"  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
    print(f"  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")


def compute_periodic_solution(lat, model_state, named_index, desc, no):
    """
    Todo:
        model_state: rename to calculation_configuration or calc_config
    """
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    stable, M, A = lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc,
                                           tpsa_order=no)
    print("\nM:\n", M)
    res = cs.compute_Twiss_A(A)
    Twiss = res[:3]
    prt_Twiss("\nTwiss:\n", Twiss)
    A_map = gtpsa.ss_vect_tpsa(desc, no)
    A_map.set_jacobian(A)
    ds = \
        lo.compute_Twiss_along_lattice \
        (n_dof, lat, model_state, A=A_map, desc=desc, mapping=named_index)

    return M, A, ds


def prt_map(str, map):
    print(str)
    map.x.print()
    map.px.print()
    map.y.print()
    map.py.print()
    map.delta.print()
    map.ct.print()


def read_lattice(t_file):
    # Read in & parse lattice file.
    print("\nread_lattice ->")
    lat = accelerator_from_config(t_file)

    # Set lattice state (Rf cavity on/off, etc.)
    model_state = tslib.ConfigType()

    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    print("-> read_lattice")
    return n_dof, lat, model_state


corresponding_types = {
    tslib.Quadrupole:        tslib.QuadrupoleTpsa,
    tslib.HorizontalSteerer: tslib.HorizontalSteererTpsa,
    tslib.VerticalSteerer:   tslib.VerticalSteererTpsa,
}


def convert_magnet_to_knobbable(a_magnet: tslib.Mpole) -> tslib.MpoleTpsa:
    config = a_magnet.config()
    corresponding_type = corresponding_types[type(a_magnet)]
    return corresponding_type(config)


def multipole_prm(elems, mult_family, n, desc, named_index):
    print("\nmultipole_prm ->")
    for k in range(len(mult_family)):
        index = mult_family[k].index
        print(mult_family[k].name, mult_family[k].index)
        elem = convert_magnet_to_knobbable(elems[index])
        elems[index] = \
            knobs.make_magnet_knobbable(
                elem, po=1, desc=desc, named_index=named_index,
                multipole_number=n, offset=True
            )
        # While the RHS pointer can be recasted to:
        #   CellVoid
        #   ElemTypeKnobbed
        #   QuadrupoleType
        # the assignment of the LHS only gives:
        #   CellVoid
        # and the lack of:
        #   ElemtypeKnobbed
        # will lead to an assert on line 158 in:
        #   thor_scsi/std_machine/accelerator.cc
        #
    print("-> multipole_prm\n")
    return elems


def prt_lat(lat):
    print("\nprt_lat:")
    for k in range(0, 3):
        print("\n", lat[k].name, "\n", lat[k], "\n", lat[k].config)


# Number of phase-space coordinates.
nv = 6
# Variables max order.
no = 2
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB")
t_file = os.path.join(t_dir, "b3_sfsf4Q_tracy_jb_3.lat")

n_dof, lat, model_state = read_lattice(t_file)

if False:
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    desc = gtpsa.desc(nv, no, nv_prm, no_prm)
    M, A, data = \
        compute_periodic_solution(lat, model_state, named_index, desc, no)

nv_prm = 3
no_prm = 1

named_index = gtpsa.IndexMapping(
    dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6, dx=7, dy=8)
)

desc = gtpsa.desc(nv, no, nv_prm, no_prm)

elems = [elem for elem in lat]
mult_family = lat.elements_with_name("uq4")
elems = multipole_prm(elems, mult_family, 2, desc, named_index)
lat_tpsa = tslib.AcceleratorTpsa(elems)

if False:
    prt_lat(lat)
    prt_lat(lat_tpsa)
    assert False

# Bug in lattice propagator.
if not False:
    M = gtpsa.ss_vect_tpsa(desc, no, nv, index_mapping=named_index)
    M.set_identity()

    index1 = mult_family[0].index
    lat_tpsa.propagate(model_state, M, 0, index1)
    lat_tpsa[index1].propagate(model_state, M)
    index2 = mult_family[1].index
    lat_tpsa.propagate(model_state, M, index1+1, index2-index1-1)
    lat_tpsa[index2].propagate(model_state, M)
    lat_tpsa.propagate(model_state, M, index2+1, len(lat_tpsa)-index2-1)

    print("\nM:\n", M)

    if False:
        prt_map("\nM:", M)

if False:
    M = gtpsa.ss_vect_tpsa(desc, no, nv, index_mapping=named_index)
    M.set_identity()
    lat.propagate(model_state, M)
    print()

    M = gtpsa.ss_vect_tpsa(desc, no, nv, index_mapping=named_index)
    M.set_identity()
    lat_tpsa.propagate(model_state, M)
    print("\nM:\n", M)

    if not False:
        prt_map("\nM:", M)
