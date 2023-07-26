"""Use Case:
     Sextupole response matrix analysis.

     J. Bengtsson 
     𝑇ℎ𝑒 𝑆𝑒𝑥𝑡𝑢𝑝𝑜𝑙𝑒 𝑆𝑐ℎ𝑒𝑚𝑒 𝑓𝑜𝑟 𝑡ℎ𝑒 𝑆𝑤𝑖𝑠𝑠 𝐿𝑖𝑔ℎ𝑡 𝑆𝑜𝑢𝑟𝑐𝑒 (𝑆𝐿𝑆) – 𝐴𝑛 𝐴𝑛𝑎𝑙𝑦𝑡𝑖𝑐 𝐴𝑝𝑝𝑟𝑜𝑎𝑐ℎ
     See SLS Tech Note 09/97

     https://ados.web.psi.ch/slsnotes/sls0997.pdf
"""

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
from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
     radiate as rad


# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt


# Configuration space coordinates.
X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]


class MultipoleIndex(enum.IntEnum):
    quadrupole = 2
    sextupole  = 3


def prt_Twiss(str, Twiss):
    """
    """
    eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    print(str, end="")
    print(f"  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
    print(f"  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
    print(f"  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")


def compute_periodic_solution(lat, model_state, named_index, desc, no):
    """
    """
    # Deegrees of freedom.
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
        lo.compute_Twiss_along_lattice(
            n_dof, lat, model_state, A=A_map, desc=desc, mapping=named_index)

    return M, A, ds


def get_phi_elem(lat, fam_name, n_kid):
    elem = lat.find(fam_name, n_kid)
    return elem.get_length() * elem.get_curvature() * 180e0 / np.pi


def compute_phi(lat):
    """Compute total bend angle.
    """
    prt = False
    phi = 0e0
    for k in range(len(lat)):
        if type(lat[k]) == tslib.Bending:
            dphi = get_phi_elem(lat, lat[k].name, 0)
            phi += dphi
            if prt:
                print("{:8s} {:5.3f} {:6.3f}".
                      format(lat[k].name, lat[k].get_length(), dphi))
    return phi


def print_map(str, map):
    """
    # 
    """
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
    energy = 2.5e9
    model_state.Energy    = energy
    model_state.radiation = False
    model_state.Cavity_on = False

    print("-> read_lattice")
    return n_dof, lat, model_state


def set_b_n_elem(lat, fam_name, kid_num, n, b_n):
    mp = lat.find(fam_name, kid_num)
    mp.get_multipoles().set_multipole(n, b_n)


def set_b_n_fam(lat, fam_name, n, b_n):
    for mp in lat.elements_with_name(fam_name):
        mp.get_multipoles().set_multipole(n, b_n)


corresponding_types = {
    tslib.Sextupole:         tslib.SextupoleTpsa,
    tslib.Quadrupole:        tslib.QuadrupoleTpsa,
    tslib.HorizontalSteerer: tslib.HorizontalSteererTpsa,
    tslib.VerticalSteerer:   tslib.VerticalSteererTpsa,
}


def convert_magnet_to_knobbable(a_magnet: tslib.Mpole) -> tslib.MpoleTpsa:
    config = a_magnet.config()
    corresponding_type = corresponding_types[type(a_magnet)]
    return corresponding_type(config)


def mult_prm(elems, mult_family, n, desc):
    print("\nmult_prm ->")
    for k in range(len(mult_family)):
        index = mult_family[k].index
        elem = convert_magnet_to_knobbable(elems[index])
        elems[index] = \
            knobs.make_magnet_knobbable(
                elem, po=1, desc=desc, named_index=named_index,
                multipole_number=n, offset=False
            )
        print("\nmult_prm:", elems[index])
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
    print("-> mult_prm\n")
    return elems


def lat_mult_prm(mult_prm_name, lat, n, desc):
    elems = [elem for elem in lat]
    mult_family = lat.elements_with_name(mult_prm_name)
    elems = mult_prm(elems, mult_family, n, desc)
    return tslib.AcceleratorTpsa(elems)


# Work-around for C++ virtual function -> Python mapping issue.
# See function above:
#   mult_prm
def propagate(lat, model_state, desc, no, nv, named_index):
    M = gtpsa.ss_vect_tpsa(desc, no, nv, index_mapping=named_index)
    M.set_identity()
    for k in range(len(lat)):
        lat[k].propagate(model_state, M)
    return M


# Number of phase-space coordinates.
nv = 6
# Variables max order.
no = 3
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

t_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "BESSY-III",
    "ipac_2023")
t_file = os.path.join(t_dir, "b3_cf425cf_thor_scsi.lat")

n_dof, lat, model_state = read_lattice(t_file)

if False:
    # Zero sextopoles.
    print("\nZeroing sextupoles.")
    set_b_n_fam(lat, "sf", MultipoleIndex.sextupole, 0e0)
    set_b_n_fam(lat, "sd", MultipoleIndex.sextupole, 0e0)
    set_b_n_fam(lat, "sd2", MultipoleIndex.sextupole, 0e0)

if False:
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    desc = gtpsa.desc(nv, no, nv_prm, no_prm)
    M, A, data = \
        compute_periodic_solution(lat, model_state, named_index, desc, no)

    # Compute linear chromaticity.
    M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=no)
    stable, nu, xi = lo.compute_nu_xi(desc, no, M)

    C = rad.compute_circ(lat)

    print("\n  C [m]     = {:9.7f}".format(C))
    print("  phi [deg] = {:5.3f}".format(compute_phi(lat)))
    print("  nu        = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))
    print("  xi        = [{:7.5f}, {:7.5f}]".format(xi[X_], xi[Y_]))

nv_prm = 1
no_prm = 1

named_index = gtpsa.IndexMapping(
    dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6)
)

desc = gtpsa.desc(nv, no, nv_prm, no_prm)

lat_ptc = lat_mult_prm("uq3", lat, 2, desc)
M = propagate(lat_ptc, model_state, desc, no, nv, named_index)

print("\nM:\n", M)
if not False:
    print_map("\nM:", M)

if False:
    M2 = M
    # Busted.
    M2.getOrder(M2, 2)
    print("\n:\n", M2)

if False:
    M_inv = tslib.inv(M)

    M_M_inv = gtpsa.ss_vect_tpsa(desc, no, nv, index_mapping=named_index)
    M_M_inv = tslib.compose(M, M_inv)
    print("\nM:", M)
    print("\nM^-1:", M_inv)
    print("\nM*M^-1:", M_M_inv)
    # print("\nM*M^-1:", M_M_inv[0])

# assert False

h = tslib.M_to_h_DF(M)
print("\nh:")
h.print()
