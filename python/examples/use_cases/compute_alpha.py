"""Use Case:
     Compute the momentum compaction, alpha, to arbitrary order.
"""


import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="ERROR")
logger = logging.getLogger("thor_scsi")


import math
import os

import copy

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
    radiate as rad, closed_orbit as co

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


def read_lattice(t_file):
    # Read in & parse lattice file.
    lat = accelerator_from_config(t_file)

    # Set lattice state (Rf cavity on/off, etc.)
    model_state = ts.ConfigType()

    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    return n_dof, lat, model_state


def plt_Twiss(ds, file_name, title):
    # Turn interactive mode off.
    plt.ioff()

    fig, (gr_1, gr_2) = plt.subplots(2)

    gr_1.set_title(title)
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel(r"$\beta_{x,y}$ [m]")
    gr_1.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), "b-",
              label=r"$\beta_x$")
    gr_1.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), "r-",
              label=r"$\beta_y$")
    gr_1.legend()

    gr_2.set_xlabel("s [m]")
    gr_2.set_ylabel(r"$\eta_x$ [m]")
    gr_2.plot(ds.s, ds.dispersion.sel(phase_coordinate="x"), label=r"$\eta_x$")
    fig.tight_layout()
    plt.savefig(file_name)


def print_Twiss(str, Twiss):
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
    print("\nM:", M)
    res = cs.compute_Twiss_A(A)
    Twiss = res[:3]
    print_Twiss("Twiss:\n", Twiss)
    A_map = gtpsa.ss_vect_tpsa(desc, no)
    A_map.set_jacobian(A)
    ds = \
        lo.compute_Twiss_along_lattice(
            n_dof, lat, model_state, A=A_map, desc=desc, mapping=named_index)

    return M, A, ds


def compute_disp(map, no, nv, desc):
    n = 4

    ind = np.zeros(nv, dtype=int)
    Id = gtpsa.ss_vect_tpsa(desc, 1)
    D = gtpsa.ss_vect_double(0e0)
    M = gtpsa.ss_vect_tpsa(desc, 1)

    Id.set_identity()
    ind[delta_] = 1
    for k in range(n):
        D.iloc[k] = map.iloc[k].get(ind)

    M.set_jacobian(map.jacobian())
    M.x.set(ind, 0e0, 0e0)
    M.px.set(ind, 0e0, 0e0)
    M.ct.set(ind, 0e0, 0e0)
    ind[delta_] = 0
    ind[x_] = 1
    M.ct.set(ind, 0e0, 0e0)
    ind[x_] = 0
    ind[px_] = 1
    M.ct.set(ind, 0e0, 0e0)
    space_sel = np.array([1, 1, 1, 1, 0, 0, 0], dtype=int)
    tmp = Id - M
    tmp.delta = Id.delta
    tmp.ct = Id.ct
    tmp.iloc[6] = Id.iloc[6]
    print("\ntmp:", tmp)
    print("\ninv:", ts.inv(tmp))
    assert False

    print("\nD:\n", D)
    print("map:", map)
    print("M:\n", M)

    # id_mat = np.identity(n)
    # disp = np.linalg.inv(id_mat - M[:n, :n]) @ D[:n]

    # print("\neta:\n", mat2txt(disp))

    # disp = np.append(disp, [1e0, 0e0, 0e0])
    # print("eta:\n", mat2txt(map.jacobian() @ disp))


def print_cod(str, cod):
    n_dof = 3
    ind = np.zeros(nv, int)
    print(str, end="")
    for k in range(2*n_dof):
        print(" {:10.3e}".format(cod.iloc[k]), end="")
    print()

    
def print_map(str, map):
    n_dof = 3
    print(str)
    for k in range(2*n_dof):
        map.iloc[k].print()


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 2
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB",
                     "BESSY-III", "ipac_2023")
t_file = os.path.join(t_dir, "b3_cf425cf_thor_scsi.lat")

print("\nlattice file:  \n", t_file)

n_dof, lat, model_state = read_lattice(t_file)

M, A, data = \
    compute_periodic_solution(lat, model_state, named_index, desc)

if not True:
    plt_Twiss(data, "lin_opt.png", "Linear Optics")
    plt.show()
    print("\nPlots saved as: before.png & after.png")

r = co.compute_closed_orbit(lat, model_state, delta=0e0, eps=1e-10, desc=desc)

map = gtpsa.ss_vect_tpsa(desc, no)
map.set_identity()
r.x0.ct = 0e0
map += r.x0
lat.propagate(model_state, map)
print("\nmap:", map)

if not True:
    print("\nmap[ct]:\n")
    map.ct.print()

C = rad.compute_circ(lat)
print(f"\nC [m] = {C:5.3f}")
ind = np.zeros(nv, int)
print("\nalpha:")
for k in range(1, no+1):
    ind[delta_] = k
    print("  {:1d} {:10.3e}".format(k, map.ct.get(ind)/C))

compute_disp(map, no, nv, desc)
