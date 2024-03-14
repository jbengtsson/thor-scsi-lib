"""Use Case:
     Parametric scans of the unit cell.
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
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
     radiate as rad

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


class lin_opt_class:
    # Private

    def __init__(self):
        self._prt_Twiss_file_name = "twiss.txt"
        self._plt_Twiss_file_name = "twiss.png"
        self._M                   = np.nan
        self._A                   = np.nan
        self._data                = np.nan
        self._nu                  = np.zeros(2, dtype=float)

    # Public.

    def get_type(self, elem):
        match type(elem):
            case ts.Marker:
                type_code = 0.0
            case ts.Drift:
                type_code = 0.0
            case ts.Bending:
                if elem.get_multipoles().get_multipole(2).real == 0e0:
                    type_code = np.sign(elem.get_curvature())*0.5
                else:
                    type_code = \
                        np.sign(elem.get_multipoles().get_multipole(2).real) \
                        *0.75
            case ts.Quadrupole:
                type_code = \
                    np.sign(elem.get_multipoles().get_multipole(2).real)*1.0
            case ts.Sextupole:
                type_code = \
                    np.sign(elem.get_multipoles().get_multipole(3).real)*1.25
            case ts.Octupole:
                type_code = \
                    np.sign(elem.get_multipoles().get_multipole(4).real)*1.5
            case _:
                type_code = np.nan
        return type_code

    def get_types(self, lat):
        n = len(lat)
        types = np.zeros(n)
        for k in range(n):
            types[k] = self.get_type(lat[k])
        return types

    def plt_Twiss(self, types, plot):
        fig, (gr_1, gr_2) = plt.subplots(2)

        fig.suptitle("Linear Optics")

        gr_1.set_title("Beta Functions")
        gr_1.set_xlabel("s [m]")
        gr_1.set_ylabel(r"$\beta_{x,y}$ [m]")
        gr_1.plot(
            self._data.s, self._data.twiss.sel(plane="x", par="beta"), "b",
            label=r"$\beta_x$")
        gr_1.plot(
            self._data.s, self._data.twiss.sel(plane="y", par="beta"), "r",
            label=r"$\beta_y$")
        gr_1.legend()

        gr_1_r = gr_1.twinx()
        gr_1_r.set_ylim([-2.0, 20.0])
        gr_1_r.set_yticks([])
        gr_1_r.step(self._data.s, types, "k")

        gr_2.set_title("Dispersion")
        gr_2.set_xlabel("s [m]")
        gr_2.set_ylabel(r"$\eta_x$ [m]")
        gr_2.plot(
            self._data.s, self._data.dispersion.sel(phase_coordinate="x"), "b")

        gr_2_r = gr_2.twinx()
        gr_2_r.set_ylim([-2.0, 20.0])
        gr_2_r.set_yticks([])
        gr_2_r.step(self._data.s, types, "k")

        fig.tight_layout()

        plt.savefig(self._plt_Twiss_file_name)
        print("\nPlot saved as:", self._plt_Twiss_file_name)

        if plot:
            plt.show()

    def prt_Twiss_param(self, str, Twiss):
        """
        """
        eta, alpha, beta = Twiss
        self._nu = lo.compute_nu_symp(2, self._M.jacobian())
        print(str, end="")
        print(f"  eta   = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
        print(f"  alpha = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
        print(f"  beta  = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")
        print(f"  nu    = [{self._nu[X_]:7.5f}, {self._nu[Y_]:7.5f}]")

    def prt_Twiss(self, lat, types):
        """
        Print Twiss parameters along the lattice.
        """
        file = open(self._prt_Twiss_file_name, 'w')
        s = 0e0
        nu = np.zeros(2, dtype=float)
        print("\n    Name         s    type  alpha_x   beta_x  nu_x    eta_x"
              "   eta'_x    alpha_y   beta_y  nu_y    eta_y   eta'_y",
              file=file)
        print("                [m]                    [m]             [m]"
              "     [m]                 [m]             [m]", file=file)
        for k in range(len(self._data.index)):
            s += lat[k].get_length()
            nu[X_] += self._data.twiss.sel(plane="x", par="dnu").values[k]
            nu[Y_] += self._data.twiss.sel(plane="y", par="dnu").values[k]
            print("{:3d} {:8s} {:7.3f} {:4.1f} {:9.5f} {:8.5f} {:7.5f} {:7.5f}"
                  " {:8.5f} {:9.5f} {:8.5f} {:7.5f} {:7.5f} {:8.5f}".
                  format(k, lat[k].name, s, self.get_type(lat[k]),
                         self._data.twiss.sel(plane="x", par="alpha").values[k],
                         self._data.twiss.sel(plane="x", par="beta").values[k],
                         nu[X_],
                         self._data.dispersion.sel(
                             phase_coordinate="x").values[k],
                         self._data.dispersion.sel(
                             phase_coordinate="px").values[k],
                         self._data.twiss.sel(plane="y", par="alpha").values[k],
                         self._data.twiss.sel(plane="y", par="beta").values[k],
                         nu[Y_],
                         self._data.dispersion.sel(
                             phase_coordinate="y").values[k],
                         self._data.dispersion.sel(
                             phase_coordinate="py").values[k]),
                  file=file)

    def comp_per_sol(self, lat, model_state, named_index, desc, prt):
        """
        Compute the periodic solution for a super period.
        Degrees of freedom - RF cavity is off; i.e., coasting beam.
        """
        n_dof = 2
        model_state.radiation = False
        model_state.Cavity_on = False

        stable, self._M, self._A = \
            lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc)
        if stable:
            if prt:
                print("\nM:\n"+mat2txt(self._M.jacobian()))
            res= cs.compute_Twiss_A(self._A)
            Twiss = res[:3]
            if prt:
                 self.prt_Twiss_param("\nTwiss:\n", Twiss)
            A_map = gtpsa.ss_vect_tpsa(desc, no)
            A_map.set_jacobian(self._A)
            self._data = \
                lo.compute_Twiss_along_lattice(
                    n_dof, lat, model_state, A=A_map, desc=desc,
                    mapping=named_index)
        else:
            self._data = np.nan
            print("\ncomp_per_sol: unstable")


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 2
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

E0 = 2.5e9

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(nv, no, nv_prm, no_prm)

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "BESSY-III",
    "ipac_2024")
t_file = os.path.join(home_dir, "6BA_0401_UC_with_RB_SVBR_tracy.lat")

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = ts.ConfigType()

# D.O.F. (Degrees-Of-Freedom) - coasting beam.
n_dof = 2
model_state.radiation = False
model_state.Cavity_on = False

lin_opt = lin_opt_class()

# Compute Twiss parameters along lattice.
lin_opt.comp_per_sol(lat, model_state, named_index, desc, True)

types = lin_opt.get_types(lat)

if not False:
   lin_opt.prt_Twiss(lat, types)

if not False:
    lin_opt.plt_Twiss(types, False)

