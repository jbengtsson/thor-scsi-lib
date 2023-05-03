"""Use Case:
     Modelling of BESSY II Transport Line: Booster to Storage Ring.
   Input:
     Lattice.
   Parameters:
     Twiss functions at entrance.
"""

import os

import numpy as np

import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as tslib

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import linear_optics as lo

# from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import prt2txt, mat2txt, vec2txt

from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv


tpsa_order = 1

# Descriptor for Truncated Power Series Algebra variables.
desc = gtpsa.desc(6, tpsa_order)


def plot_Twiss(ds, file_name, title):
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


t_dir = \
    os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "BESSY-II")
t_file = os.path.join(t_dir, "b2_transferline.lat")

# Read in & parse lattice file.
lat = accelerator_from_config(t_file)
# Set lattice state (Rf cavity on/off, etc.)
model_state = tslib.ConfigType()

n_dof = 4

# Compute the Poincaré map.
M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=tpsa_order)
print("\nM:\n", M)

# Twiss functions at the entrance.
eta   = np.array([0.0, 0.0, 0.0, 0.0])
alpha = np.array([1.0, 1.0])
beta  = np.array([1.0, 1.0])

A = lo.compute_A(eta, alpha, beta, desc)
print("\nA:\n", A)

data = lo.compute_Twiss_along_lattice(n_dof, lat, model_state, A=A, desc=desc)
plot_Twiss(data, "twiss.png", "Transport Line")
df = twiss_ds_to_df(data)
with open("twiss.tsf", "wt") as fp:
    fp.write(df_to_tsv(df))
df.to_json("twiss.json")

plt.show()
