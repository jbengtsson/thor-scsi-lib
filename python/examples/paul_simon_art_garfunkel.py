import logging
logging.basicConfig(level="WARNING")
# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.

import os
import numpy as np
import matplotlib.pyplot as plt

import thor_scsi.lib as tslib

from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.linear_optics import (
    compute_map,
    compute_nu_symp,
    compute_M_diag,
    compute_twiss_along_lattice
)
from thor_scsi.utils.radiate import compute_radiation
from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import mat2txt, vec2txt
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv


X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]


def plt_twiss(ds):
    fig, (gr1, gr2) = plt.subplots(2)

    gr1.set_xlabel("s [m]")
    gr1.set_ylabel(r"$\beta_x, \beta_y$ [m]")
    gr1.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), label=r"$\beta_x$")
    gr1.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), label=r"$\beta_y$")
    gr1.legend()

    gr2.set_xlabel("s [m]")
    gr2.set_ylabel(r"$\eta_x [m]")
    # gr2.plot(ds.s, ds.twiss.sel(plane="x", par="eta"), label=r"$\eta_x$")
    gr2.legend()

    plt.show()


def prt_fam(acc, fam_name):
    print("\nprt_fam:")
    for q in acc.elementsWithName(fam_name):
        print("  {:4s} {:3d} {:6.3f} {:4.2f}".
              format(q.name, q.index, q.getMultipoles().getMultipole(2).real,
                     q.getLength()))


def set_db_2_fam(acc, fam_name, db_2):
    for q in acc.elementsWithName(fam_name):
        b_2 = q.getMultipoles().getMultipole(2).real
        index = q.index
        q.getMultipoles().setMultipole(2, b_2+db_2)


def set_db_2L_fam(acc, fam_name, db_2L):
    q = acc.find(fam_name, 0)
    L = q.getLength()
    set_db_2_fam(acc, q.name, db_2L/L)


def compute_dnu_db_2(acc, calc_config, fam_name, db_2L):
    dof = 2
    set_db_2L_fam(acc, fam_name, db_2L)
    M = map2numpy(compute_map(acc, calc_config))[:6, :6]
    nup = compute_nu_symp(dof, M)
    set_db_2L_fam(acc, fam_name, -db_2L)
    M = map2numpy(compute_map(acc, calc_config))[:6, :6]
    num = compute_nu_symp(dof, M)
    dnu_db_2 = (nup-num)/(2e0*db_2L)
    print("  {:5s} nu = {:7.5f} {:7.5f}".
          format(fam_name, dnu_db_2[X_], dnu_db_2[Y_]))
    return dnu_db_2


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "alsu-7ba-20180503c.lat")
# t_file = os.path.join(t_dir, "b3_sf_40Grad_JB.lat")

acc = accelerator_from_config(t_file)
calc_config = tslib.ConfigType()

cav = acc.find("cav", 0)
print("\nCavity", repr(cav))
print(f"""\nCavity info:
  f [MHz] {1e-6*cav.getFrequency()}
  V [MV]  {1e-6*cav.getVoltage()}
  h       {cav.getHarmonicNumber()}
  phi     {cav.getPhase()}   (initial)
""", end="")
# Set RF phase for negative alpha_c.
cav.setPhase(180e0)
print("  phi     {:} (set to)".format(cav.getPhase()))

dof = 2

E = 2.0e9

M = map2numpy(compute_map(acc, calc_config))[:6, :6]
print("\nM:\n", mat2txt(M))
A, A_inv, _ = compute_M_diag(dof, M)

ds = compute_twiss_along_lattice(acc, calc_config, A)
df = twiss_ds_to_df(ds)
# print("\nTwiss - ds:\n", ds)
# print("\nTwiss - df:\n", df)

with open("twiss.tsf", "wt") as fp:
    fp.write(df_to_tsv(df))
df.to_json("twiss.json")

if not True:
    plt_twiss(ds)

compute_radiation(acc, calc_config, E, 1e-15)


# Tweak the tune by using: [q1, q2, q3, q4, mb3h, qfh, mqfh].
dnu_db_2 = compute_dnu_db_2(acc, calc_config, "q1",   1e-4)
print(dnu_db_2)
np.append(dnu_db_2, compute_dnu_db_2(acc, calc_config, "q2",   1e-4))
print(dnu_db_2)
compute_dnu_db_2(acc, calc_config, "q3",   1e-4)
compute_dnu_db_2(acc, calc_config, "q4",   1e-4)
compute_dnu_db_2(acc, calc_config, "mb3h", 1e-4)
compute_dnu_db_2(acc, calc_config, "qfh",  1e-4)
compute_dnu_db_2(acc, calc_config, "mqfh", 1e-4)

# SVD

# Tweak the tune.
