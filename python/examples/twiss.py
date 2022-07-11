"""Read lattice file and calculate radiation
"""
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv
from thor_scsi.utils.linear_optics import compute_twiss_along_lattice
import thor_scsi.lib as tslib
import matplotlib.pyplot as plt

import os

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

t_file = os.path.join("lattices", "tme.lat")
# t_file = os.path.join("lattices", "tme_rb.lat")
acc = accelerator_from_config(t_file)
calc_config = tslib.ConfigType()

ds = compute_twiss_along_lattice(acc, calc_config)
df = twiss_ds_to_df(ds)
# print(df)

with open("twiss.tsf", "wt") as fp:
    fp.write(df_to_tsv(df))
df.to_json("twiss.json")

fig, ax = plt.subplots(1, 1, sharex=True)
ax.set_xlabel("s [m]")
ax.set_ylabel(r"$\beta_x, \beta_y$ [m]")
ax.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), label=r"$\beta_x$")
ax.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), label=r"$\beta_y$")
ax.legend()
plt.show()
