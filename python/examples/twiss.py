"""Read lattice file and calculate radiation
"""
import logging
#logging.basicConfig(level=logging.DEBUG)
import matplotlib.pyplot as plt
import os.path

from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv
from thor_scsi.utils.linear_optics import compute_Twiss_along_lattice
import thor_scsi.lib as tslib
from pathlib import Path
import gtpsa

desc = gtpsa.desc(6, 1)

t_file = Path(os.path.dirname(__file__)) /  "lattices" / "tme_rb.lat"
print("twiss from file", t_file)
acc = accelerator_from_config(t_file)
calc_config = tslib.ConfigType()
acc.set_log_level(tslib.accelerator_log_level.warning)

if True:
    n_dof = 2
    calc_config.radiation = False
    calc_config.Cavity_on = False
else:
    n_dof = 3
    calc_config.radiation = True
    calc_config.Cavity_on = True

ds = compute_Twiss_along_lattice(n_dof, acc, calc_config, desc=desc, mapping=gtpsa.IndexMapping())
df = twiss_ds_to_df(ds)
# print(df)
# exit()

with open("twiss.tsf", "wt") as fp:
    fp.write(df_to_tsv(df))
df.to_json("twiss.json")

logging.basicConfig(level=logging.WARNING)
fig, ax = plt.subplots(1, 1, sharex=True)
ax.set_xlabel("s [m]")
ax.set_ylabel(r"$\beta_x, \beta_y$ [m]")
ax.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), label=r"$\beta_x$")
ax.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), label=r"$\beta_y$")
ax.legend()
plt.show()
