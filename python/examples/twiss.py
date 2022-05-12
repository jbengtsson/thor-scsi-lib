"""Read lattice file and calculate radiation
"""
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv
from thor_scsi.utils.linear_optics import compute_twiss_along_lattice
import thor_scsi.lib as tslib

import os

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)
calc_config = tslib.ConfigType()

ds = compute_twiss_along_lattice(acc, calc_config)
df = twiss_ds_to_df(ds)
# print(df)

with open("twiss.tsf", "wt") as fp:
    fp.write(df_to_tsv(df))

df.to_json("twiss.json")
