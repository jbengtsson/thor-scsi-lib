"""Read lattice file and calculate radiation
"""
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils import linear_optics
from thor_scsi.utils.extract_info import accelerator_info

from thor_scsi.lib import (
    ConfigType,
    ss_vect_tps,
    ss_vect_double,
    RadiationDelegate,
    RadiationDelegateKick,
    phase_space_ind,
    ObservedState
)

import xarray as xr
import os

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)

calc_config = ConfigType()
twiss = linear_optics.compute_twiss_parameters(acc, calc_config)
twiss.name = "twiss_parameters"
md = accelerator_info(acc)
md.attrs = dict(calc_config=calc_config)
twiss_with_md = xr.merge([twiss, md])

# print(res)

# combine
