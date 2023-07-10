# Find the lattice file
from pathlib import Path

import os

import tqdm

from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils import knobs
import thor_scsi.lib as tslib
from thor_scsi.utils.accelerator import instrument_with_standard_observers
import gtpsa.utils
import numpy as np
import pandas as pd
import xarray as xr
from typing import Sequence

named_index_d = dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6)
named_index = gtpsa.IndexMapping(named_index_d)
# default_desc = gtpsa.desc(nv=6+3, mo=3)
# maximum order of variables
mo = 3
# maximum order of parameters ... at least one
po = 2
default_desc = gtpsa.desc(nv=6, mo=mo, np=3, po=po)
# default_desc = gtpsa.desc(nv=6, mo=3, np=144*3, po=1)

def response_of_one_magnet(
    acc, magnet, *, desc: gtpsa.desc, mo=1, nv=6, calc_config = tslib.ConfigType(), **kwargs
):
    magnet = knobs.make_magnet_knobbable(magnet, po=po, named_index=named_index, desc=desc, offset=False, **kwargs)
    ps = gtpsa.ss_vect_tpsa(desc, mo, nv, index_mapping=named_index)
    ps.set_identity()
    acc.propagate(calc_config, ps)
    knobs.make_magnet_unknobbable(magnet)
    return ps


def combine_responses(observers: Sequence[tslib.StandardObserver]):
    def extract_x_y(ob):
        ps = ob.get_truncated_power_series_a()
        data = [
            (ps.x.get(d), ps.y.get(d))
            for d in (dict(K=1), dict(K=2))
        ]

        ob_name = ob.get_observed_name()
        da = xr.DataArray(
            data=[data],
            dims=["bpm", "dep", "plane"],
            coords=[
                [ob_name],
                ["K", "K2"],
                ["x", "y"],
            ],
            name="bpm_observation",
        )

        idx = ob.get_observed_index()
        index = xr.DataArray(name="index", data=[idx], dims=["bpm"], coords=[[ob_name]])
        ds = xr.Dataset(dict(index=index, effect=da))
        return ds

    data = [extract_x_y(ob) for ob in observers]
    ds = xr.concat(data, dim="bpm")
    return ds



def process_one_magnet(magnet, calc_config=tslib.ConfigType()):
    nps = response_of_one_magnet(acc_tpsa, magnet, desc=default_desc)
    acc_tpsa.propagate(calc_config, nps)
    response_one_magnet = combine_responses(observers)
    response_one_magnet["name"] = magnet.name
    # list comprehension as currently iloc not active?
    ps = [
        np.array(
            [
                nps.x.get(d),
                nps.px.get(d),
                nps.y.get(d),
                nps.py.get(d),
                nps.delta.get(d),
                nps.ct.get(d),
            ]
        )
        for d in [dict(K=1), dict(K=2)]
    ]
    phase_space = xr.DataArray(
        data=ps,
        name="phase_space",
        dims=["dep", "coor"],
        coords=[["K", "K2",], ["x", "px", "y", "py", "delta", "ct"]],
    )
    response_one_magnet["phase_space"] = phase_space
    return response_one_magnet


t_file = (
    Path(os.environ["HOME"])
    / "Devel"
    / "gitlab"
    / "dt4acc"
    / "lattices"
    / "b2_stduser_beamports_blm_tracy_corr.lat"
)

print(f"Reading lattice file {t_file}")
acc = accelerator_from_config(t_file)
acc_tpsa =  tslib.AcceleratorTpsa([knobs.convert_if_steerer(elem) for elem in acc])

# assuming that the fixed point corresponds to 0, 0,
# add observers before we propagate
# select bpms
bpm_elems = [elem for elem in acc if isinstance(elem, tslib.BPM)]
observers = instrument_with_standard_observers(bpm_elems, mapping=named_index)


all_responses = {
    elem.name: process_one_magnet(elem)
    for elem in tqdm.tqdm(acc_tpsa, total=len(acc_tpsa))
    if isinstance(elem, (tslib.HorizontalSteererTpsa, tslib.VerticalSteererTpsa))
}
tmp = [(key, item) for key, item in all_responses.items()]
db = xr.concat([t[1] for t in tmp], pd.Index([t[0] for t in tmp], name="steerer"))
db.to_netcdf("steerer_response_matrix.nc")
