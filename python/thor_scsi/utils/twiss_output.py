import xarray as xr
import pandas as pd
import numpy as np


def twiss_ds_to_df(ds: xr.Dataset) -> pd.DataFrame:
    """Export twiss parameters from a dataset to a dataframe

    Simplifies further export to comma separated values (if that is
    still required).

    For storing the dataset have a look to
    :any:`xarray.DataSet.to_netcdf`
    """
    columns = [
        "name",
        "s",
        "alpha_x",
        "beta_x",
        "nu_x",
        "eta_x",
        "etap_x",
        "alpha_y",
        "beta_y",
        "nu_y",
        "eta_y",
        "etap_y",
    ]

    df = pd.DataFrame(columns=columns, index=ds.coords["index"], data=np.nan)
    df.loc[:, "name"] = ds.names
    df.loc[:, "s"] = ds.s
    for p in "x", "y":
        df.loc[:, f"alpha_{p}"] = ds.twiss.sel(par="alpha", plane=p)
        df.loc[:, f"beta_{p}"] = ds.twiss.sel(par="beta", plane=p)

        dnu = ds.twiss.sel(par="dnu", plane=p).values
        df.loc[:, f"nu_{p}"] = np.add.accumulate(dnu)

    disp = ds.dispersion
    df.loc[:, "eta_x"] = disp.sel(phase_coordinate="x")
    df.loc[:, "etap_x"] = disp.sel(phase_coordinate="px")
    df.loc[:, "eta_y"] = disp.sel(phase_coordinate="y")
    df.loc[:, "etap_y"] = disp.sel(phase_coordinate="py")

    return df


__all__ = ["twiss_ds_to_df"]
