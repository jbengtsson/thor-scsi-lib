import xarray as xr
import pandas as pd
import numpy as np
import functools

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

format_dict = {
    "index": "{:4d}",
    "name": "{:10s}",
    "s": None,
    "alpha_x": None,
    "beta_x": None,
    "nu_x": None,
    "eta_x": None,
    "etap_x": None,
    "alpha_y": None,
    "beta_y": None,
    "nu_y": None,
    "eta_y": None,
    "etap_y": None,
}


def f(txt, fmt=None):
    return fmt.format(txt)


formatters = {key: functools.partial(f, fmt=fmt) for key, fmt in format_dict.items()}
_default_format = "{:9.5f}"


def find_format(col_name: str, default_format: str = None) -> str:
    if default_format is None:
        def_fmt = _default_format
    else:
        def_fmt = default_format

    t_fmt = format_dict[col_name]
    if t_fmt is None:
        t_fmt = def_fmt

    return t_fmt


def create_formatters(columns):
    def f(txt, fmt=None):
        return fmt.format(txt)

    formatters = {
        key: functools.partial(f, fmt=find_format(key)) for key in columns
    }
    return formatters


def df_to_lists(df: pd.DataFrame,) -> list:
    """Returns list of lists containing formatted entries
    """

    formats = [format_dict["index"]] + [find_format(col) for col in df.columns]

    def format_line(index, t_series):
        fi = formats[0]
        txt_idx = fi.format(index)
        r = [txt_idx]
        r += [fmt.format(v) for fmt, v in zip(formats[1:], t_series)]
        return r

    lines = [format_line(idx, ser) for idx, ser in df.iterrows()]
    return lines


def df_header():
    """Create header for standard format dictionary
    """
    txt = "index "
    txt += format_dict["name"].format("name")
    txt += "    "
    d = format_dict.copy()
    del d["index"]
    del d["name"]
    tl = ["{:9s}".format(k) for k in d.keys()]
    txt += "  ".join(tl)
    return txt


def df_to_tsv(
    df: pd.DataFrame,
    # sep: str = "  ", line_terminator: str = "\n", **kwargs
) -> str:
    """
    see :func:`df_to_lists` for keyword documentation

    Args:
        sep: column separator
        line_terminator :
    Todo:
        check if using pandas style is not sufficient
    """

    formatters = create_formatters(df.columns)
    return df.to_string(formatters=formatters, justify="end")

    lines = df_to_lists(df, **kwargs)

    txt_lines = [sep.join(l) for l in lines]
    txt = line_terminator.join(txt_lines)
    return txt


def twiss_ds_to_df(ds: xr.Dataset) -> pd.DataFrame:
    """Export twiss parameters from a dataset to a dataframe

    Simplifies further export to comma separated values (if that is
    still required).

    For storing the dataset have a look to
    :any:`xarray.DataSet.to_netcdf`
    """

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

    num_vals = set(df.columns).difference(set("name"))
    num_vals = tuple(num_vals)

    return df


__all__ = ["twiss_ds_to_df"]
