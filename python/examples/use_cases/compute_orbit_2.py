"""Use Case:
     Compute the orbit shift from reverse bends.
"""

import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="ERROR")
logger = logging.getLogger("thor_scsi")

import os
import sys
from dataclasses import dataclass
from typing import ClassVar

import math
import numpy as np

from scipy.interpolate import CubicSpline

import gtpsa

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, index_class as ind


ind = ind.index_class()


@dataclass
class gtpsa_prop:
    # GTPSA properties.
    # Number of variables - phase-space coordinates & 1 for parameter
    #dependence.
    nv: ClassVar[int] = 6 + 1
    # Max order.
    no: ClassVar[int] = 1
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]


def get_lat(lat_file_name, file_name, E_0):
    lat_prop = lp.lattice_properties_class(
        gtpsa_prop, lat_file_name, E_0, cod_eps)
    lat_prop.prt_lat(file_name)
    return lat_prop


def get_phi(el):
    if (type(el) == ts.Bending) or (type(el) == ts.Quadrupole):
        if el.get_curvature() is not None:
            L = el.get_length()
            phi = math.degrees(L*el.get_curvature())
        else:
            phi = 0e0
    else:
        phi = 0e0
    return phi


def compute_layout(lat_prop):
    s_buf = []
    X_buf = []
    Y_buf = []
    p_x_buf = []
    s = 0e0
    s_ref = s
    X = 0e0
    Y = 0e0
    p_x = 0e0
    s_buf.append(s)
    X_buf.append(X)
    Y_buf.append(Y)
    p_x_buf.append(p_x)
    for k, elem in enumerate(lat_prop._lattice):
        elem = lat_prop._lattice[k]
        L = elem.get_length()
        phi = math.radians(get_phi(elem))
        if phi == 0e0:
            X += L*np.cos(p_x)
            Y += L*np.sin(p_x)
        else:
            rho = L/phi
            X += rho*(np.sin(p_x+phi)-np.sin(p_x))
            Y += rho*(np.cos(p_x)-np.cos(p_x+phi))
        s += L
        p_x += phi
        # Cubic spline fit requires strictly monotonic series.
        if s > s_ref:
            s_buf.append(s)
            X_buf.append(X)
            Y_buf.append(Y)
            p_x_buf.append(p_x)
            s_ref = s
    s_buf = np.array(s_buf)
    X_buf = np.array(X_buf)
    Y_buf = np.array(Y_buf)
    p_x_buf = np.array(p_x_buf)
    X_cs = CubicSpline(s_buf, X_buf, bc_type="natural")
    Y_cs = CubicSpline(s_buf, Y_buf, bc_type="natural")
    return s_buf, X_cs, Y_cs, p_x_buf


def compute_orbit(s_ref, X_ref, Y_ref, p_x_ref, X, Y, k):
    Dx = X(s_ref[k]) - X_ref(s_ref[k])
    Dy = Y(s_ref[k]) - Y_ref(s_ref[k])
    dx = np.sqrt((Dx*np.sin(p_x_ref[k]))**2+(Dy*np.cos(p_x_ref[k]))**2)
    return Dx, Dy, dx


def prt_layout(file_name, s, X, Y, p_x):
    file = open(file_name, "w")

    print("# k      s            X          Y         p_x\n"
          "#                    [m]        [m]       [rad]", file=file)
    for k in range(len(s)):
        print(f"{k:4d}  {s[k]:9.5f}  {X(s[k]):9.5f}  {Y(s[k]):9.5f}"
              f"  {p_x[k]:9.5f}", file=file)


def prt_orbit(s_ref, X_ref, Y_ref, p_x_ref, X, Y):
    file_name = "orbit.txt"
    file = open(file_name, "w")

    print("# k      s            DX            DY            dx\n"
          "#                     [m]           [m]           [m]", file=file)
    for k in range(len(s_ref)):
        Dx, Dy, dx = compute_orbit(s_ref, X_ref, Y_ref, p_x_ref, X, Y, k)
        print(f"{k:4d}  {s_ref[k]:9.5f}  {Dx:12.5e}  {Dy:12.5e}  {dx:12.5e}",
              file=file)


# TPSA max order.
gtpsa_prop.no = 1

cod_eps = 1e-15
E_0     = 3.0e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")

file_name_ref = os.path.join(home_dir, "max_iv", "max_iv_baseline_2.lat")
file_name = os.path.join(home_dir, sys.argv[1]+".lat")

lat_ref = get_lat(file_name_ref, "lattice_ref_lat.txt", E_0)
s_ref, X_ref, Y_ref, p_x_ref = compute_layout(lat_ref)
prt_layout("lattice_ref_layout.txt", s_ref, X_ref, Y_ref, p_x_ref)

lat_prop = get_lat(file_name, "lattice_lat.txt", E_0)
s, X, Y, p_x = compute_layout(lat_prop)
prt_layout("lattice_layout.txt", s, X, Y, p_x)

prt_orbit(s_ref, X_ref, Y_ref, p_x_ref, X, Y)
