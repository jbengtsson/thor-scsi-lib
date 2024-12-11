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

import copy as _copy

import math
import numpy as np

import gtpsa

import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, nonlin_dyn as nld_cl, \
    closed_orbit as co, index_class as ind


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


def compute_optics(lat_prop):
    try:
        # Compute Twiss parameters along lattice.
        if not lat_prop.comp_per_sol():
            print("\ncomp_per_sol: unstable")
            raise ValueError

        # Compute radiation properties.
        stable, stable_rad = lat_prop.compute_radiation()
        print(stable, stable_rad)
        if not stable:
            print("\ncompute_radiation: unstable")
            raise ValueError
    except ValueError:
        assert False


def get_lat(home_dir, lat_name, E_0):
    file_name = os.path.join(home_dir, lat_name+".lat")

    lat_prop = lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)
    lat_prop.prt_lat(lat_name+"_lat.txt")

    print("\n{:s}".format(lat_name))
    print("Circumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
    print("Total bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))

    lat_prop.prt_lat(lat_name+"_lat.txt")

    # Computes element s location: lat_prop._Twiss.s[].
    compute_optics(lat_prop)
    lat_prop._model_state.radiation = False
    lat_prop._model_state.Cavity_on = False

    if False:
        lat_prop.prt_lat_param()
        lat_prop.prt_rad()
        lat_prop.prt_M()
        lat_prop.prt_M_rad()

    if False:
        lat_prop.plt_Twiss("compute_orbit_Twiss.png", not False)

    return lat_prop


def get_phi(el):
    L = el.get_length()
    if el.get_curvature() is not None:
        phi = math.degrees(L*el.get_curvature())
    else:
        phi = 0e0
    return phi


def prt_bend(lat_ref, lat_prop):
    file_name = "prt_bend.txt"
    file = open(file_name, "w")

    phi_tot_ref = 0e0
    phi_tot = 0e0
    dphi_tot = 0e0
    b_1xL_sum = 0e0
    print("#  k      s        name                  L              phi"
          "          dphi    b_1xL  b_1xL_sum\n"
          "#        [m]                            [m]            [deg]"
          "        [mrad]  [mrad]   [mrad]", file=file)
    for k in range(len(lat_ref._lattice)):
        el_ref = lat_ref._lattice[k]
        if (type(el_ref) == ts.Bending) or (type(el_ref) == ts.Quadrupole):
            fam_name_ref = el_ref.name
            s = lat_ref._Twiss.s[k]
            L_ref = el_ref.get_length()
            phi_ref = get_phi(el_ref)
            phi_tot_ref += phi_ref

            el = lat_prop._lattice[k]
            fam_name = el.name
            L = el.get_length()
            phi = get_phi(el)
            dphi = phi - phi_ref
            phi_tot += phi
            dphi_tot += dphi

            b_1xL = math.radians(dphi)
            b_1xL_sum += b_1xL
            el.get_multipoles().set_multipole(1, b_1xL/L)

            print("  {:3d}  {:6.3f}  {:8s} ({:8s})  {:5.3f} ({:5.3f})" \
                  "  {:6.3f} ({:6.3f})  {:6.3f}  {:6.3f}   {:6.3f}".
                  format(k, s, fam_name_ref, fam_name, L_ref, L,
                         phi_ref, phi, 1e3*dphi*np.pi/180e0, 1e3*b_1xL,
                         1e3*b_1xL_sum), file=file)
    print("\n  phi_tot = {:8.5f} ({:8.5f}) dphi = {:9.3e}".
          format(phi_tot_ref, phi_tot, dphi_tot), file=file)
    file.close()


def compute_dx(lat_prop):
    eps = 1e-10
    M_0 = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    M_1 = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    p_x = 0e0
    print("\ncompute_dx:")
    while True:
        M_0.set_zero()
        M_0.px += p_x
        M_1 = _copy.copy(M_0)
        lat_prop._lattice.propagate(lat_prop._model_state, M_1)
        dp_x = (M_1+M_0).cst().px/2e0
        p_x -= dp_x
        print("  p_x = {:10.3e}".format(p_x))
        if abs(dp_x) < eps:
            break
    return p_x


def prt_orbit(lat_prop):
    file_name = "prt_orbit.txt"
    file = open(file_name, "w")

    dp_x = compute_dx(lat_prop)

    M_0 = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    M_1 = gtpsa.ss_vect_tpsa(lat_prop._desc, 1)
    M_1.set_zero()
    M_1.px += dp_x
    print("# k               s    type     x       p_x      dp_x\n"
          "#                [m]           [mm]    [mrad]   [mrad]",
          file=file)
    for k in range(len(lat_prop._lattice)):
        M_0 = _copy.copy(M_1)
        lat_prop._lattice.propagate(lat_prop._model_state, M_1, k, 1)
        s = lat_ref._Twiss.s[k]
        el = lat_prop._lattice[k]
        print("{:3d}  {:8s}  {:6.3f}  {:4.1f}  {:7.3f}  {:7.3f}  {:7.3f}".
              format(
                  k, el.name, s, lat_prop._type_code[k],
                  1e3*M_1.cst().x, 1e3*M_1.cst().px,
                  1e3*(M_1-M_0).cst().px), file=file)


# TPSA max order.
gtpsa_prop.no = 2

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
delta_max = 3e-2
beta_inj  = np.array([3.0, 3.0])

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
home_dir_1 = os.path.join(home_dir, "max_iv")
home_dir_2 = os.path.join(home_dir, "max_4u")
# lat_ref = get_lat(home_dir_1, "max_iv_baseline", E_0)
# lat_prop = get_lat(home_dir_2, "max_4u_sp_jb_5", E_0)
lat_ref = get_lat(home_dir_1, "max_iv_baseline_3", E_0)
lat_prop = get_lat(home_dir_2, "m4U_240831_f01_08_06_01_lattice_7", E_0)

prt_bend(lat_ref, lat_prop)

b_n_list = [ "s1_f1", "s2_f1", "s3_f1", "s4_f1", "o1_f1", "o2_f1", "o3_f1"]
nld = nld_cl.nonlin_dyn_class(lat_prop, A_max, beta_inj, delta_max, b_n_list)
if True:
    nld.zero_mult(lat_prop, 2)
nld.zero_mult(lat_prop, 3)
nld.zero_mult(lat_prop, 4)

prt_orbit(lat_prop)
