"""Use Case:
     Compute, print, and display the global linear optics properties for a
     lattice.
"""


import os
import sys
from dataclasses import dataclass
from typing import ClassVar

import gtpsa

import numpy as np

from thor_scsi.utils import lattice_properties as lp, nonlin_dyn as nld_cl, \
    index_class as ind


ind = ind.index_class()


@dataclass
class gtpsa_prop:
    # GTPSA properties.
    # Number of variables - phase-space coordinates & 1 for parameter
    #dependence.
    nv: ClassVar[int] = 6 + 1
    # Max order.
    no: ClassVar[int] = 1
    # Truncation order.
    to: ClassVar[int] = 1
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]


def compute_optics(gtpsa_prop, lat_prop):
    try:
        # Compute Twiss parameters along lattice.
        lat_prop._desc.truncate(gtpsa_prop.no)
        if not lat_prop.comp_per_sol(gtpsa_prop):
            print("\ncomp_per_sol: unstable")
            raise ValueError

        # Adjust RF phase for sign of alpha_c.
        print("\ncompute_optics:")
        cav_loc = lat_prop._lattice.find("cav", 0)
        if lat_prop._alpha_c[1] >= 0e0:
            cav_loc.set_phase(0.0)
        else:
            print("  phi_rf  = 180 deg") 
            cav_loc.set_phase(180.0)

        # Compute radiation properties.
        stable, stable_rad = lat_prop.compute_radiation()
        print("\n", stable, stable_rad)
        if not stable:
            print("\ncompute_radiation: unstable")
            raise ValueError
    except ValueError:
        assert False


class new():
    def tpsa():
        return gtpsa.tpsa(gtpsa_prop.desc, gtpsa_prop.no)
    def ctpsa():
        return gtpsa.ctpsa(gtpsa_prop.desc, gtpsa_prop.no)
    def ss_vect_tpsa():
        return gtpsa.ss_vect_tpsa(gtpsa_prop.desc, gtpsa_prop.no)


def compute_dispersion(M):
    n = 4
    I = new.ss_vect_tpsa()
    ImM = new.ss_vect_tpsa()
    # Workaround for compose.
    D = new.ss_vect_tpsa()
    eta = new.ss_vect_tpsa()
    ImM_inv = new.ss_vect_tpsa()

    I.set_identity()
    ind_delta = [0, 0, 0, 0, 1]
    D.x.set([1, 0, 0, 0, 1], 0e0, M.x.get(ind_delta))
    D.px.set([1, 0, 0, 0, 1], 0e0, M.px.get(ind_delta))
    ImM_inv.pinv(I-M, [1, 1, 1, 1, 0, 0, 0])
    eta.compose(ImM_inv, D)
    eta.x.print("\neta.x:\n")
    eta.px.print("\neta.px:\n")
    eta_1, _, _, _ = lat_prop.get_Twiss(0)
    print("\neta^(1): [{:10.3e}, {:10.3e}]".format(eta_1[0], eta_1[1]))
    return


# TPSA max order.
gtpsa_prop.no = 3
gtpsa_prop.desc = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)

cod_eps = 1e-15
E_0     = 3.0e9

A_max     = np.array([6e-3, 3e-3])
beta_inj  = np.array([3.0, 3.0])
delta_max = 3e-2

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

lat_prop.prt_lat("lat_prop_lat.txt")

b_3_list = ["s1", "s2", "s3", "s4"]
nld = nld_cl.nonlin_dyn_class(
    gtpsa_prop, lat_prop, A_max, beta_inj, delta_max, b_3_list)
if False:
    nld.zero_mult(lat_prop, 3)
    nld.zero_mult(lat_prop, 4)

compute_optics(gtpsa_prop, lat_prop)
lat_prop._M.x.print("\nM.x:\n")

compute_dispersion(lat_prop._M)
exit

lat_prop.prt_lat_param(gtpsa_prop)
lat_prop.prt_rad()
lat_prop.prt_M()
lat_prop.prt_M_rad()

if not False:
    lat_prop.prt_Twiss("lat_prop_Twiss.txt")
    lat_prop.plt_Twiss("lat_prop_Twiss.png", not False)
    lat_prop.plt_chrom("lat_prop_chrom.png", not False)

if False:
    nld.compute_map(lat_prop, no)
    nld.compute_nl(lat_prop)
    nld.prt_nl(lat_prop)
