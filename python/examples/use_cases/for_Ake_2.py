"""Use Case:
     Compute, print, and plot the global lattice properties.
"""

import os
from typing import Tuple
from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from scipy import optimize as opt

import gtpsa

from thor_scsi.utils import lattice_properties as lp, index_class as ind
from thor_scsi.utils.output import vec2txt


ind = ind.index_class()


@dataclass
class gtpsa_prop:
    # GTPSA properties.
    # Number of phase-space coordinates.
    nv: ClassVar[int] = 7
    # Max order for Poincaré map.
    no: ClassVar[int] = 1
    # Number of parameters.
    nv_prm: ClassVar[int] = 0
    # Parameters max order.
    no_prm: ClassVar[int] = 0
    # Index.
    named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))
    # Descriptor
    desc : ClassVar[gtpsa.desc]


class new():
    def tpsa():
        return gtpsa.tpsa(gtpsa_prop.desc, gtpsa_prop.no)
    def ss_vect_tpsa():
        return gtpsa.ss_vect_tpsa(gtpsa_prop.desc, gtpsa_prop.no)


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


# Variables max order.
gtpsa_prop.no = 5
gtpsa_prop.desc = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)

cod_eps = 1e-15
E_0     = 3.0e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV", "max_iv")
lat_name = "ake_0"
file_name = os.path.join(home_dir, lat_name+".lat")

lat_prop = \
    lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

print("\nTotal bend angle [deg] = {:7.5f}".format(lat_prop.compute_phi_lat()))
print("Circumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))

# Compute the linear optics.
compute_optics(lat_prop)

lat_prop.prt_lat_param()
lat_prop.prt_M()
lat_prop.prt_rad()
# lat_prop.prt_M_rad()
# lat_prop.prt_Twiss(lat_name+"_Twiss.txt")

if not False:
    lat_prop.plt_Twiss(lat_name+"_Twiss.png", not False)
    lat_prop.plt_chrom(lat_name+"_chrom.png", not False)
