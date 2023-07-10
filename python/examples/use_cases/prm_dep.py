import enum
import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

import os

from thor_scsi.pyflame import Config

import gtpsa
import thor_scsi.lib as tslib

from thor_scsi.utils import knobs


def prt_map(str, map):
    print(str)
    map.x.print()
    map.px.print()
    map.y.print()
    map.py.print()
    map.delta.print()
    map.ct.print()


def quad_dep_b_2_dx_dy():
    """
    """
    # Number of phase-space coordinates.
    nv = 6
    # Variables max order.
    no = 2
    # Number of parameters.
    nv_prm = 3
    # Parameters max order.
    no_prm = 1

    length = 0.1
    K      = 1e0
    N      = 4

    name = "QF"
    C = Config()
    C.setAny("L", length)
    C.setAny("name", name)
    C.setAny("K", K)
    C.setAny("N", N)
    quad = tslib.QuadrupoleTpsa(C)
    quad.set_number_of_integration_steps(N)

    model_state = tslib.ConfigType()

    named_index = gtpsa.IndexMapping(
        dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6, dx=7, dy=8)
    )

    desc = gtpsa.desc(nv, no, nv_prm, no_prm)

    knobs.make_magnet_knobbable(
        quad, po=2, offset=True, desc=desc, named_index=named_index
    )

    M = gtpsa.ss_vect_tpsa(desc, no, nv, named_index)
    M.set_identity()
    quad.propagate(model_state, M)
    print("\nM:\n", M)
    prt_map("\nM:", M)


def sext_dep_b_3():
    """
    """
    # Number of phase-space coordinates.
    nv = 6
    # Variables max order.
    no = 3
    # Number of parameters.
    nv_prm = 3
    # Parameters max order.
    no_prm = 1

    length = 0.1
    K      = 1e0
    N      = 4

    name = "SF"
    C = Config()
    C.setAny("L", length)
    C.setAny("name", name)
    C.setAny("K", K)
    C.setAny("N", N)
    sext = tslib.SextupoleTpsa(C)
    sext.set_number_of_integration_steps(N)

    model_state = tslib.ConfigType()

    named_index = gtpsa.IndexMapping(
        dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6, dx=7, dy=8)
    )

    desc = gtpsa.desc(nv, no, nv_prm, no_prm)

    knobs.make_magnet_knobbable(
        sext, po=2, offset=True, desc=desc, named_index=named_index
    )

    M = gtpsa.ss_vect_tpsa(desc, no, nv, named_index)
    M.set_identity()
    sext.propagate(model_state, M)
    print("\nM:\n", M)
    prt_map("\nM:", M)


quad_dep_b_2_dx_dy()

sext_dep_b_3()
