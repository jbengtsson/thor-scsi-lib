"""Use Case:
     Compute momentum compaction to arbitrary order.
"""


import os
import sys
from dataclasses import dataclass
from typing import ClassVar

import numpy as np
import matplotlib.pyplot as plt

import gtpsa
import thor_scsi.lib as ts

from thor_scsi.utils import lattice_properties as lp, closed_orbit as co, \
    index_class as ind

from thor_scsi.utils.output import mat2txt


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


class new():
    def tpsa():
        return gtpsa.tpsa(gtpsa_prop.desc, gtpsa_prop.no)
    def ss_vect_double():
        return gtpsa.ss_vect_double(0e0)
    def ss_vect_tpsa():
        return gtpsa.ss_vect_tpsa(gtpsa_prop.desc, gtpsa_prop.no)


class beam_dyn_class:
    # Private.

    def __init__(self, no):
        gtpsa_prop.no = no
        gtpsa_prop.desc   = gtpsa.desc(gtpsa_prop.nv, gtpsa_prop.no)
        self._M           = np.nan
        self._alpha_c     = np.nan
        self._delta       = np.nan
        self._ct          = np.nan
        self._circ        = lat_prop.compute_circ()
        self._phi         = lat_prop.compute_phi_lat()

        print(f"\nTotal bend angle [deg] = {self._phi:7.5f}")
        print(f"Circumference [m]      = {self._circ:7.5f}")

    # Public.

    def compute_map(self, delta):
        # Compute Poincaré Taylor M to arbitrary order.
        self._M = new.ss_vect_tpsa()
        # Compute closed orbit.
        cod = co.compute_closed_orbit(
            lat_prop._lattice, lat_prop._model_state, delta=delta, eps=1e-10,
            desc=gtpsa_prop.desc)
        cod.x0.iloc[ind.ct] = 0e0
        cod.x0.iloc[6] = 0e0
        print("\ncompute_map\n  cod = "+
              "[{:9.3e}, {:9.3e}, {:9.3e}, {:9.3e}, {:9.3e}, {:9.3e}]".format
              (cod.x0.x, cod.x0.px, cod.x0.y, cod.x0.py, cod.x0.delta,
               cod.x0.ct))
        # Initialise map to identity.
        self._M.set_identity()
        # Add closed orbit.
        self._M += cod.x0
        # Track one turn.
        lat_prop._lattice.propagate(lat_prop._model_state, self._M)
        self._M.ct.set([0, 0, 0, 0, 0, 0, 0], 0e0, 0e0)

    def compute_ct_tpsa(self):
        index = np.zeros(gtpsa_prop.nv, int)

        print("\nalpha_c:")
        alpha_c_buf = []
        # Set constant term to zero.
        alpha_c_buf.append(0e0)
        for k in range(1, gtpsa_prop.no+1):
            index[ind.delta] = k
            alpha_c = self._M.ct.get(index)/self._circ
            alpha_c_buf.append(alpha_c)
            print("  {:1d} {:10.3e}".format(k, alpha_c))
        self._alpha_c = np.array(alpha_c_buf)

        n_alpha_c = len(self._alpha_c)
        if n_alpha_c >= 3:
            print("\nFixed points to O(3) [%]: {:5.2f}, {:5.2f}".
                  format(0e0, -1e2*self._alpha_c[1]/self._alpha_c[2]))

        if n_alpha_c >= 4:
            po2 = self._alpha_c[2]/(2e0*self._alpha_c[3])
            q = self._alpha_c[1]/self._alpha_c[3]
            if po2**2-q > 0e0:
                pm = np.sqrt(po2**2-q)
                print("Fixed points to O(4) [%]: {:5.2f}, {:5.2f}, {:5.2f}".
                      format(0e0, -1e2*(po2+pm), -1e2*(po2-pm)))
            else:
                print("Fixed points to O(4) [%]: complex solution")

    def compute_ct_num(self, n_step, delta_max):
        # Compose of a Taylor map with floating point phase space vector not
        # (yet) supported.
        ps_fp   = new.ss_vect_double()
        ps_tpsa = new.ss_vect_tpsa()

        ddelta = delta_max/n_step
        delta = -delta_max
        print("\n    k     delta      ct TPSA      ct FP       diff")
        print("           [%]         [m]         [m]")
        for k in range(-n_step, n_step+1):
            ps_tpsa.set_zero()
            ps_tpsa.delta = delta
            ps_tpsa.compose(self._M, ps_tpsa)
            ct_tpsa = ps_tpsa.cst().ct

            ps_fp.set_zero()
            ps_fp.delta = delta
            lat_prop._lattice.propagate(lat_prop._model_state, ps_fp)
            ct_fp = ps_fp.ct

            diff = ct_tpsa - ct_fp
            print(f"  {k:3d}  {delta:10.3e}  {ct_tpsa:10.3e}  {ct_fp:10.3e}"
                  f"  {diff:10.3e}")

            delta += ddelta

    def graph_ct(self, n_step, delta_max):
        # Compose of a Taylor map with floating point phase space vector not
        # (yet) supported.
        ps  = new.ss_vect_tpsa()
        ps1 = new.ss_vect_tpsa()

        delta_buf= []
        ct_buf   = []
        ddelta   = delta_max/n_step
        delta    = -delta_max
        ps.set_zero()
        print()
        for k in range(-n_step, n_step+1):
            ps.delta = delta
            ps1.compose(self._M, ps)
            delta_buf.append(delta)
            ct_buf.append(ps1.cst().ct)
            print(" {:10.3e} {:10.3e}".
                  format(delta_buf[k+n_step], ct_buf[k+n_step]))
            delta += ddelta
        delta_buf = np.array(delta_buf)
        ct_buf = np.array(ct_buf)
        self. plt_ct_vs_delta(
            "ct_vs_delta", delta_buf, ct_buf, True)

    def plt_ct_vs_delta(self, file_name, delta, ct, plot):
        fig, gr = plt.subplots()

        fig.suptitle("Momentum Compaction")

        gr.set_title(r"$m_{56}(\delta)$")
        gr.set_xlabel(r"$\delta$ [%]")
        gr.set_ylabel("")
        gr.plot(1e2*delta, ct, "b")
        gr.legend("")

        fig.tight_layout()

        plt.savefig(file_name)
        print("\nplt_m_56_vs_delta - plot saved as:", file_name)

        if plot:
            plt.show()


cod_eps = 1e-15
E_0     = 3.0e9

# Read in lattice file.
home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
lat_name = sys.argv[1]
file_name = os.path.join(home_dir, lat_name+".lat")

# Process lattice file & instantiate lattice object.
lat_prop = \
    lp.lattice_properties_class(gtpsa_prop, file_name, E_0, cod_eps)

# Instantiate beam dynamics class & set TPSA max order - after processing
# lattice file.
bd = beam_dyn_class(5)

bd.compute_map(0e0)
# Print linear part.
print("\nM:", bd._M)

if False:
    bd.graph_ct(10, 5e-2)

# Print m_56, m_566, m_5666, etc.:
print("\nm_56, m_566, m_5666,...:")
print("  no    m_56...")
index = np.zeros(gtpsa_prop.nv, int)
for k in range(1, gtpsa_prop.no+1):
    index[ind.delta] = k
    print("  {:1d} {:10.3e}".format(k, bd._M.ct.get(index)))

# Print out momentum compaction & fixed points in the longitudinal (aka
# the "ignored") plane.
bd.compute_ct_tpsa()

bd.compute_ct_num(10, 5e-2)
