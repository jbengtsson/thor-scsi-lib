"""Read lattice file and calculate radiation
"""
import logging
import copy
import xarray as xr

logging.basicConfig(level="DEBUG")
logger = logging.getLogger("thor_scsi")
print("\nlogger:\n", logger)
logger.setLevel("DEBUG")
# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.

from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.accelerator import instrument_with_radiators
from thor_scsi.utils.courant_snyder import compute_twiss_A
from thor_scsi.utils.phase_space_vector import ps_jac2ss_vect_tps

from thor_scsi.lib import (
    ConfigType,
    ss_vect_tps,
    ss_vect_double,
    RadiationDelegate,
    RadiationDelegateKick,
    ObservedState
)

from thor_scsi.lib import phase_space_index_internal as phase_space_ind

import os
import numpy as np

import thor_scsi.lib as tslib

import numpy as np
import scipy as sp

import thor_scsi.lib as tslib

# from thor_scsi.utils.linalg import match_eigenvalues_to_plane_orig
from thor_scsi.utils.closed_orbit import compute_closed_orbit
from thor_scsi.utils.output import vec2txt, mat2txt, chop_array
from thor_scsi.utils.linear_optics import compute_M_diag, calculate_nu_symp, \
    acos2


X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]

x_, px_, y_, py_, ct_, delta_ = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta
]


def extract_diffusion_coefficient():
    dD_rad = \
        np.array([rk.getDiffusionCoefficientsIncrements()
                  for rk in rad_del_kicks])
    print("\nextract_diffusion_coefficient:")
    print(dD_rad)
    D_rad = np.add.accumulate(dD_rad, axis=0)
    return D_rad


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
# t_file = os.path.join(t_dir, "b3_tst.lat")
t_file = os.path.join(t_dir, "b3_sf_40Grad_JB.lat")

acc = accelerator_from_config(t_file)
print(" ".join([elem.name for elem in acc]))
print("\nC = ", np.sum([elem.getLength() for elem in acc]))

# b2 = acc.find("b2", 0)

energy = 2.5e9
# Just to test diffusion
# energy = 2.5

# Planes x y z
# from thor_scsi.lib import spatial_index
# print(spatial_index)

cav = acc.find("cav", 0)
# cav.setVoltage(cav.getVoltage() * 1./2.)
# cav.setVoltage(0)
print("\nCavity", repr(cav))
txt = f"""\nCavity info:
  f [MHz] {1e-6*cav.getFrequency()}",
  V [MV]  {1e-6*cav.getVoltage()}
  h       {cav.getHarmonicNumber()}
"""
print(txt)

mbb = acc.find("mbb", 0)
print("{:s}: N = {:d}".
      format(mbb.name, mbb.getNumberOfIntegrationSteps()))

# Install radiators that radiation is calculated
rad_del_kicks = instrument_with_radiators(acc, energy=energy)

calc_config = tslib.ConfigType()

calc_config.radiation = True
calc_config.emittance = False
calc_config.Cavity_on = True

print("\ncalc_config: [radiation, emittance, Cavity_on] = ",
      calc_config.radiation, calc_config.emittance, calc_config.Cavity_on)

calc_config.Energy = energy

if calc_config.Cavity_on == True:
    dof = 3
else:
    dof = 2

ps = tslib.ss_vect_double()
ps.set_zero()
ps[x_]     =  1e-10
ps[y_]     = -1e-10
ps[delta_] =  1e-10
print("\nps_0 = ", ps, end="")
# acc.propagate(calc_config, ps, 8-1, 8)
acc.propagate(calc_config, ps)
print("ps_1 = ", ps)

#exit()

r = compute_closed_orbit(acc, calc_config, delta=0e0)
M = r.one_turn_map[:6, :6]
print("M:\n" + mat2txt(M))

#exit()

# r = calculate_radiation(
#     acc, energy=2.5e9, calc_config=calc_config, install_radiators=True
# )

A, A_inv = compute_M_diag(dof, M)

calc_config.emittance = True

A_cpy = ps_jac2ss_vect_tps(r.x0, A)
acc.propagate(calc_config, A_cpy)
extract_diffusion_coefficient()

exit()


use_tpsa = True
if not use_tpsa:
    ps = tslib.ss_vect_double()
    ps.set_zero()
    ps[tslib.phase_space_ind.x_] = 1e-3
else:
    ps = tslib.ss_vect_tps()
    ps.set_identity()


# First step:
#
# use closed orbit
# 1. calculate fix point and Poincarè Map M with damped system (i.e. radiation
#    on and cavity on (without dispersion in a second case)
# 2. diagonalise M = A $\Gamma$ A$^{-1}$
# 3. eigenvalues:
#        - complex part: tunes,
#        - real part: damping times  (refer equation)
#    use eigen values of symplectic matrix to identify the planes
# 4. propagate A, thin kick will create diffusion coeffs (don't forget to zero
#    them before calculation starts (sum it up afterwards




# print("Poincaré map calulation .... ")
# print("radiation OFF")
# print("start point")
# ps_start = copy.copy(ps)
# print(ps)
# print(ps.cst())



# calc_config.radiation = False
# acc.propagate(calc_config, ps, 0, 2000)
# print(ps)


print("\n\nradiation ON")
calc_config.radiation = True
calc_config.emittance = True
calc_config.Cavity_on = True
# ps_wr = copy.copy(ps_start)

r = compute_closed_orbit(acc, calc_config, delta=0e0)
M = r.one_turn_map[:6, :6]

# print("\nM:\n" + mat2txt(M))
print("\n\nFixed point:\n", vec2txt(r.x0))
tune_x, tune_y, tune_long = calculate_nu_symp(3, M)
print(f"\n{tune_x=:.16f} {tune_y=:.16f} {tune_long=:.16f}")

exit()


print(ps_wr)
acc.propagate(calc_config, ps_wr, 0, 2000)
print(ps_wr.cst())
print(ps_wr)
print("Effect of radiation")
print(ps_wr - ps)
print(ps_wr.cst() - ps.cst())


exit()
print(extract_diffusion_coefficient())
use_tpsa = False
if use_tpsa:
    # Inspect curly_H in
    for a_del in rad_del:
        name = a_del.getDelegatorName()
        idx = a_del.getDelegatorIndex()
        curly_H_x = a_del.getCurlydHx()
        txt = f"{name:10s} {idx:4d} curly_H_x {curly_H_x:5f}"
        print(txt)

    I = np.array([a_del.getSynchrotronIntegralsIncrements() for a_del in rad_del_kick])

    for a_del in rad_del_kick:
        name = a_del.getDelegatorName()
        idx = a_del.getDelegatorIndex()
        curly_H_x = a_del.getCurlydHx()
        dI = a_del.getSynchrotronIntegralsIncrements()
        D_rad = a_del.getDiffusionCoefficientsIncrements()

        txt = f"{name:10s} {idx:4d} curly_H_x {curly_H_x: 10.6e}"
        txt += "    dI " + ",".join(["{: 10.6e}".format(v) for v in dI])
        txt += "   "
        txt += "    D_rad" + ",".join(["{: 10.6e}".format(v) for v in D_rad])
        txt += "   "
        print(txt)
