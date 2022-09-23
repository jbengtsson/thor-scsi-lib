"""Read lattice file and calculate radiation
"""
import logging
import xarray as xr

logging.basicConfig(level="DEBUG")
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.accelerator import instrument_with_radiators
from thor_scsi.utils.radiate import calculate_radiation
import os
import numpy as np

import thor_scsi.lib as tslib

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)
print(" ".join([elem.name for elem in acc]))
print("Length", np.sum([elem.getLength() for elem in acc]))

b2 = acc.find("b2", 0)

energy = 2.5e9


# cav.setVoltage(cav.getVoltage() * 1./2.)
# cav.setVoltage(0)
cav = acc.find("cav", 0)
print("acc cavity", repr(cav))
txt=\
    f"""Cavity info
frequency         {cav.getFrequency()/1e6} MHz",
voltage           {cav.getVoltage()/1e6} MV
harmonic number   {cav.getHarmonicNumber()}
    """
print(txt)

radiate = True
if radiate:
    calc_config = tslib.ConfigType()
    calc_config.radiation = True
    # is this used anywhere?
    calc_config.emittance = False
    calc_config.Cavity_on = True

    print(
        "calc_config",
        calc_config.radiation,
        calc_config.emittance,
        calc_config.Cavity_on,
    )

    r = calculate_radiation(
        acc, energy=2.5e9, calc_config=calc_config, install_radiators=True
    )

exit()

use_tpsa = True
if not use_tpsa:
    ps = ss_vect_double()
    ps.set_zero()
    ps[phase_space_ind.x_] = 1e-3
else:
    ps = ss_vect_tps()
    ps.set_identity()


# First step:
#
# use closed orbit
# 1. calculate fix point and Poincarè Map M with damped system (i.e. radiation on
#    and cavity on (without dispersion in a second case)
# 2. diagonalise M = A $\Gamma$ A$^{-1}$
# 3. eigenvalues:
#        - complex part: tunes,
#        - real part: damping times  (refer equation)
#    use eigen values of symplectic matrix to identify the planes
# 4. propagate A, thin kick will create diffusion coeffs (don't forget to zero
#    them before calculation starts (sum it up afterwards


print(ps)
acc.propagate(calc_config, ps, 0, 2000)
print(ps)


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
