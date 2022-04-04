"""Read lattice file and calculate radiation
"""
from thor_scsi.factory import accelerator_from_config
from thor_scsi.lib import (
    ConfigType,
    ss_vect_tps,
    ss_vect_double,
    RadiationDelegate,
    RadiationDelegateKick,
    phase_space_ind,
    ObservedState
)
import os

t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)

calc_config = ConfigType()
calc_config.radiation = False
calc_config.Cavity_on = False
calc_config.emittance = False

radiate = True
energy = 2.5e9


class RK(RadiationDelegateKick):
    def __init__(self):
        RadiationDelegateKick.__init__(self)

    def view(self, fk, ps, state, cnt):
        name = fk.name
        txt = f"fk.view '{name}'; state {state} "
        print(txt)
        return RadiationDelegateKick.view(self, fk, ps, state, cnt)


if radiate:
    calc_config.radiation = True
    calc_config.emittance = True

    # Should we add radiate delegates by default ?

    # Add radiators to elements that can radiate
    # Radiation delegate for marker
    # calc_config.emittance = True
    type_name = "Marker"
    ps_zero = ss_vect_double()
    rad_del = [RadiationDelegate() for elem in acc.elementsWithNameType(type_name)]
    for a_del, elem in zip(rad_del, acc.elementsWithNameType(type_name)):
        elem.setRadiationDelegate(a_del)
        # Just use that that the marker knows who is calling him
        a_del.view(elem, ps_zero, ObservedState.start, 0)

    # Not used any more better to clean up the name space
    del ps_zero

    # Radiation delegate for field kick
    type_name = "Bending"
    radiators = [elem for elem in acc.elementsWithNameType(type_name)]
    rad_del_kick = [
        # RK()
        RadiationDelegateKick()
        for elem in radiators
    ]
    for a_del, elem in zip(rad_del_kick, radiators):
        # Should be set from accelerator
        a_del.setEnergy(energy)
        elem.setRadiationDelegate(a_del)

for e in radiators:
    print(repr(e))
    break

use_tpsa = True
if not use_tpsa:
    ps = ss_vect_double()
    ps.set_zero()
    ps[phase_space_ind.x_] = 1e-3
else:
    ps = ss_vect_tps()
    ps.set_identity()

print(ps)
acc.propagate(calc_config, ps,  0, 2000)
print(ps)


if use_tpsa:
    # Inspect curly_H in
    for a_del in rad_del:
        name = a_del.getDelegatorName()
        idx = a_del.getDelegatorIndex()
        curly_H_x = a_del.getCurlydHx()
        txt = f"{name:10s} {idx:4d} curly_H_x {curly_H_x:5f}"
        print(txt)

    for a_del in rad_del_kick:
        name = a_del.getDelegatorName()
        idx = a_del.getDelegatorIndex()
        curly_H_x = a_del.getCurlydHx()
        dI = a_del.getSynchrotronIntegralsIncrements()
        D_rad = a_del.getDiffusionCoefficientsIncrements()
        txt = f"{name:10s} {idx:4d} curly_H_x {curly_H_x:5f}"
        txt += "    dI " + ",".join(["{: 10.6f}".format(v) for v in dI])
        txt += "   "
        txt += "    D_rad" + ",".join(["{: 10.6f}".format(v) for v in D_rad])
        txt += "   "
        print(txt)
