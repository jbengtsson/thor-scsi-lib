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
energy = 1.7e9


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
    # calc_config.emittance = True

    # Should we add radiate delegates by default ?

    # Add radiators to elements that can radiate
    # Radiation delegate for marker
    # calc_config.emittance = True
    # type_name = "Marker"
    # rad_del = [RadiationDelegate() for elem in acc.elementsWithNameType(type_name)]
    # for a_del, elem in zip(rad_del, acc.elementsWithNameType(type_name)):
    #    elem.setRadiationDelegate(a_del)

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

# ps = ss_vect_tps()
# ps.set_identity()
ps = ss_vect_double()
ps.set_zero()
ps[phase_space_ind.x_] = 1e-3
print(ps)
acc.propagate(calc_config, ps,  0, 19)
print(ps)
acc.propagate(calc_config, ps, 19, 20)
print(ps)
for a_del in rad_del_kick:
    print(repr(a_del))
