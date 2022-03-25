"""Configuring and using a quadrupole
"""
from thor_scsi.flame import Config
from thor_scsi.lib import Quadrupole, ConfigType

C = Config()
C.setAny("L", 2)
C.setAny("name", "qd")
C.setAny("K", 1.2)
C.setAny("N", 1)
quad = Quadrupole(C)

print(quad)
print(repr(quad))
quad.getMultipoles().setMultipole(6, 4.2e-4)

# Will fail now: need to implement that the coefficients array can be resized
# quad.getMultipoles().setMultipole(10, 3.552e-4)
print(repr(quad))
print(dir(quad))
import thor_scsi.lib
print(dir(thor_scsi.lib))
# print(dir
from thor_scsi.lib import ss_vect_tps

ps = ss_vect_tps()
ps.set_identity()
calc_config = ConfigType()
quad.propagate(calc_config, ps)
print("One step\n", ps)
