"""Configuring and using a quadrupole
"""
from thor_scsi.pyflame import Config
import thor_scsi.lib as tslib
print(dir(tslib))
from thor_scsi.lib import Quadrupole, ConfigType
import gtpsa

C = Config()
C.setAny("L", 2)
C.setAny("name", "qd")
C.setAny("K", 1.2)
C.setAny("N", 1)
quad = Quadrupole(C)

print(quad)
print(repr(quad))
quad.get_multipoles().set_multipole(6, 4.2e-4)

# Will fail now: need to implement that the coefficients array can be resized
# quad.getMultipoles().setMultipole(10, 3.552e-4)
print(repr(quad))
print(dir(quad))

desc = gtpsa.desc(6, 2)
ps = gtpsa.ss_vect_tpsa(desc, 1)
ps.set_identity()
calc_config = ConfigType()
quad.propagate(calc_config, ps)
print("One step\n", ps)


# Quadrupole with TPSA


C = Config()
C.setAny("L", 2)
C.setAny("name", "qd")
C.setAny("K", 1.2)
C.setAny("N", 1)
quad = tslib.QuadrupoleTpsa(C)

mo = 4
desc = gtpsa.desc(8, mo)

desc = gtpsa.desc(8, mo)
ps = gtpsa.ss_vect_tpsa(desc, 4)
ps.set_identity()
z = gtpsa.ctpsa(desc, mo)
muls = quad.get_multipoles()
print(muls)
K = quad.get_multipoles().get_multipole(1).real()
z.set(0, K)
