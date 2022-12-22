"""Configuring and using a mpole
"""
from thor_scsi.pyflame import Config
from thor_scsi.lib import Mpole, ConfigType
from gtpsa import ss_vect_double

C = Config()
C.setAny("L", 2)
C.setAny("name", "mpole")
C.setAny("N", 4)
mpole = Mpole(C)

print(mpole)
print(repr(mpole))
mpole.get_field_interpolator().set_multipole(6, 4.2e-4)

# Will fail now: need to implement that the coefficients array can be resized
# quad.getMultipoles().setMultipole(10, 3.552e-4)

# The argument passed to the state space vector is required for
# compatability with gtpsa demanding a description object
ps = ss_vect_double(0.0)
ps.set_zero()
ps[0] = 1e-3
ps[2] = -1e-3
calc_config = ConfigType()
mpole.propagate(calc_config, ps)
print("One step\n", ps)
