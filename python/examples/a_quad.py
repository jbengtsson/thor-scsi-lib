"""Configuring and using a quadrupole
"""
from thor_scsi.pyflame import Config
import thor_scsi.lib as tslib
print(dir(tslib))
from thor_scsi.lib import Quadrupole, ConfigType
import gtpsa
import gtpsa.utils

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
quadK = tslib.QuadrupoleTpsa(C)

index_dict = dict(x=0, px=1, y=2, py=3, delta=4, ct=5,  K=6, dx=7, dy=8)
index_mapping = gtpsa.IndexMapping(index_dict, "quadrupole study")

mo = 4
desc = gtpsa.desc(len(index_dict), mo)

# Index mapping should not be necessary
ps = gtpsa.ss_vect_tpsa(desc, 4, len(index_dict), index_mapping)
ps.set_identity()

# Make the gradient a knob
# 1. get the value that was set in the initalisation
# 2. convert it to a base python object, variant not
#    required in python
# 3. stuff it in the place
k = quadK.get_main_multipole_strength()
#
k = k.to_object()
print("K", k, type(k))

# Make K a knob
print(gtpsa.ctpsa.__init__.__doc__)
b2 = gtpsa.ctpsa(desc, mo, mapping=index_mapping)
b2.set_variable(k + 0j, "K", 0e0j)
b2.print("K")

# this wrapping needs to be avoided ... but not now
# just requires a bit of rework of the c++ wrapper
tmp = gtpsa.CTpsaOrComplex(b2)
quadK.set_main_multipole_strength(tmp)

quadK.propagate(calc_config, ps)
print(ps)

for key in index_dict:
    getattr(ps, key).print(key, 1e-7)

coeffs = gtpsa.utils.tpsa_extract_coefficients_to_nrec(ps.x)
import pandas as pd
print("coeffs\n", pd.DataFrame(coeffs).set_index("index"))

ttmp = ps.x
print("ttmp", type(ttmp))
print(ttmp.get(dict(x=3)))

print(ps.jacobian())
