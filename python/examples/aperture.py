"""Simple aperture check
"""

from thor_scsi.flame import GLPSParser
from thor_scsi.lib import Accelerator, RectangularAperture, ss_vect_double, ConfigType

mini_lat = """
Nquad = 12;
q4m2d1r: Quadrupole, L = 0.5, K = 1.4, N = Nquad, Method = 4;
mini_cell : LINE = (q4m2d1r);
"""

width = 100e-3
height = 50e-3
ap = RectangularAperture(width, height)

acc = Accelerator(GLPSParser().parse_byte(mini_lat, None))
elem = acc.find("q4m2d1r", 0)
elem.setAperture(ap)

print(f"Using aperture {ap} for element {elem}")

ps = ss_vect_double()
ps.set_zero()

calc_config = ConfigType()

calc_config.lossplane = 0
last_elem_idx = acc.propagate(calc_config, ps)
print(f"Propagate from center: last element index {last_elem_idx} loss plane {calc_config.lossplane}")

calc_config.lossplane = 0
ps.set_zero()
ps[0] = width + 1e-3
last_elem_idx = acc.propagate(calc_config, ps)
print(f"Propagate outside: last element index {last_elem_idx} loss plane {calc_config.lossplane}")
