import thor_scsi
from thor_scsi.pyflame import Config
from thor_scsi.lib import Accelerator, AcceleratorTpsa
from thor_scsi.lib import Bending, Cavity, Drift, Marker, Quadrupole
from thor_scsi.lib import BendingTpsa, DriftTpsa, QuadrupoleTpsa
import gtpsa

use_tpsa = True
if use_tpsa:
    Bending = BendingTpsa
    Drift =  DriftTpsa
    Quadrupole = QuadrupoleTpsa

# Transcript of TME lattice
Energy = 3e9
h_rf = 1320
c0 = 2.99892458e8
C = 10.1586800
Brho = Energy / c0

config = Config()
config.setAny("name", "Start")
start = Marker(config)

config = Config()
config.setAny("HarmonicNumber", h_rf)
config.setAny("name", "CAV")
config.setAny("Frequency", c0 / C * h_rf)
config.setAny("Voltage", 5e6 / 2)

cav = Cavity(config)

config = Config()
config.setAny("L", 0.9)
config.setAny("name", "D1")

d1 = Drift(config)

config.setAny("name", "D2")
config.setAny("L", 0.4)

d2 = Drift(config)

config.setAny("name", "D3")
config.setAny("L", 3.02934)
d3 = Drift(config)

elems = [start, d3, cav]
elems = [start, d3, cav]

if use_tpsa:
    acc = AcceleratorTpsa(elems)
else:
    acc = Accelerator(elems)

for e in acc:
    print(e)


gtpsa
