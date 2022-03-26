'''Propagate turncated power series and observe it's results
'''
from collections import OrderedDict
from thor_scsi.flame import GLPSParser
from thor_scsi.lib import (Accelerator, ConfigType, ss_vect_tps, ss_vect_double,
                           ss_vect_tps_to_lists)
from thor_scsi.utils import ss_vect_to_masked_matrix
from thor_scsi.observer import Observer

import numpy as np
import os.path

#  Find the file
t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

print(f'Reading lattice file {t_file}')
# I had issues reading the lattice file directly thus I preseent
# the string stream to the parser
with open(t_file) as fp:
    text = fp.read()

# Create parser, and create configuration
C = GLPSParser().parse_byte(text, t_dir)

# The machine creates it then based on the configuration
m = Accelerator(C)

calc_config = ConfigType()

# Register an observer with each element
# it seems that python observers can not be simply handed over
# to python. In this case the propagate does not work
observers = [Observer() for elem in m]
for elem, ob in zip(m, observers):
    elem.setObserver(ob)

# Start the calculation
ps = ss_vect_tps()
ps.set_identity()
m.propagate(calc_config, ps, 0, 2000)

# Results look a bit strange
np.set_printoptions(precision=2)
for cnt, ob in enumerate(observers):
    print(f'Observed {ob.index:3d}: {ob.name:20s}')
    print(ss_vect_to_masked_matrix(ob.ps).data)
    print()
    if cnt >2:
        break
