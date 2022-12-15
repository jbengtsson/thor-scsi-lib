'''Propagate turncated power series and observe it's results
'''
from collections import OrderedDict
from thor_scsi.factory import accelerator_from_config
from thor_scsi.lib import ConfigType
import gtpsa
from thor_scsi.observer import Observer
import numpy as np
import os.path

#  Find the file
t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)

calc_config = ConfigType()

# Register an observer with each element
# it seems that python observers can not be simply handed over
# to python. In this case the propagate does not work
observers = [Observer() for elem in acc]
for elem, ob in zip(acc, observers):
    elem.set_observer(ob)

# Start the calculation
desc = gtpsa.desc(6, 2)
ps = gtpsa.ss_vect_tpsa(desc, 1)
ps.set_identity()
print(ps)
acc.propagate(calc_config, ps, 0, 100)
print(ps)

# Results look a bit strange
np.set_printoptions(precision=2)
for cnt, ob in enumerate(observers[:3]):
    print(f'Observed {ob.index:3d}: {ob.name:20s}')
    print(ob.ps)
    for r in ob.jac:
        txt = '\t'.join(["{:10.6e}".format(x) for x in r])
        print(txt)

    print("\n\n")
