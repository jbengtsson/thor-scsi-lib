'''Propagate turncated power series and observe it's results
'''
from collections import OrderedDict
from thor_scsi.factory import accelerator_from_config
from thor_scsi.lib import ConfigType, ss_vect_tps, ss_vect_double
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
    elem.setObserver(ob)

# Start the calculation
ps = ss_vect_tps()
ps.set_identity()
print(ps)
acc.propagate(calc_config, ps, 0, 3)

print(ps)

# Results look a bit strange
np.set_printoptions(precision=2)
for cnt, ob in enumerate(observers[:1]):
    print(f'Observed {ob.index:3d}: {ob.name:20s}')
    # print(ob.ps)

    for r in ob.res:

        txt = '\t'.join(["{:10.6e}".format(x) for x in r])
        print(txt)

    print("\n\n")

    # ps, mat = ss_vect_to_masked_matrix(ob.res)
    # print(ps)
    # print()
    # print(mat)
    # print()
    if cnt > 3:
        break
    if cnt >1234:
        break
