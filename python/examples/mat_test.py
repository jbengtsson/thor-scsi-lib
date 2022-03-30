from thor_scsi.lib import ss_vect_tps, ss_vect_tps_to_mat
import numpy as np

ps = ss_vect_tps()
ps.set_identity()
print(ps)

tmp = ss_vect_tps_to_mat(ps)
tmp = np.array(tmp)
np.set_printoptions(precision=4)
print(tmp)

for r in tmp:
    txt = '\t'.join(["{:10.6e}".format(x) for x in r])
    print(txt)
