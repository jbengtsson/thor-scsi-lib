from thor_scsi.lib import (ss_vect_tps, ss_vect_tps_to_mat, Matrix,
                           mat_to_ss_vect_tps)
import numpy as np
7,
ps = ss_vect_tps()
ps.set_identity()
print(ps)

rtmp = ss_vect_tps_to_mat(ps)
print(rtmp)
tmp = np.array(rtmp)
np.set_printoptions(precision=4)
print(tmp)
for r in tmp:
    txt = '\t'.join(["{:10.6e}".format(x) for x in r])
    print(txt)

x = np.fromfunction(lambda x, y: x * 100 + y, (7, 7))
x = x.T
print(x)
m2 = Matrix(x)
print(m2)
ss2 = mat_to_ss_vect_tps(m2)
print(ss2)
