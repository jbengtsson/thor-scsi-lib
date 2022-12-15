"""Test the selection of matrix subpart ...
"""
import numpy as np
import thor_scsi.lib as tslib


tmp = np.zeros([7, 7], np.double)

idx = np.arange(6)
tmp[idx, idx] = .8 + idx/10
tmp[6,6] = 1

tmp[0, 1] = .12
tmp[1, 0] = .21
tmp[2, 3] = .23
tmp[3, 2] = .32
tmp[4, 5] = .45
tmp[5, 4] = .54
ss2 = tslib.mat_to_ss_vect_tps(tslib.Matrix(tmp))
jj = np.ones(6, np.int_)
jj[4:] = 0
#jj = list(jj)

print(tmp)
print(ss2)
r = tslib.select_subpart(ss2, jj)
print(r)