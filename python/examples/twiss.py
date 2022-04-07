"""Read lattice file and calculate radiation
"""
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils import linear_optics
from thor_scsi.utils.extract_info import accelerator_info

from thor_scsi.lib import (
    ConfigType,
    ss_vect_tps,
    ss_vect_double,
    RadiationDelegate,
    RadiationDelegateKick,
    phase_space_ind,
    ObservedState,
    ss_vect_tps_to_mat
)

import xarray as xr
import os
import numpy as np


def prt_np_cmplx_vec(str, vec):
    print(str, end="")
    for k in range(len(vec)):
        print("{0.real:12.3e} {1} {0.imag:9.3e}i"
              .format(vec[k].real, "+-"[int(vec[k].imag < 0)], abs(vec[k].imag)),
              end="")
    print()


def prt_np_vec(str, vec):
    print(str, end="")
    for k in range(len(vec)):
        print("{:15.6e}".format(vec[k]), end="")
    print()


def prt_np_cmplx_mat(str, mat):
    print(str, end="")
    for k in range(len(mat)):
        prt_np_cmplx_vec("", mat[k])


def prt_np_mat(str, mat):
    print(str, end="")
    for k in range(len(mat)):
        prt_np_vec("", mat[k])


def chop_np_cmplx_vec(vec, eps):
    for k in range(len(vec)):
        if (np.abs(vec[k].real) < eps):
            vec[k] = complex(0e0, vec[k].imag)
        if (np.abs(vec[k].imag) < eps):
            vec[k] = complex(vec[k].real, 0e0)


def chop_np_cmplx_mat(mat, eps):
    for k in range(len(mat)):
        chop_np_cmplx_vec(mat[k], eps)


def compute_nu(mat):
    tr = mat.trace()
    if tr >= 2e0:
        print("Unstable Poincaré map: {:10.3e}".format(tr))
        exit(1)
    else:
        nu = np.arccos(tr/2e0)/(2e0*np.pi)
        if mat[0][1] < 0e0:
            nu = 1e0 - nu
        print("nu = {:7.5f}".format(nu))


def compute_nus(dof, mat):
    for k in range(dof):
        compute_nu(mat[2*k:2*k+2, 2*k:2*k+2])


def sign(x):
    if x > 0e0:
        return 1
    elif x < 0e0:
        return -1
    else:
        return 0


def compute_J(n):
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            J[i][j] = 0e0;

    for j in range(n):
        if j & 1:
            J[j-1][j] =  1e0
            J[j][j-1] = -1e0
    prt_np_mat("\nJ:\n", J)
    return J


def compute_A_inv(dof, V):
    n = 2*dof
    print("\ncompute_A_inv:\nn = {:d}".format(n))
    A_inv = np.identity(6)

    J = compute_J(n)
    for i in range(n):
        if (i+1) & 1:
            z = 0e0
            for j in range(n):
                for k in range(n):
                    z += V[j][i].real*J[j][k]*V[k][i].imag
            sgn = sign(z)
            z = np.sqrt(np.abs(1e0/z))
            for j in range(n):
                V[j][i] = z*complex(V[j][i].real, sgn*V[j][i].imag)
    for i in range(n):
        A_inv[0][i] = sign(V[0][0].real)*V[i][0].real
        A_inv[1][i] = sign(V[0][0].real)*V[i][0].imag
        A_inv[2][i] = sign(V[2][2].real)*V[i][2].real
        A_inv[3][i] = sign(V[2][2].real)*V[i][2].imag
        if n > 4:
            A_inv[4][i] = sign(V[4][4].real)*V[i][4].real
            A_inv[5][i] = sign(V[4][4].real)*V[i][4].imag
    prt_np_mat("\nA^-1:\n", A_inv)
    return A_inv


def compute_poincare_map(acc, calc_config):
    map = ss_vect_tps()

    map.set_identity()
    acc.propagate(calc_config, map,  0, len(acc))
    print("\nMap:\n", map)
    mat = ss_vect_tps_to_mat(map)
    print("Matrix:\n", mat)
    mat = np.array(mat)
    prt_np_mat("Numpy Matrix:\n", mat)

    return mat


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc         = accelerator_from_config(t_file)
calc_config = ConfigType()

M = compute_poincare_map(acc, calc_config)

prt_np_mat("\nPoincaré Map:\n", M)

print()
dof = 2
compute_nus(dof, M)

rank = np.linalg.matrix_rank(M[0:2*dof, 0:2*dof])
print("\nrank = {:d}".format(rank))

[w, v] = np.linalg.eig(M[0:2*dof, 0:2*dof])
# chop_np_cmplx_mat(v, 1e-13)
prt_np_cmplx_vec("\nEigenvalues:\n", w)
prt_np_cmplx_mat("\nEigenvectors:\n", v)
compute_A_inv(dof, v)

# A = v.real[0:4, 0:4]
# prt_np_mat("\nA:\n", A)
# prt_np_mat("\nA:\n", A*M*np.linalg.inv(A))

exit()

twiss = linear_optics.compute_twiss_parameters(acc, calc_config)
twiss.name = "twiss_parameters"
md = accelerator_info(acc)
md.attrs = dict(calc_config=calc_config)
twiss_with_md = xr.merge([twiss, md])

# print(res)

# combine
