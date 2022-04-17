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

[X_, Y_, Z_]                    = [0, 1, 2]
[x_, px_, y_, py_, ct_, delta_] = [0, 1, 2, 3, 5, 4]


def prt_np_cmplx_vec(str, vec):
    print(str, end="")
    for k in range(len(vec)):
        print("{:12.3e} {:} {:9.3e}i".
              format(vec[k].real, "+-"[int(vec[k].imag < 0)], abs(vec[k].imag)),
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


def compute_nu(mat):
    tr = mat.trace()
    if tr < 2e0:
        stable = True
    else:
        printf("\ncompute_nu: unstable\n")
        stable = False
    nu = np.arccos(tr/2e0)/(2e0*np.pi)
    if mat[0][1] < 0e0:
        nu = 1e0 - nu
    return [nu, stable]


def compute_nus(n_dof, mat):
    nus = np.zeros(n_dof)
    stable = [False, False]
    for k in range(n_dof):
        [nus[k], stable[k]] = compute_nu(mat[2*k:2*k+2, 2*k:2*k+2])
    return [nus, stable]


def compute_nus_symp_mat(n_dof, M):
    # Compute nu for a general symplectic periodic transport matrix;
    # i.e., not assuming mid-plane symmetry.
    n = 2*n_dof
    nu = np.zeros(2)
    tr = np.zeros(2)
    M1 = M
    for i in range(n):
        M1[i][i] -= 1e0
    detp = np.linalg.det(M1)
    for i in range(n):
        M1[i][i] += 2e0
    detm = np.linalg.det(M1)
    for i in range(n):
        M1[i][i] -= 1e0
    for i in range(2):
        tr[i] = M1[2*i][2*i] + M1[2*i+1][2*i+1]
    sgn = 1e0 if tr[X_] > tr[Y_] else -1e0
    b = (detp-detm)/16e0
    c = (detp+detm)/8e0 - 1e0
    b2mc = b**2 - c
    if b2mc < 0e0:
        stable = False
        nu[X_] = nu[Y_] = NAN
        printf("\nGetNu: unstable\n")
        return [nu, stable]

    for i in range(2):
        x = -b + sgn*np.sqrt(b2mc) if i == 0 else -b - sgn*np.sqrt(b2mc)
        if abs(x) <= 1e0:
            nu[i] = np.arccos(x)/(2e0*np.pi)
            if M1[2*i][2*i+1] < 0e0:
                stable = True
                nu[i] = 1e0 - nu[i]
            else:
                stable = False
                nu[i] = NAN
                print("\ncompute_nu: unstable {:%s} plane {:%10.3e}\n".
                      format("hor" if i == 0 else "ver", x))
                return [nu, stable]

    return [nu, stable]


def sign(x):
    if x > 0e0:
        return 1
    elif x < 0e0:
        return -1
    else:
        return 0


def compute_S(n_dof):
    n = 2*n_dof
    S = np.zeros((n, n))
    for k in range(n_dof):
        S[2*k][2*k+1] =  1e0
        S[2*k+1][2*k] = -1e0
    return S


def swap(x, y):
   return [y, x]


def swap_mat(n, A, i, j):
    print("\nswap_mat:\n", n, i, j)
    if True:
        for k in range(n):
            x = A[k][i-1]
            A[k][i-1] = A[k][j-1]
            A[k][j-1] = x
    else:
        v[:, [j-1, i-1]] = v[:, [i-1, j-1]]
    return A


def closest(x, x1, x2, x3):
    dx1 = np.abs(x-x1)
    dx2 = np.abs(x-x2)
    dx3 = np.abs(x-x3)
    if (dx1 < dx2) and (dx1 < dx3):
        k = 1
    elif (dx2 < dx1) and (dx2 < dx3):
        k = 2
    else:
        k = 3
    return k;


def sort_eigen(n_dof, M, w, v):
    n = 2*n_dof
    sin_M = np.zeros(n_dof)
    cos_M = np.zeros(n_dof)
    nu1_M = np.zeros(3)
    nu2_M = np.zeros(3)
    nu1   = np.zeros(3)
    nu2   = np.zeros(3)
    for i in range(n_dof):
        j = (i+1)*2 - 1
        cos_M[i] = (M[j-1][j-1]+M[j][j])/2
        if np.abs(cos_M[i]) > 1e0:
            print("geigen: unstable |cos_M[nu_%d]-1e0| = %10.3e\n",
	          i+1, fabs(cos_M[i]-1e0))
            stable = False
        sin_M[i] = sign(M[j-1][j])*np.sqrt(1e0-cos_M[i]**2)
        nu1_M[i] = np.arctan2(sin_M[i], cos_M[i])/(2e0*np.pi)
        if nu1_M[i] < 0e0:
            nu1_M[i] += 1e0
        if nu1_M[i] <= 0.5e0:
            nu2_M[i] = nu1_M[i]
        else:
            nu2_M[i] = 1e0 - nu1_M[i];
        nu1[i] = np.arctan2(w[j-1].imag, w[j-1].real)/(2e0*np.pi)
        if nu1[i] < 0e0:
            nu1[i] += 1e0
        if nu1[i] <= 0.5e0:
            nu2[i] = nu1[i]
        else:
            nu2[i] = 1e0 - nu1[i]
    for i in range(n_dof):
        c = closest(nu2_M[i], nu2[0], nu2[1], nu2[2])
        if c != i+1:
            j = c*2 - 1
            k = i*2 + 1
            v = swap_mat(n, v, j,   k)
            v = swap_mat(n, v, j+1, k+1)
            [w[j-1], w[k-1]] = swap(w[j-1], w[k-1])
            [w[j], w[k]]     = swap(w[j],   w[k])
            [nu1[i+1-1], nu1[c-1]] = swap(nu1[i+1-1], nu1[c-1])
            [nu2[i+1-1], nu2[c-1]] = swap(nu2[i+1-1], nu2[c-1])
    for i in range(n_dof):
        if (0.5e0-nu1_M[i])*(0.5e0-nu1[i]) < 0e0:
            j = i*2
            for k in range(n):
                c = v[k, j+1]
                v[k][j+1] = complex(v[k][j+1].real, v[k][j].imag)
                v[k][j]   = complex(v[k][j].real,   c.imag)
            c = w[j+1]
            w[j+1] = complex(w[j+1].real, w[j].imag)
            w[j]   = complex(w[j].real,   c.imag)
    return [w, v]


def compute_A_inv(n_dof, eta, v):
    n = 2*n_dof
    v1 = np.array(v)
    A_inv = np.identity(6)
    S = compute_S(n_dof)

    for i in range(n):
        if (i+1) % 2 != 0:
            z = 0e0
            for j in range(n):
                for k in range(n):
                    z += v1[j][i].real*S[j][k]*v1[k][i].imag
            sgn = sign(z)
            z = np.sqrt(np.abs(1e0/z))
            for j in range(n):
                v1[j][i] = z*complex(v1[j][i].real, sgn*v1[j][i].imag)

    for i in range(n):
        A_inv[0][i] = sign(v1[0][0].real)*v1[i][0].real
        A_inv[1][i] = sign(v1[0][0].real)*v1[i][0].imag
        A_inv[2][i] = sign(v1[2][2].real)*v1[i][2].real
        A_inv[3][i] = sign(v1[2][2].real)*v1[i][2].imag
        if n > 4:
            A_inv[4][i] = sign(v1[4][4].real)*v1[i][4].real
            A_inv[5][i] = sign(v1[4][4].real)*v1[i][4].imag

    B = np.identity(6)
    B[x_, delta_]  = eta[x_]
    B[px_, delta_] = eta[px_]
    B[ct_, x_]     = eta[px_]
    B[ct_, px_]    = -eta[x_]

    A_inv = np.dot(A_inv, np.linalg.inv(B))
    return [A_inv, v1]


def compute_dnu(n_dof, A):
    eps = 1e-15
    dnu = np.zeros(n_dof)
    for k in range(n_dof):
        if k < 2:
            dnu[k] = np.arctan2(A[2*k][2*k+1], A[2*k][2*k])/(2e0*np.pi)
        else:
            dnu[k] = -np.arctan2(A[ct_][delta_], A[ct_][ct_])/(2e0*np.pi)
    if dnu[k] < -eps:
        dnu[k] += 1e0
    return dnu


def compute_A_CS(n_dof, A):
    dnu = np.zeros(n_dof)
    R   = np.identity(6)
    dnu = compute_dnu(n_dof, A)
    for k in range(n_dof):
        c = np.cos(2e0*np.pi*dnu[k])
        s = np.sin(2e0*np.pi*dnu[k])
        R[2*k][2*k]     = c
        R[2*k][2*k+1]   = -s
        R[2*k+1][2*k]   = s
        R[2*k+1][2*k+1] = c
    return [np.dot(A, R), dnu]


def compute_twiss_A(A):
    n_dof = 2
    eta   = np.zeros(n_dof)
    etap  = np.zeros(n_dof)
    alpha = np.zeros(n_dof)
    beta  = np.zeros(n_dof)
    for k in range(n_dof):
        eta[k]   = A[2*k][delta_]
        etap[k]  = A[2*k+1][delta_]
        alpha[k] = -A[2*k][2*k]*A[2*k+1][2*k] + A[2*k][2*k+1]*A[2*k+1][2*k+1]
        beta[k]  = A[2*k][2*k]**2 + A[2*k][2*k+1]**2
    print("\neta   = [{:6.3f}, {:6.3f}]\netap  = [{:6.3f}, {:6.3f}]\n"
          "alpha = [{:6.3f}, {:6.3f}]\nbeta  = [{:6.3f}, {:6.3f}]".
          format(eta[X_], eta[Y_], etap[X_], etap[Y_], alpha[X_], alpha[Y_],
                 beta[X_], beta[Y_]))
    return [eta, etap, alpha, beta]


def compute_twiss_A2(A):
    n_dof = 2
    n     = 2*n_dof
    eta   = np.zeros(n_dof)
    etap  = np.zeros(n_dof)
    alpha = np.zeros(n_dof)
    beta  = np.zeros(n_dof)
    dnu   = np.zeros(n_dof)
    A_A_tp = np.dot(A[0:n, 0:n], np.transpose(A[0:n, 0:n]))
    for k in range(n_dof):
        eta[k]   = A[2*k][delta_]
        etap[k]  = A[2*k+1][delta_]
        alpha[k] = -A_A_tp[2*k][2*k+1]
        beta[k] = A_A_tp[2*k][2*k]
    dnu = compute_dnu(n_dof, A)
    print("\neta   = [{:6.3f}, {:6.3f}]\netap  = [{:6.3f}, {:6.3f}]\n"
          "alpha = [{:6.3f}, {:6.3f}]\nbeta  = [{:6.3f}, {:6.3f}]".
          format(eta[X_], eta[Y_], etap[X_], etap[Y_], alpha[X_], alpha[Y_],
                 beta[X_], beta[Y_]))
    return [eta, etap, alpha, beta, dnu]


def compute_disp(M):
    n = 2
    I = np.identity(n)
    D = np.array((M[x_][delta_], M[px_][delta_]))
    return np.dot(np.linalg.inv(I-M[0:n, 0:n]), D)


def compute_twiss_M(M):
    n_dof = 2
    eta   = np.zeros(n_dof)
    alpha = np.zeros(n_dof)
    beta  = np.zeros(n_dof)
    eta = compute_disp(M)
    for k in range(n_dof):
        cos = M[2*k:2*k+2, 2*k:2*k+2].trace()/2e0
        if abs(cos) >= 1e0:
            print("\ncompute_twiss_M: {:5.3f}\n".format(cos))
            exit(1)
        sin      = np.sqrt(1e0-cos**2)*sign(M[2*k][2*k+1])
        alpha[k] = (M[2*k][2*k]-M[2*k+1][2*k+1])/(2e0*sin)
        beta[k]  = M[2*k][2*k+1]/sin
    [nu, stable] = compute_nus_symp_mat(n_dof, M)
    print("\neta   = [{:6.3f}, {:6.3f}]\nalpha = [{:6.3f}, {:6.3f}]\n"
          "beta  = [{:6.3f}, {:6.3f}]\nnu    = [{:6.3f}, {:6.3f}]".
          format(eta[0], eta[1], alpha[X_], alpha[Y_], beta[X_], beta[Y_],
                 nu[X_], nu[Y_]))
    return [eta, alpha, beta, nu]


def compute_map(acc, calc_config):
    map = ss_vect_tps()
    map.set_identity()
    acc.propagate(calc_config, map,  0, len(acc))
    return map


def get_mat(map):
    mat = ss_vect_tps_to_mat(map)
    mat = np.array(mat)
    return mat


def get_map(M):
    Id  = ss_vect_tps()
    map = ss_vect_tps()
    map.set_zero()
    Id.set_identity()
    for j in range(len(M)):
        for k in range(len(M[0])):
            map[j] += M[j][k]*Id[k] 
    return map


def compute_twiss_lat(acc, calc_config, A):
    A1 = A
    for k in range(len(acc)):
        acc.propagate(calc_config, A1, k, k+1)
        A_mat = get_mat(A1)
        prt_np_mat("\nA:\n", A_mat)
        compute_twiss_A(A1)


t_dir  = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc         = accelerator_from_config(t_file)
calc_config = ConfigType()

if not True:
    map = compute_map(acc, calc_config)
    M = get_mat(map)
else:
    n = 15
    map = ss_vect_tps()
    map.set_identity()
    acc.propagate(calc_config, map,  n, len(acc))
    acc.propagate(calc_config, map,  0, n)
    M = get_mat(map)

# Reduce matrix from 7x7 to 6x6.
M = M[0:6, 0:6]
prt_np_mat("\nPoincaré Map:\n", M)

n_dof = 2
n     = 2*n_dof

[nus, stable] = compute_nus(n_dof, M)
print("\nnu = [{:5.3f}, {:5.3f}]".format(nus[X_], nus[Y_]))

[eta, alpha, beta, nu] = compute_twiss_M(M)

M_tp = np.transpose(M)
[w, v] = np.linalg.eig(M_tp[0:n, 0:n])
prt_np_cmplx_vec("\neigenvalues:\n", w)
prt_np_cmplx_mat("\neigenvectors:\n", v)
[w, v] = sort_eigen(n_dof, M, w, v)
prt_np_cmplx_vec("\neigenvalues:\n", w)
prt_np_cmplx_mat("\neigenvectors:\n", v)
[A_inv, v1] = compute_A_inv(n_dof, eta, v)

A = np.linalg.inv(A_inv)
A = compute_A_CS(2, A)[0]
prt_np_mat("\nA:\n", A)

compute_twiss_A2(A)
compute_twiss_A(A)

# Cross check.
prt_np_mat("\nA^-1*M*A:\n",
           np.linalg.multi_dot([np.linalg.inv(A), M, A]))

A_map = get_map(A)
print("\nA:\n", A_map)
compute_twiss_lat(acc, calc_config, A_map)

# Swap columns.
# v[:, [1, 0]] = v[:, [0, 1]]
# Swap rows.
# v[[1, 0]] = v[[0, 1]]

# mat = np.zeros((6, 6))

exit()

twiss = linear_optics.compute_twiss_parameters(acc, calc_config)
twiss.name = "twiss_parameters"
md = accelerator_info(acc)
md.attrs = dict(calc_config=calc_config)
twiss_with_md = xr.merge([twiss, md])

# print(res)

# combine
