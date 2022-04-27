"""Read lattice file and calculate radiation
"""
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils import linear_optics
from thor_scsi.utils.extract_info import accelerator_info
from thor_scsi.utils.output import vec2txt, mat2txt

from thor_scsi.lib import (
    ConfigType,
    Accelerator,
    ss_vect_tps,
    ss_vect_double,
    RadiationDelegate,
    RadiationDelegateKick,
    phase_space_ind,
    spatial_ind,
    ObservedState,
    StandardObserver,
    ss_vect_tps_to_mat,
)

import os
import sys
import xarray as xr
import numpy as np
import copy
from typing import Sequence

sign = np.sign

X_, Y_, Z_ = spatial_ind.X_, spatial_ind.Y_, spatial_ind.Z_

x_ = phase_space_ind.x_
px_ = phase_space_ind.px_
y_ = phase_space_ind.y_
py_ = phase_space_ind.py_
delta_ = phase_space_ind.delta_
ct_ = phase_space_ind.ct_


def prt_np_vec(name, vec):
    print(name, end="")
    print(vec2txt(vec))


def prt_np_cmplx_vec(name, vec):
    prt_np_vec(name, vec)


def prt_np_mat(name, mat):
    print(name, end="")
    print(mat2txt(mat))


def prt_np_cmplx_mat(name, mat):
    prt_np_mat(name, mat)


def get_mat(t_map):
    mat = ss_vect_tps_to_mat(t_map)
    mat = np.array(mat)
    return mat


def get_map(M):
    [Id, t_map] = [ss_vect_tps(), ss_vect_tps()]
    Id.set_identity()
    t_map.set_zero()
    for j in range(len(M)):
        for k in range(len(M[0])):
            t_map[j] += M[j][k] * Id[k]
    return t_map


def compute_nu(M):
    """
    """
    tr = M.trace()
    if tr < 2e0:
        nu = np.arccos(tr / 2e0) / (2e0 * np.pi)
        if M[0][1] < 0e0:
            nu = 1e0 - nu
        return [nu, True]
    else:
        print("\ncompute_nu: unstable\n")
        return [np.nan, False]


def compute_nus(n_dof, M):
    """
    """
    nus = np.zeros(n_dof)
    stable = [False, False]
    for k in range(n_dof):
        [nus[k], stable[k]] = compute_nu(M[2 * k : 2 * k + 2, 2 * k : 2 * k + 2])
    return [nus, stable]


def compute_nus_symp_mat(n_dof: int, M: np.ndarray) -> (float, bool):
    """Compute nu for a general symplectic periodic transport matrix

    i.e., not assuming mid-plane symmetry.
    """
    n = 2 * n_dof

    nu, tr = np.zeros(2), np.zeros(2)

    # Should be a local copy.
    M1 = copy.copy(M[:n, :n])
    I = np.identity(n)
    detp = np.linalg.det(M1 - I)
    detm = np.linalg.det(M1 + I)
    for i in range(2):
        s = 2 * i
        e = 2 * i + 1
        tr[i] = M[s:e, s:e].trace()

    if tr[X_] > tr[Y_]:
        sgn = 1e0
    else:
        sgn = -1e0

    b = (detp - detm) / 16e0
    c = (detp + detm) / 8e0 - 1e0

    b2mc = b ** 2 - c
    if b2mc < 0e0:
        nu[X_] = nu[Y_] = np.nan
        print("\ncompute_nus_symp_mat: unstable\n")
        return nu, False

    for i in range(2):
        x = -b + sgn * np.sqrt(b2mc) if i == 0 else -b - sgn * np.sqrt(b2mc)
        if abs(x) <= 1e0:
            nu[i] = np.arccos(x) / (2e0 * np.pi)
            if M1[2 * i][2 * i + 1] < 0e0:
                stable = True
                nu[i] = 1e0 - nu[i]
            else:
                nu[i] = np.nan
                plane = ["hor", "vert"][i]
                txt = "fcompute_nus_symp_mat: unstable {plane:%s} plane: x = {x:%10.3e}"
                print("\n" + txt + "\n")

                return nu, False

    return nu, stable


def compute_S(n_dof):
    """Compute the omega matrix

    Required to test the symplectisism for a matrix in block style
    """
    n = 2 * n_dof
    S = np.zeros((n, n))
    for k in range(n_dof):
        s = 2 * k
        e = 2 * k + 1
        S[s, e] = 1e0
        S[e, s] = -1e0
    return S


def swap(w, i, j):
    w[i], w[j] = w[j], w[i]


def swap_imag_part(a: complex, b: complex) -> (complex, complex):
    return a.real + b.imag * 1j, b.real + a.imag * 1j


def swap_imag(w: Sequence, i: int, j: int):
    w[i], w[j] = swap_imag_part(w[i], w[j])


def swap_mat(A, i, j):
    A[:, [i, j]] = A[:, [j, i]]


def swap_mat_imag(A, i, j):
    n = len(A)
    for k in range(n):
        c = A[k][i].imag
        A[k][i] = A[k][i].real + A[k][j].imag * 1j
        A[k][j] = A[k][j].real + c * 1j


def closest(x: float, x1: float, x2: float, x3: float) -> int:
    """Find the value that is closest to x

    Todo:
        check if corner cases need to be checked
        e.g. dx1 == dx2 ...
    """

    dx1, dx2, dx3 = [np.abs(x - tx) for tx in (x1, x2, x3)]
    if dx1 < dx2 and dx1 < dx3:
        k = 1
    elif dx2 < dx1 and dx2 < dx3:
        k = 2
    else:
        k = 3
    return k


def sort_eigen(n_dof, M, w, v):
    """

    Todo:
        can n_dof be derived from the shape of M, w or v

        Check that shape of M, w, v correspond
    """
    n = 2 * n_dof

    sin_M, cos_M = np.zeros(n_dof), np.zeros(n_dof)
    nu1_M, nu2_M, nu1, nu2 = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

    for i in range(n_dof):
        j = (i + 1) * 2 - 1
        cos_M[i] = (M[j - 1][j - 1] + M[j][j]) / 2
        if np.abs(cos_M[i]) > 1e0:
            print(
                "sort_eigen: unstable |cos_M[nu_%d]-1e0| = %10.3e\n",
                i + 1,
                fabs(cos_M[i] - 1e0),
            )
            stable = False
        sin_M[i] = sign(M[j - 1][j]) * np.sqrt(1e0 - cos_M[i] ** 2)
        nu1_M[i] = np.arctan2(sin_M[i], cos_M[i]) / (2 * np.pi)
        if nu1_M[i] < 0.0:
            nu1_M[i] += 1.0
        if nu1_M[i] <= 0.5:
            nu2_M[i] = nu1_M[i]
        else:
            nu2_M[i] = 1.0 - nu1_M[i]
        nu1[i] = np.arctan2(w[j - 1].imag, w[j - 1].real) / (2 * np.pi)
        if nu1[i] < 0.0:
            nu1[i] += 1.0
        if nu1[i] <= 0.5:
            nu2[i] = nu1[i]
        else:
            nu2[i] = 1.0 - nu1[i]

    for i in range(n_dof):
        c = closest(nu2_M[i], nu2[0], nu2[1], nu2[2])
        if c != i + 1:
            j = c * 2 - 1
            k = i * 2 + 1
            swap_mat(v, j - 1, k - 1)
            swap_mat(v, j, k)
            swap(w, j - 1, k - 1)
            swap(w, j, k)
            swap(nu1, i - 1, c - 1)
            swap(nu2, i - 1, c - 1)
    for i in range(n_dof):
        if (0.5 - nu1_M[i]) * (0.5 - nu1[i]) < 0e0:
            j = i * 2 + 1
            swap_mat_imag(v, j - 1, j)
            swap_imag(w, j - 1, j)


def compute_A_inv(n_dof, eta, v):
    """

    Todo:
       whant is special about this inverse
    """
    n = 2 * n_dof

    # Should be a local copy.
    # v1 = np.array(v)
    v1 = copy.copy(v)
    A_inv = np.identity(6)
    S = compute_S(n_dof)

    for i in range(n):
        if (i + 1) % 2 != 0:
            z = np.linalg.multi_dot([v1[:, i].real, S, v1[:, i].imag])
            sgn = sign(z)
            z = np.sqrt(np.abs(1e0 / z))
            for j in range(n):
                v1[j][i] = z * complex(v1[j][i].real, sgn * v1[j][i].imag)

    for i in range(n):
        A_inv[0, i] = sign(v1[0][0].real) * v1[i][0].real
        A_inv[1, i] = sign(v1[0][0].real) * v1[i][0].imag
        A_inv[2, i] = sign(v1[2][2].real) * v1[i][2].real
        A_inv[3, i] = sign(v1[2][2].real) * v1[i][2].imag
        if n > 4:
            A_inv[4, i] = sign(v1[4][4].real) * v1[i][4].real
            A_inv[5, i] = sign(v1[4][4].real) * v1[i][4].imag

    B = np.identity(6)
    B[x_, delta_], B[px_, delta_] = eta[x_], eta[px_]
    B[ct_, x_], B[ct_, px_] = eta[px_], -eta[x_]

    A_inv = np.dot(A_inv, np.linalg.inv(B))
    return [A_inv, v1]


def compute_dnu(n_dof, A):
    """

    a lot of special treatment for the internal representation of the matrix
    """
    eps = 1e-15
    dnu = np.zeros(n_dof)
    for k in range(n_dof):
        if k < 2:
            k2 = 2 * k
            dnu[k] = np.arctan2(A[k2][k2 + 1], A[k2][k2]) / (2e0 * np.pi)
        else:
            dnu[k] = -np.arctan2(A[ct_][delta_], A[ct_][ct_]) / (2e0 * np.pi)
    if dnu[k] < -eps:
        dnu[k] += 1e0
    return dnu


def compute_A_CS(n_dof, A):
    """compute Courant Snyder form of A

    Todo:
        check if dnu and R should not match in shape
    """
    dnu = np.zeros(n_dof)
    # Should that not be n_dof * 2 ?
    R = np.identity(6)

    dnu = compute_dnu(n_dof, A)

    # Build required rotation submatrices
    for k in range(n_dof):
        arg = 2 * np.pi * dnu[k]
        c = np.cos(arg)
        s = np.sin(arg)
        # Rotation matrix
        tmp = np.array([[c, -s], [s, c]])
        k2 = 2 * k
        R[k2 : k2 + 2, k2 : k2 + 2] = tmp
        continue
        R[k2, k2], R[k2][k2 + 1] = c, -s
        R[2 * k + 1][2 * k], R[2 * k + 1][2 * k + 1] = s, c

    Ar = np.dot(A, R)
    return Ar, dnu


def compute_twiss_A(A):
    """
    """
    n_dof = 2
    n = 2 * n_dof

    [eta, alpha, beta] = [np.zeros(n), np.zeros(n_dof), np.zeros(n_dof)]

    for k in range(n_dof):
        k2 = 2 * k
        eta[k2] = A[k2][delta_]
        eta[k2 + 1] = A[k2 + 1][delta_]
        alpha[k] = -(A[k2][k2] * A[k2 + 1][k2] + A[k2][k2 + 1] * A[k2 + 1][k2 + 1])
        beta[k] = A[k2][k2] ** 2 + A[k2][k2 + 1] ** 2
    dnu = compute_dnu(n_dof, A)

    return eta, alpha, beta, dnu


def compute_twiss_A_A_tp(A):
    """

    Todo:
      difference ti compute_twiss_A?
    """
    n_dof = 2
    n = 2 * n_dof

    eta, alpha, beta = np.zeros(n), np.zeros(n_dof), np.zeros(n_dof)

    A_A_tp = np.dot(A[0:n, 0:n], np.transpose(A[0:n, 0:n]))
    for k in range(n_dof):
        k2 = 2 * k
        eta[k2], eta[k2 + 1] = A[k2][delta_], A[k2 + 1][delta_]
        alpha[k] = -A_A_tp[k2][k2 + 1]
        beta[k] = A_A_tp[k2][k2]
    dnu = compute_dnu(n_dof, A)

    return eta, alpha, beta, dnu


def compute_dispersion(M):
    """
    """
    n = 4
    I = np.identity(n)
    D = M[:n, delta_]
    return np.dot(np.linalg.inv(I - M[:n, :n]), D)


def compute_twiss_M(M):
    n_dof = 2

    alpha, beta, nu = np.zeros(n_dof), np.zeros(n_dof), np.zeros(n_dof)

    eta = compute_dispersion(M)

    stable = [True, True]
    for k in range(n_dof):
        k2 = 2 * k
        cos = M[k2 : k2 + 2, k2 : k2 + 2].trace() / 2e0
        if abs(cos) >= 1e0:
            print("\ncompute_twiss_M: {:5.3f}\n".format(cos))
            stable[k] = False
        sin = np.sqrt(1e0 - cos ** 2) * sign(M[k2][k2 + 1])
        alpha[k] = (M[k2][k2] - M[k2 + 1][k2 + 1]) / (2e0 * sin)
        beta[k] = M[k2][k2 + 1] / sin
        nu[k] = np.arctan2(sin, cos) / (2e0 * np.pi)
        if nu[k] < 0e0:
            nu[k] += 1e0

    return [eta, alpha, beta, nu, stable]


def compute_map(acc, calc_config):
    """Propaagate an identity map through the accelerator
    """
    t_map = ss_vect_tps()
    t_map.set_identity()
    # acc.propagate(calc_config, t_map, 0, len(acc))
    acc.propagate(calc_config, t_map)
    return t_map


def _extract_tps(elem):
    """Extract tps data from the elments' observer
    """
    ob = elem.getObserver()
    assert ob.hasTruncatedPowerSeries()
    return ob.getTruncatedPowerSeries()


def tps2twiss(tps) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
    """

    returns eta, alpha, beta, nu

    """
    A_mat = get_mat(tps)[:6, :6]
    return compute_twiss_A(A_mat)


def compute_twiss_along_lattice(
    acc: Accelerator, calc_config: ConfigType, A: ss_vect_tps
) -> xr.Dataset:
    """

    Args:
        acc : an :class:`tlib.Accelerator` instance
        calc_config : an :class:`tlib.Config` instance
        A :


    returns xr.Dataset

    Todo:
        Consider if tps should be returned too
    """

    # Instrument it with observers ... I guess I have to keep them here so that
    # internally these are not weak references ... to be tested
    #
    observers = [StandardObserver() for elem in acc]
    # Register them to the elments
    for ob, elem in zip(observers, acc):
        elem.setObserver(ob)

    # Propaagate through the accelerator
    acc.propagate(calc_config, A)

    # Retrieve information

    tps_tmp = [_extract_tps(elem) for elem in acc]
    data = [tps2twiss(t) for t in tps_tmp]
    twiss_pars = [d[1:] for d in data]
    disp = [d[0] for d in data]
    indices = [elem.index for elem in acc]

    tps_tmp = np.array(tps_tmp)
    phase_space_coords_names = ["x", "px", "y", "py", "delta", "ct"]

    tps = xr.DataArray(
        data=tps_tmp,
        name="tps",
        dims=["index", "phase_coordinate"],
        coords=[indices, phase_space_coords_names],
    )
    dispersion = xr.DataArray(
        data=disp,
        name="dispersion",
        dims=["index", "phase_coordinate"],
        coords=[indices, ["x", "px", "y", "py"]],
    )
    twiss_parameters = xr.DataArray(
        twiss_pars,
        name="twiss",
        dims=["index", "par", "plane"],
        coords=[indices, ["alpha", "beta", "dnu"], ["x", "y"]],
    )
    res = xr.merge([twiss_parameters, dispersion, tps])
    return res


def _compute_twiss_lattice(acc, calc_config, A):
    n_dof = 2

    nu = np.zeros(n_dof)

    s = 0
    # Should be a local copy (ss_vect<tps>).
    Ak = A
    print(
        "\n     name            s       alpha_x   beta_x     nu_x"
        "      eta_x    etap_x    alpha_y   beta_y     nu_y      eta_y"
        "    etap_y"
    )
    print(
        "                    [m]                 [m]                  [m]"
        "                          [m]                  [m]"
    )
    for k in range(len(acc)):
        acc.propagate(calc_config, Ak, k, k + 1)
        A_mat = get_mat(Ak)[0:6, 0:6]
        [eta, alpha, beta, dnu] = compute_twiss_A(A_mat)
        nu += dnu
        A_mat = compute_A_CS(2, A_mat)[0]
        Ak = get_map(A_mat)

        print(
            "{:4d} {:10s} {:9.5f} {:9.5f} {:9.5f} {:9.5f} {:9.5f} {:9.5f}"
            " {:9.5f} {:9.5f} {:9.5f} {:9.5f} {:9.5f}".format(
                k,
                acc[k].name,
                s,
                alpha[X_],
                beta[X_],
                nu[X_],
                eta[x_],
                eta[px_],
                alpha[Y_],
                beta[Y_],
                nu[Y_],
                eta[y_],
                eta[py_],
            )
        )


def compute_twiss_lattice(file_name, acc, calc_config, A):
    """
    """
    compute_twiss_along_lattice(acc, calc_config, A)

    stdout = sys.stdout
    try:
        sys.stdout = open(file_name, "w")
        r = _compute_twiss_lattice(acc, calc_config, A)
    finally:
        sys.stdout = stdout

    return r


def compute_ring_twiss(M):
    n_dof = 2
    n = 2 * n_dof

    M_tp = np.transpose(M[0:n, 0:n])
    [w, v] = np.linalg.eig(M_tp)
    sort_eigen(n_dof, M_tp, w, v)
    eta = compute_dispersion(M)
    [A_inv, v1] = compute_A_inv(n_dof, eta, v)

    A = np.linalg.inv(A_inv)
    A = compute_A_CS(2, A)[0]

    return A


def prt_twiss(eta, alpha, beta, nu):
    txt = (
        f"\n  eta   = [{eta[x_]:5.3f}, {eta[px_]:5.3f}, {eta[y_]:5.3f}, {eta[py_]:5.3f}]"
        f"\n  alpha = [{alpha[X_]:5.3f}, {alpha[Y_]:5.3f}]"
        f"\n  beta  = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]"
        f"\n  nu    = [{nu[X_]:5.3f}, {nu[Y_]:5.3f}]"
    )
    print(txt)


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc = accelerator_from_config(t_file)
calc_config = ConfigType()

if True:
    # propagate through the accelerator
    t_map = compute_map(acc, calc_config)
    M = get_mat(t_map)
else:
    # the lattice given is the lattice for a ring
    # so one can start at any point and pass it around
    n = 15
    t_map = ss_vect_tps()
    t_map.set_identity()
    acc.propagate(calc_config, t_map, n)
    acc.propagate(calc_config, t_map, 0, n)
    M = get_mat(t_map)

# Reduce matrix from 7x7 to 6x6.
M = M[:6, :6]
prt_np_mat("\nPoincaré T_Map:\n", M)

A = compute_ring_twiss(M)
prt_np_mat("\nA:\n", A)

n_dof = 2
n = 2 * n_dof

nus, stable = compute_nus(n_dof, M)
print("\ncompute_nus:\n  nu    = [{:5.3f}, {:5.3f}]".format(nus[X_], nus[Y_]))

nus, stable = compute_nus_symp_mat(n_dof, M)
print("\ncompute_nus_symp_mat:\n  nu    = [{:5.3f}, {:5.3f}]".format(nus[X_], nus[Y_]))

eta, alpha, beta, nu, stable = compute_twiss_M(M)
print("\ncompute_twiss_M:")
prt_twiss(eta, alpha, beta, nu)

eta, alpha, beta, nu = compute_twiss_A(A)
print("\ncompute_twiss_A:")
prt_twiss(eta, alpha, beta, nu)
eta, alpha, beta, nu = compute_twiss_A_A_tp(A)
print("\ncompute_twiss_A_A_tp:")
prt_twiss(eta, alpha, beta, nu)

# Cross check.
prt_np_mat("\nA^-1*M*A:\n", np.linalg.multi_dot([np.linalg.inv(A), M, A]))

compute_twiss_lattice("linlat.out", acc, calc_config, get_map(A))

# Swap columns.
# v[:, [1, 0]] = v[:, [0, 1]]
# Swap rows.
# v[[1, 0]] = v[[0, 1]]

# mat = np.zeros((6, 6))

# exit()

## twiss = linear_optics.compute_twiss_parameters(acc, calc_config)
## twiss.name = "twiss_parameters"
## md = accelerator_info(acc)
## md.attrs = dict(calc_config=calc_config)
## twiss_with_md = xr.merge([twiss, md])

# print(res)

# combine
