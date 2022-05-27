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
#    phase_space_indicator,
    ObservedState,
    ss_vect_tps_to_mat,
    partialInverse,
    xabs
)

import os
import sys
import xarray as xr
import numpy as np

[X_, Y_, Z_]                    = [0, 1, 2]
[x_, px_, y_, py_, ct_, delta_] = [0, 1, 2, 3, 5, 4]
[Dip, Quad, Sext]               = [1, 2, 3]

ss_dim = 6


def get_mat(map):
    mat = ss_vect_tps_to_mat(map)
    mat = np.array(mat)
    return mat


def prt_np_vec(str, fmt, vec):
    print(str, end="")
    for k in range(len(vec)):
        print(f"{vec[k]:{fmt}}", end="")
    print()


def prt_np_cmplx_vec(str, vec):
    print(str, end="")
    for k in range(len(vec)):
        print("{:12.3e} {:} {:9.3e}i".
              format(vec[k].real, "+-"[int(vec[k].imag < 0)], abs(vec[k].imag)),
              end="")
    print()


def prt_np_mat(str, fmt, mat):
    print(str, end="")
    for k in range(len(mat)):
        prt_np_vec("", fmt, mat[k])


def prt_np_cmplx_mat(str, mat):
    print(str, end="")
    for k in range(len(mat)):
        prt_np_cmplx_vec("", mat[k])


def prt_ps_vec(str, ps):
    n = 6
    print(str, end="")
    for k in range(n):
        print(f"{ps[k]:11.3e}", end="")
    print()


def prt_map(str, fmt, map):
    prt_np_mat(str, fmt, get_mat(map))


def get_map(M):
    [Id, map] = [ss_vect_tps(), ss_vect_tps()]
    Id.set_identity()
    map.set_zero()
    for j in range(len(M)):
        for k in range(len(M[0])):
            map[j] += M[j][k]*Id[k]
    return map


def compute_nu(M):
    tr = M.trace()
    if tr < 2e0:
        nu = np.arccos(tr/2e0)/(2e0*np.pi)
        if M[0][1] < 0e0:
            nu = 1e0 - nu
        return [nu, True]
    else:
        printf("\ncompute_nu: unstable\n")
        return [NAN, False]


def compute_nus(n_dof, M):
    nus    = np.zeros(n_dof)
    stable = [False, False]
    for k in range(n_dof):
        [nus[k], stable[k]] = compute_nu(M[2*k:2*k+2, 2*k:2*k+2])
    return [nus, stable]


def compute_nus_symp_mat(n_dof, M):
    # Compute nu for a general symplectic periodic transport matrix;
    # i.e., not assuming mid-plane symmetry.
    n = 2*n_dof

    [nu, tr] = [np.zeros(2), np.zeros(2)]

    # Should be a local copy.
    M1 = M[0:n, 0:n]
    I = np.identity(n)
    detp = np.linalg.det(M1-I)
    detm = np.linalg.det(M1+I)
    for i in range(2):
        tr[i] = M[2*i:2*i+2, 2*i:2*i+2].trace()
    sgn = 1e0 if tr[X_] > tr[Y_] else -1e0
    b = (detp-detm)/16e0
    c = (detp+detm)/8e0 - 1e0
    b2mc = b**2 - c
    if b2mc < 0e0:
        nu[X_] = nu[Y_] = NAN
        printf("\ncompute_nus_symp_mat: unstable\n")
        return [nu, False]

    for i in range(2):
        x = -b + sgn*np.sqrt(b2mc) if i == 0 else -b - sgn*np.sqrt(b2mc)
        if abs(x) <= 1e0:
            nu[i] = np.arccos(x)/(2e0*np.pi)
            if M1[2*i][2*i+1] < 0e0:
                stable = True
                nu[i] = 1e0 - nu[i]
            else:
                nu[i] = NAN
                print("\ncompute_nus_symp_mat: unstable {:%s} plane {:%10.3e}\n"
                      .format("hor" if i == 0 else "ver", x))
                return [nu, False]

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


def swap(w, i, j):
    [w[i], w[j]] = [w[j], w[i]]


def swap_imag(w, i, j):
    [w[i], w[j]] = \
        [complex(w[i].real, w[j].imag), complex(w[j].real, w[i].imag)]


def swap_mat(A, i, j):
    A[:, [i, j]] = A[:, [j, i]]


def swap_mat_imag(A, i, j):
    n = len(A)
    for k in range(n):
        c = A[k][i].imag
        A[k][i] = complex(A[k][i].real, A[k][j].imag)
        A[k][j] = complex(A[k][j].real, c)


def closest(x, x1, x2, x3):
    [dx1, dx2, dx3] = [np.abs(x-x1), np.abs(x-x2), np.abs(x-x3)]
    if (dx1 < dx2) and (dx1 < dx3):
        k = 1
    elif (dx2 < dx1) and (dx2 < dx3):
        k = 2
    else:
        k = 3
    return k;


def sort_eigen(n_dof, M, w, v):
    n = 2*n_dof

    [sin_M, cos_M]           = [np.zeros(n_dof), np.zeros(n_dof)]
    [nu1_M, nu2_M, nu1, nu2] = \
        [np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)]

    for i in range(n_dof):
        j = (i+1)*2 - 1
        cos_M[i] = (M[j-1][j-1]+M[j][j])/2
        if np.abs(cos_M[i]) > 1e0:
            print("sort_eigen: unstable |cos_M[nu_%d]-1e0| = %10.3e\n",
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
            swap_mat(v, j-1, k-1)
            swap_mat(v, j,   k)
            swap(w, j-1, k-1)
            swap(w, j,   k)
            swap(nu1, i-1, c-1)
            swap(nu2, i-1, c-1)
    for i in range(n_dof):
        if (0.5e0-nu1_M[i])*(0.5e0-nu1[i]) < 0e0:
            j = i*2 + 1
            swap_mat_imag(v, j-1, j)
            swap_imag(w, j-1, j)


def compute_A_inv(n_dof, eta, v):
    n = 2*n_dof

    # Should be a local copy.
    v1 = np.array(v)
    A_inv = np.identity(6)
    S = compute_S(n_dof)

    # Normalize the eigenvectors.
    res = np.linalg.multi_dot([np.transpose(v1), S, v1])
    for k in range(2):
        res_12 = res[2*k][2*k+1].imag
        scl = sign(res_12)*np.sqrt(np.abs(res_12/2e0))
        v1[:, 2*k] /= scl
        v1[:, 2*k+1] /= scl

    # Construct A^-1 from the real & imaginary parts of the eigenvectors.
    for k in range(n_dof):
        A_inv[2*k, 0:n]   = sign(v1[2*k][2*k].real)*v1[0:n, 2*k].real
        A_inv[2*k+1, 0:n] = sign(v1[2*k][2*k].real)*v1[0:n, 2*k].imag

    # Include dispersion; i.e., transform to the momentum dependent fix point.
    B = np.identity(6)
    [B[x_, delta_], B[px_, delta_]] = [eta[x_],   eta[px_]]
    [B[ct_, x_],    B[ct_, px_]]    = [eta[px_], -eta[x_]]

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
    [dnu, R] = [np.zeros(n_dof), np.identity(6)]

    dnu = compute_dnu(n_dof, A)

    for k in range(n_dof):
        [c, s] = [np.cos(2e0*np.pi*dnu[k]), np.sin(2e0*np.pi*dnu[k])]
        [R[2*k][2*k], R[2*k][2*k+1]]     = [c, -s]
        [R[2*k+1][2*k], R[2*k+1][2*k+1]] = [s,  c]

    return [np.dot(A, R), dnu]


def compute_twiss_A(A):
    n_dof = 2
    n     = 2*n_dof

    [eta, alpha, beta] = [np.zeros(n), np.zeros(n_dof), np.zeros(n_dof)]

    for k in range(n_dof):
        [eta[2*k], eta[2*k+1]] = [A[2*k][delta_], A[2*k+1][delta_]]
        alpha[k] = -(A[2*k][2*k]*A[2*k+1][2*k]+A[2*k][2*k+1]*A[2*k+1][2*k+1])
        beta[k]  = A[2*k][2*k]**2 + A[2*k][2*k+1]**2
    dnu = compute_dnu(n_dof, A)

    return [eta, alpha, beta, dnu]


def compute_twiss_A_A_tp(A):
    n_dof = 2
    n     = 2*n_dof

    [eta, alpha, beta] = [np.zeros(n), np.zeros(n_dof), np.zeros(n_dof)]

    A_A_tp = np.dot(A[0:n, 0:n], np.transpose(A[0:n, 0:n]))
    for k in range(n_dof):
        [eta[2*k], eta[2*k+1]] = [A[2*k][delta_], A[2*k+1][delta_]]
        alpha[k] = -A_A_tp[2*k][2*k+1]
        beta[k]  = A_A_tp[2*k][2*k]
    dnu = compute_dnu(n_dof, A)

    return [eta, alpha, beta, dnu]


def compute_disp(M):
    n = 4
    I = np.identity(n)
    D = M[0:n, delta_]
    return np.dot(np.linalg.inv(I-M[0:n, 0:n]), D)


def compute_twiss_M(M):
    n_dof = 2

    [alpha, beta, nu] = [np.zeros(n_dof), np.zeros(n_dof), np.zeros(n_dof)]

    eta = compute_disp(M)

    stable = [True, True]
    for k in range(n_dof):
        cos = M[2*k:2*k+2, 2*k:2*k+2].trace()/2e0
        if abs(cos) >= 1e0:
            print(f"\ncompute_twiss_M: {cos:5.3f}")
            stable[k] = False
        sin      = np.sqrt(1e0-cos**2)*sign(M[2*k][2*k+1])
        alpha[k] = (M[2*k][2*k]-M[2*k+1][2*k+1])/(2e0*sin)
        beta[k]  = M[2*k][2*k+1]/sin
        nu[k] = np.arctan2(sin, cos)/(2e0*np.pi)
        if nu[k] < 0e0:
            nu[k] += 1e0

    return [eta, alpha, beta, nu, stable]


def compute_map(acc, calc_config):
    map = ss_vect_tps()
    map.set_identity()
    acc.propagate(calc_config, map,  0, len(acc))
    return map


def compute_twiss_lat(file_name, acc, calc_config, A):
    n_dof = 2

    nu = np.zeros(n_dof)

    stdout     = sys.stdout
    sys.stdout = open(file_name, 'w')

    s = 0
    # Should be a local copy (ss_vect<tps>).
    Ak = A
    print("\n     name            s       alpha_x   beta_x     nu_x"
          "      eta_x    etap_x    alpha_y   beta_y     nu_y      eta_y"
          "    etap_y")
    print("                    [m]                 [m]                  [m]"
          "                          [m]                  [m]")
    for k in range(len(acc)):
        s += acc[k].getLength()
        acc.propagate(calc_config, Ak, k, k+1)
        A_mat = get_mat(Ak)[0:6, 0:6]
        [eta, alpha, beta, dnu] = compute_twiss_A(A_mat)
        nu += dnu
        A_mat = compute_A_CS(2, A_mat)[0]
        Ak = get_map(A_mat)

        print(f"{k:4d} {acc[k].name:10s} {s:9.5f} {alpha[X_]:9.5f}"
              " {beta[X_]:9.5f} {nu[X_]:9.5f} {eta[x_]:9.5f} {eta[px_]:9.5f}"
              " {alpha[Y_]:9.5f} {beta[Y_]:9.5f} {nu[Y_]:9.5f} {eta[y_]:9.5f}"
              " {eta[py_]:9.5f}")

    sys.stdout = stdout


def compute_ring_twiss(M):
    n_dof = 2
    n     = 2*n_dof

    M_tp = np.transpose(M[0:n, 0:n])
    [w, v] = np.linalg.eig(M_tp)
    sort_eigen(n_dof, M_tp, w, v)

    eta = compute_disp(M)
    [A_inv, v1] = compute_A_inv(n_dof, eta, v)

    A = np.linalg.inv(A_inv)
    A = compute_A_CS(2, A)[0]

    return A


def prt_twiss(str, eta, alpha, beta, nu):
    print(f"{str:s}"
          "  eta   = [{eta[x_]:5.3f}, {eta[px_]:5.3f}"
          ", {eta[y_]:5.3f}, {eta[py_]:5.3f}]"
          "\n  alpha = [{alpha[X_]:5.3f}, {alpha[Y_]:5.3f}]"
          "\n  beta  = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]"
          "\n  nu    = [{nu[X_]:5.3f}, {nu[Y_]:5.3f}]")


def test_stuff(n):
    map = ss_vect_tps()

    map.set_identity()
    acc.propagate(calc_config, map,  n, len(acc))
    acc.propagate(calc_config, map,  0, n)
    M = get_mat(map)
    # Reduce matrix from 7x7 to 6x6.
    M = M[0:6, 0:6]
    prt_np_mat("\nPoincaré Map:\n", "14.6e", M)

    m51 = M[x_][x_]*M[px_][delta_] - M[px_][x_]*M[x_][delta_]
    m52 = M[x_][px_]*M[px_][delta_] - M[px_][px_]*M[x_][delta_]
    print(f"\nm_15 = {m51:21.14e} m_16 = {m52:21.14e}")
    print(f"m_15 = {M[ct_][x_]:21.14e} m_16 = {M[ct_][px_]:21.14e}")

    A = compute_ring_twiss(M)
    prt_np_mat("\nA:\n", "14.6e", A)

    n_dof = 2
    n     = 2*n_dof

    [nus, stable] = compute_nus(n_dof, M)
    print(f"\ncompute_nus:\n  nu    = [{nus[X_]:5.3f}, {nus[Y_]:5.3f}]")

    [nus, stable] = compute_nus_symp_mat(n_dof, M)
    print(f"\ncompute_nus_symp_mat:\n  nu    = "
          "[{nus[X_]:5.3f}, {nus[Y_]:5.3f}]")

    [eta, alpha, beta, nu, stable] = compute_twiss_M(M)
    prt_twiss("\ncompute_twiss_M:\n", eta, alpha, beta, nu)

    [eta, alpha, beta, nu] = compute_twiss_A(A)
    prt_twiss("\ncompute_twiss_A:\n", eta, alpha, beta, nu)

    [eta, alpha, beta, nu] = compute_twiss_A_A_tp(A)
    prt_twiss("\ncompute_twiss_A_A_tp:\n", eta, alpha, beta, nu)

    # Cross check.
    prt_np_mat("\nA^-1*M*A:\n", "14.6e",
               np.linalg.multi_dot([np.linalg.inv(A), M, A]))


def compute_closed_orbit(acc, conf, delta, n_max, eps):
    jj  = np.zeros(ss_dim, int)
    x0  = ss_vect_double()
    x1  = ss_vect_double()
    dx  = ss_vect_double()
    dx0 = ss_vect_tps()
    I   = ss_vect_tps()
    M   = ss_vect_tps()

    debug = True

    n_loc = len(acc)

    if conf.Cavity_on:
        n = 6
    else:
        n = 4

    first = True

    if debug:
        print("\ncompute_closed_orbit:")

    if (conf.mat_meth and (first or (delta != conf.dPparticle))):
        # Recompute transport matrices.
        if (debug):
            print(f"  recomputing transport matrices:  delta = {delta:9.3e}"
                  " ({conf.dPparticle:9.3e}) first = {first:1d")
        get_lin_maps(delta)
        conf.dPparticle = delta
        first = False

    conf.dPparticle = delta

    x0.set_zero()
    x0[delta_] = delta

    for k in range(ss_dim):
        jj[k] = 1 if k < n else 0

    if False:
        # For 2.5 D.O.F. initial COD estimate is: eta*delta.
        for k in range(2):
            x0[2*k]   = Cell[2*k].Eta[X_]*delta
            x0[2*k+1] = Cell[2*k+1].Etap[X_]*delta

    if debug:
        print(f"  {0:d}", end="")
        prt_ps_vec("                  ", x0)

    n_iter = 0
    I.set_identity()
    dx_abs = 1e30
    while (dx_abs >= eps) and (n_iter <= n_max):
        n_iter += 1
        M.set_identity()
        M += x0

#        acc.propagate(conf, M, s_loc)
        s_loc = n_loc
        acc.propagate(conf, M)

        if (True or (s_loc == n_Loc)):
            x1 = M.cst()
            dx = x0 - x1
            tmp = M-I-x1
            dx0 = partialInverse(tmp, jj) * dx
            dx_abs = xabs(n, dx)
            x0 += dx0.cst()
        else:
            dx_abs = NAN
            break

        if debug:
            print(f"{n_iter:3d} {dx_abs:7.1e} ({eps:7.1e})", end="")
            prt_np_vec("", "11.3e", x0)

    cod = dx_abs < eps

    if cod:
        conf.CODvect = x0
        conf.OneTurnMat = get_mat(M)
#        acc.propagate(conf, x0, s_loc)
        acc.propagate(conf, x0)
        if debug:
            prt_np_mat("\nPoincaré Map:\n", "14.6e", conf.OneTurnMat)
    else:
        print(f"\ncompute_closed_orbit: failed to converge after {n_iter:d}"
              "  delta = {delta:12.5e}, particle lost at element {s_loc:3d}"
              "  x_0   = {x0:13.5e}  x_k-1 = {Cell[s_loc-1].BeamPos:13.5e}"
              "  x_k   = {M.cst():13.5e}")

    return [x0, cod, s_loc]


t_dir  = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc  = accelerator_from_config(t_file)
conf = ConfigType()

if False:
    test_stuff(15)

cod, found, s_loc = compute_closed_orbit(acc, conf, 0e0, 10, 1e-10)
elem = acc.find("uq1", 1)
#print(repr(elem.getMultipoles()))
muls = elem.getMultipoles()
muls.setMultipole(Dip, 1e-3 - 1e-3j)
print("\n", repr(elem))
cod, found, s_loc = compute_closed_orbit(acc, conf, 0e0, 10, 1e-10)
exit()

t_map = compute_map(acc, conf)
M = get_mat(t_map)

# Reduce matrix from 7x7 to 6x6.
M = M[0:6, 0:6]
prt_np_mat("\nPoincaré Map:\n", "14.6e", M)

[eta, alpha, beta, nu, stable] = compute_twiss_M(M)
prt_twiss("\ncompute_twiss_M:\n", eta, alpha, beta, nu)

A = compute_ring_twiss(M)
prt_np_mat("\nA:\n", "14.6e", A)

# Cross check.
prt_np_mat("\nA^-1*M*A:\n", "14.6e",
           np.linalg.multi_dot([np.linalg.inv(A), M, A]))

compute_twiss_lat("linlat.out", acc, conf, get_map(A))

# ds = ds.drop(["elements", "tps"])
# ds.to_netcdf("twiss.nc")
# df = twiss_ds_to_df(ds)
# #print(df)
# df.to_csv("twiss.csv")


# twiss = linear_optics.compute_twiss_parameters(acc, conf)
# twiss.name = "twiss_parameters"
# md = accelerator_info(acc)
# md.attrs = dict(calc_config=calc_config)
# twiss_with_md = xr.merge([twiss, md])

# print(res)
