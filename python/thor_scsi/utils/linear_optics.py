"""Functionallity for computing linear optics
"""
import thor_scsi.lib as tslib
from .phase_space_vector import ss_vect_tps2ps_jac, array2ss_vect_tps
from .courant_snyder import compute_A_CS
from .extract_info import accelerator_info
from .accelerator import instrument_with_standard_observers
from .phase_space_vector import omega_block_matrix, map2numpy
from .output import mat2txt, vec2txt, chop_array

import xarray as xr
import numpy as np
import copy as _copy
import logging


logger = logging.getLogger("thor_scsi")


[X_, Y_, Z_] = [0, 1, 2]

[x_, px_, y_, py_, ct_, delta_] = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta,
]


def acos2(sin, cos):
    # Calculate the normalised phase advance from the Poincaré map:
    #   Tr{M} = 2 cos(2 pi nu)
    # i.e., assuming mid-plane symmetry.
    # The sin part is used to determine the quadrant.
    mu = np.arccos(cos)
    if sin < 0e0:
        mu = 2e0 * np.pi - mu
    return mu


def calculate_nu(M):
    tr = M.trace()
    # Check if stable.
    if tr < 2e0:
        calculate_nu(tr / 2e0, M[0][1]) / (2e0 * np.pi)
        return nu
    else:
        print("\ncalculate_nu: unstable\n")
        return float("nan")


def calculate_nus(n_dof, M):
    nus = np.zeros(n_dof, float)
    for k in range(n_dof):
        nus[k] = calculate_nu(M[2 * k : 2 * k + 2, 2 * k : 2 * k + 2]) / (2e0 * np.pi)
        if n_dof == 3:
            nus[2] = 1e0 - nus[2]
    return nus


def calculate_nu_symp(n_dof, M):
    """Calculate normalised phase advance from a symplectic periodic matrix.
    """
    n = 2 * n_dof
    I = np.identity(4)
    tr = np.zeros(n_dof, float)
    for k in range(n_dof):
        tr[k] = np.trace(M[2 * k : 2 * k + 2, 2 * k : 2 * k + 2])
    M4b4 = M[0:4, 0:4]
    p1  = np.linalg.det(M4b4 - I)
    pm1 = np.linalg.det(M4b4 + I)
    po2, q = (p1 - pm1) / 16e0, (p1 + pm1) / 8e0 - 1e0
    if tr[X_] > tr[Y_]:
        sgn = 1
    else:
        sgn = -1

    radix = sgn * np.sqrt(po2 ** 2 - q)
    x, y = -po2 + radix, -po2 - radix
    nu = []
    nu.extend([
        acos2(M[0][1], x) / (2e0 * np.pi),
        acos2(M[2][3], y) / (2e0 * np.pi)
    ])
    if n_dof == 3:
        nu.append(1e0 - acos2(M[4][5], tr[Z_] / 2e0) / (2e0 * np.pi))
    return np.array(nu)


def find_closest_nu(nu, w):
    min = 1e30
    for k in range(w.size):
        nu_k = acos2(w[k].imag, w[k].real) / (2e0 * np.pi)
        diff = np.abs(nu_k - nu)
        if diff < min:
            [ind, min] = [k, diff]
    return ind


def sort_eigen_vec(dof, nu, w):
    """
    Args:
        nu : eigenvalues / computed tunes
        dof: degrees of freedom
        w :  complex eigen values (from eigen vector computation
             procedure)

    Todo:
        Check that vectors are long enough ...
    """
    order = np.zeros(2 * dof, int)
    for k in range(dof):
        order[2 * k] = find_closest_nu(nu[k], w)
        order[2 * k + 1] = find_closest_nu(1e0 - nu[k], w)
    return order


def compute_A(n_dof, eta, u):
    """

    Original version thanks to Johan

    Todo:
       what's is special about this inverse
    """
    sign = np.sign

    n = 2 * n_dof

    # Should be a local copy.
    # u1 = np.array(u)
    u1 = _copy.copy(u)
    A = np.identity(6)
    S = omega_block_matrix(n_dof)

    # Normalise eigenvectors: A^T.omega.A = omega.
    for i in range(n_dof):
        z = u1[:, 2 * i].real @ S @ u1[:, 2 * i].imag
        sgn_im = sign(z)
        scl = np.sqrt(np.abs(z))
        sgn_vec = sign(u1[2 * i][2 * i].real)
        [u1[:, 2 * i], u1[:, 2 * i + 1]] = [
            sgn_vec * (u1[:, 2 * i].real + sgn_im * u1[:, 2 * i].imag * 1j) / scl,
            sgn_vec
            * (u1[:, 2 * i + 1].real + sgn_im * u1[:, 2 * i + 1].imag * 1j)
            / scl,
        ]

    for i in range(n_dof):
        [A[:n, 2 * i], A[:n, 2 * i + 1]] = [u1[:, 2 * i].real, u1[:, 2 * i].imag]

    if n_dof == 2:
        # For coasting beam translate to momentum dependent fix point.

        B = np.identity(6)
        B[x_, delta_], B[px_, delta_] = eta[x_], eta[px_]
        B[ct_, x_], B[ct_, px_] = eta[px_], -eta[x_]

        A = B @ A

    A_inv = np.linalg.inv(A)

    return A, A_inv, u1


def compute_dispersion(M: np.ndarray) -> np.ndarray:
    """

    Todo:
        Define which standard M should follow standard rules or
        thor_scsi / tracy internal rules
    """
    n = 4
    id_mat = np.identity(n)
    D = M[:n, tslib.phase_space_index_internal.delta]

    return np.linalg.inv(id_mat - M[:n, :n]) @ D


#: scale arctan2 (a12/a11) to Floquet coordinates (correct?)
scale_nu = 1.0 / (2 * np.pi)


def jac2twiss(A: np.ndarray) -> (float, float, float):
    """extract twiss parameters from array for one plane

    returns: alpha, beta, dnu
    See XXX eq (172)
    """
    eps = 1e-15

    A = np.atleast_2d(A)
    nr, nc = A.shape
    if nr != 2:
        raise AssertionError(f"expected 2 rows got {nr}")

    if nc != 2:
        raise AssertionError(f"expected 2 cols got {nc}")

    del nr, nc

    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1, 1]

    alpha = -(a11 * a21 + a12 * a22)
    beta = a11 ** 2 + a12 ** 2

    dnu = np.arctan2(a12, a11) * scale_nu
    # epsilon must be negative here ... ensure that the value is
    # not accidentially zero
    if dnu < -eps:
        dnu += 1.0

    return alpha, beta, dnu


def _extract_tps(elem: tslib.Observer) -> tslib.ss_vect_tps:
    """Extract tps data from the elments' observer
    """
    ob = elem.getObserver()
    assert ob.hasTruncatedPowerSeries()
    return ob.getTruncatedPowerSeries()


def transform_matrix_extract_twiss(A: np.ndarray) -> (np.ndarray, np.ndarray):
    """Extract twiss parameters from matrix rotating towards Floquet space

    Returns:
        eta, alpha, beta, nu

    """
    # njac, dnus = compute_A_CS(2, jac)
    tmp = [jac2twiss(A[p : p + 2, p : p + 2]) for p in range(0, 4, 2)]
    twiss = np.array(tmp)

    # only sensible for the first two planes ...
    delta = tslib.phase_space_index_internal.delta
    etas = A[:4, delta]

    return etas, twiss


def tps2twiss(tps: tslib.ss_vect_tps) -> (np.ndarray, np.ndarray):
    """Extract twiss information from turncated power series

    Warning:
        It is thue users responsibility to ensure that the matrix
        is a matrix rotationg towards to Floquet space
    """
    ps, Aj = ss_vect_tps2ps_jac(tps)
    return transform_matrix_extract_twiss(Aj)


def find_phase_space_origin(M: np.ndarray) -> np.ndarray:
    """Transform to energy dependent fix point

    * Diagonalise M = A R A^(-1)

    This gives:
    * Transfrom to the energy dependent fix point
    * Diagonalising in [x, px, y, py]

       R : Block diagonal
       :math:`A^(-1)` : From phase space ellipse to Floquet space circle

    * Transfrom from energy dependent fix point math::`\eta` math::`\delta` to
      origin of phase space


    Todo:
        Revisit naming: fix point instead of origin
        Add reference to Johan's Tech Note

    """

    n_dof = 2
    n = 2 * n_dof

    M_tp = np.transpose(M[:n, :n])

    # Finds eigen values and eigen vectors
    w, v = np.linalg.eig(M_tp)

    # Sort eigen values to match the trace of the plane
    # print(f"before sort w:{mat2txt(w)}  ")
    # print(f"before sort v:{mat2txt(v)}  ")
    w, v = match_eigenvalues_to_plane_orig(M_tp, w, v, n_dof=n_dof)
    # wabs = np.absolute(w)
    # print(f"sort  wabs   {mat2txt(wabs)}")
    # print(f"after sort w:{mat2txt(w)}  ")
    # print(f"after sort v:{mat2txt(v)}  ")

    eta = compute_dispersion(M)
    A_inv, v1 = compute_A_inv_with_dispersion(v, eta, n_dof=n_dof)
    A = np.linalg.inv(A_inv)
    A, _ = compute_A_CS(2, A)

    return A


def compute_map(
    acc: tslib.Accelerator, calc_config: tslib.ConfigType
) -> tslib.ss_vect_tps:
    """Propagate an identity map through the accelerator
    """
    t_map = tslib.ss_vect_tps()
    t_map.set_identity()
    # acc.propagate(calc_config, t_map, 0, len(acc))
    acc.propagate(calc_config, t_map)
    return t_map


def propagate_and_find_phase_space_orgin(acc, calc_config):
    """propagate once around ring. use this map to find phase space origin

    Todo:
         Rename phase space origin fix point
    """
    t_map = compute_map(acc, calc_config)
    M = map2numpy(t_map)

    Aj = find_phase_space_origin(M)
    Atmp = np.zeros([7, 7], dtype=np.float)
    Atmp[:6, :6] = Aj
    A = array2ss_vect_tps(Atmp)

    return A


def compute_twiss_along_lattice(
    acc: tslib.Accelerator,
    calc_config: tslib.ConfigType = None,
    A: tslib.ss_vect_tps = None,
) -> xr.Dataset:
    """

    Args:
        acc :         an :class:`tlib.Accelerator` instance
        calc_config : an :class:`tlib.Config` instance
        A : see :func:`compute_ring_twiss` for output

    returns xr.Dataset

    Todo:
        Consider if tps should be returned too
        Split up functionallity (not everyone wants output as
        xarrays. But what then to use?)

    """
    if not calc_config:
        calc_config = tslib.ConfigType()

    if A is None:
        A = propagate_and_find_phase_space_orgin(acc, calc_config)

    # Not really required ... but used for convienience
    instrument_with_standard_observers(acc)

    # Propagate through the accelerator
    for k in range(len(acc)):
        acc.propagate(calc_config, A, k, k + 1)
        # Extract the first order term from the map representation of A
        # thus called Aj as it is similar to the Jacobian
        _, Aj = ss_vect_tps2ps_jac(A)
        # A will be rotated ... so rotate it back so that the
        # next step starts at a Courant Snyder form
        rjac, _ = compute_A_CS(2, Aj)
        Atmp = np.zeros([7, 7], dtype=np.float)
        Atmp[:6, :6] = rjac
        A = array2ss_vect_tps(Atmp)

    indices = [elem.index for elem in acc]
    tps_tmp = [_extract_tps(elem) for elem in acc]
    data = [tps2twiss(t) for t in tps_tmp]
    tps_tmp = np.array(tps_tmp)
    twiss_pars = [d[1] for d in data]
    disp = [d[0] for d in data]

    # Stuff information into xarrays ..
    dispersion = xr.DataArray(
        data=disp,
        name="dispersion",
        dims=["index", "phase_coordinate"],
        coords=[indices, ["x", "px", "y", "py"]],
    )
    twiss_parameters = xr.DataArray(
        twiss_pars,
        name="twiss",
        dims=["index", "plane", "par"],
        coords=[indices, ["x", "y"], ["alpha", "beta", "dnu"]],
    )
    phase_space_coords_names = ["x", "px", "y", "py", "ct", "delta"]
    tps = xr.DataArray(
        data=tps_tmp,
        name="tps",
        dims=["index", "phase_coordinate"],
        coords=[indices, phase_space_coords_names],
    )
    info = accelerator_info(acc)
    res = info.merge(dict(twiss=twiss_parameters, dispersion=dispersion, tps=tps))
    return res


def compute_M_diag(
    dof: int, M: np.ndarray, debug_prt: bool = False
) -> [np.ndarray, np.ndarray]:
    """

    Args:
        M: transfer matrix
        dof: degrees of freedom


    Return:

    See xxx reference
    """
    n = 2 * dof

    nu_symp = calculate_nu_symp(dof, M)
    # Diagonalise M.
    [w, u] = np.linalg.eig(M[:n, :n])

    # nu_eig = acos2(w.imag, w.real) / (2e0 * np.pi)
    nu_eig = np.zeros(n)
    for k in range(n):
        nu_eig[k] = acos2(w[k].imag, w[k].real) / (2e0 * np.pi)

    if debug_prt:
        print("\nu:\n" + mat2txt(u))
        print("\nnu_symp:\n" + vec2txt(nu_symp))
        print("\nnu_eig:\n" + vec2txt(nu_eig))
        print("\nlambda:\n" + vec2txt(w))
        # print("\nu:\n"+mat2txt(u))

    order = sort_eigen_vec(dof, nu_symp, w)

    w_ord = np.zeros(n, complex)
    u_ord = np.zeros((n, n), complex)
    nu_eig_ord = np.zeros(n, float)
    for k in range(n):
        w_ord[k] = w[order[k]]
        u_ord[:, k] = u[:, order[k]]
        nu_eig_ord[k] = acos2(w_ord[k].imag, w_ord[k].real) / (2e0 * np.pi)

    if debug_prt:
        print("\norder:\n", order)
        print("\nnu_eig_ord:\n" + vec2txt(nu_eig_ord))
        print("\nlambda_ord:\n" + vec2txt(w_ord))
        # print("\nu_ord:\n"+mat2txt(u_ord))
        print(
            "\nu_ord^-1.M.u_ord:\n"
            + mat2txt(chop_array(np.linalg.inv(u_ord) @ M[:n, :n] @ u_ord, 1e-15))
        )

    eta = compute_dispersion(M)

    if debug_prt:
        print("\neta:\n", eta)

    [A, A_inv, u1] = compute_A(dof, eta, u_ord)

    if debug_prt:
        print("\nu1:\n" + mat2txt(u1))
        print(
            "\nu1^-1.M.u1:\n"
            + mat2txt(chop_array(np.linalg.inv(u1) @ M[:n, :n] @ u1, 1e-13))
        )
        print(
            "\nu1^T.omega.u1:\n"
            + mat2txt(chop_array(u1.T @ omega_block_matrix(dof) @ u1, 1e-13))
        )
        print("\nA:\n" + mat2txt(chop_mat(A_inv, 1e-13)))
        print(
            "\nA^T.omega.A:\n"
            + mat2txt(
                chop_array(A[:n, :n].T @ omega_block_matrix(dof) @ A[:n, :n], 1e-13)
            )
        )

    print("\nA_CS:\n" + mat2txt(chop_array(compute_A_CS(dof, A)[0], 1e-10)))
    print("\nR:\n" + mat2txt(chop_array(A_inv @ M @ A, 1e-10)))

    return A, A_inv


__all__ = ["compute_twiss_along_lattice", "jac2twiss", "compute_M_diag"]
