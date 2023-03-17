"""Functionality for computing linear optics
"""

import thor_scsi.lib as tslib
# from .phase_space_vector import ss_vect_tps2ps_jac, array2ss_vect_tps
from .courant_snyder import compute_A_CS
from .extract_info import accelerator_info
from .accelerator import instrument_with_standard_observers
from .phase_space_vector import omega_block_matrix, map2numpy
from .output import prt2txt, mat2txt, vec2txt, chop_array

import gtpsa

import math
import xarray as xr
import numpy as np
import copy as _copy
import logging


logger = logging.getLogger("thor_scsi")


X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]

[x_, px_, y_, py_, ct_, delta_] = [
    tslib.phase_space_index_internal.x,
    tslib.phase_space_index_internal.px,
    tslib.phase_space_index_internal.y,
    tslib.phase_space_index_internal.py,
    tslib.phase_space_index_internal.ct,
    tslib.phase_space_index_internal.delta,
]


def compute_map(
    acc: tslib.Accelerator,
    calc_config: tslib.ConfigType,
    *,
    delta: float = 0e0,
    t_map: gtpsa.ss_vect_tpsa = None,
    desc: gtpsa.desc = None,
    tpsa_order: int = 1
) -> gtpsa.ss_vect_tpsa:
    """Propagate an identity map through the accelerator
    """
    if t_map is None:
        assert desc is not None
        t_map = gtpsa.ss_vect_tpsa(desc, tpsa_order)
        t_map.set_identity()

    t_map[delta_] += delta
    acc.propagate(calc_config, t_map)
    return t_map


def acos2(sin, cos):
    """arcos(phi): 0 - 2 * pi.
       The sin part is used to determine the quadrant.
    """
    if abs(cos) > 1e0:
        raise Exception(f"arg = {cos:22.15e}")
    phi = np.arccos(cos)
    if sin < 0e0:
        phi = 2e0 * np.pi - phi
    return phi


def compute_nu(M):
    """

    """
    tr = M.trace()
    # Check if stable.
    if tr < 2e0:
        sin2 = M[0][1]*M[1][0]
        sgn = np.sign(sin2)
        sin = np.sqrt(np.fabs(sin2))
        return acos2(sgn*sin, tr / 2e0) / (2e0 * np.pi)
    else:
        print("\ncompute_nu: unstable")
        return np.nan


def compute_nus(n_dof, M):
    nus = np.zeros(n_dof, float)
    for k in range(n_dof):
        nus[k] = compute_nu(M[2 * k : 2 * k + 2, 2 * k : 2 * k + 2])
        if n_dof == 3:
            nus[2] = 1e0 - nus[2]
    return nus


def compute_nu_symp(n_dof, M):
    """Calculate normalised phase advance from a symplectic periodic matrix.
    """
    n = 2 * n_dof
    I = np.identity(4)
    tr = np.zeros(n_dof, float)
    for k in range(n_dof):
        tr[k] = np.trace(M[2 * k : 2 * k + 2, 2 * k : 2 * k + 2])
    M4b4 = M[0:4, 0:4]
    p1 = np.linalg.det(M4b4 - I)
    pm1 = np.linalg.det(M4b4 + I)
    po2, q = (p1 - pm1) / 16e0, (p1 + pm1) / 8e0 - 1e0
    if tr[X_] > tr[Y_]:
        sgn = 1
    else:
        sgn = -1

    radix = sgn * np.sqrt(po2 ** 2 - q)
    x, y = -po2 + radix, -po2 - radix
    nu = []
    nu.extend([acos2(M[0][1], x) / (2e0 * np.pi), acos2(M[2][3], y)
               / (2e0 * np.pi)])
    if n_dof == 3:
        nu.append(1e0 - acos2(M[4][5], tr[Z_] / 2e0) / (2e0 * np.pi))
    return np.array(nu)


def check_if_stable_1D(dim, M):
    # Dim is [X_, Y_] = [0, 1].
    return math.fabs(M[2*dim:2*dim+2, 2*dim:2*dim+2].trace()) <= 2e0


def check_if_stable_2D(M):
    return check_if_stable_1D(0, M), check_if_stable_1D(1, M)


def check_if_stable_3D(M):
    return \
        check_if_stable_1D(0, M), check_if_stable_1D(1, M), \
        check_if_stable_1D(2, M)


def check_if_stable(n_dof, M):
    if n_dof == 2:
        stable = check_if_stable_2D(M)
        stable = stable[0] and stable[1]
    elif n_dof == 3:
        stable = check_if_stable_3D(M)
        stable = stable[0] and stable[1] and stable[2]
    else:
        raise Exception("check_if_stable: undef. n_dof")
    return stable


def ind_1(k):
    """Index for linear TPSA term.
    """
    dof = 3
    ind = np.zeros(2*dof, int)
    ind[k] = 1
    return ind


def ind_2(j, k):
    """Index for quadratic TPSA term.
    """
    dof = 3
    ind = np.zeros(2*dof, int)
    ind[j] = 1
    ind[k] = 1
    return ind


def acos2_tpsa(sin, cos):
    """arcos(phi): 0 - 2 * pi.
       The sin part is used to determine the quadrant.
    """
    phi = gtpsa.acos(cos)
    if sin < 0e0:
        phi = 2e0 * np.pi - phi
    return phi


def compute_nu_xi(desc, tpsa_order, M):
    """Compute tune & linear chromaticity from trace of Poincaré map:
          nu + xi * delta = arccos( Trace{M} / 2 ) / ( 2 * pi ).
    """
    stable = check_if_stable(2, M.jacobian())
    if stable:
        nu, xi = [np.zeros(2), np.zeros(2)]
        M_delta = gtpsa.tpsa(desc, tpsa_order)
        for k in range(2):
            M_delta.clear()
            # m_11 + delta * m_16.
            M_delta += M[2*k].get(ind_1(2*k))
            M_delta.set(ind_1(delta_), 0e0, M[2*k].get(ind_2(2*k, delta_)))
            # m_22 + delta * m_26.
            M_delta += M[2*k+1].get(ind_1(2*k+1))
            M_delta.set(ind_1(delta_), 1e0, M[2*k+1].get(ind_2(2*k+1, delta_)))
            # M_delta = m_11 + m_22.
            nu_tpsa = \
                acos2_tpsa(M.jacobian()[2*k][2*k+1], M_delta/2e0)/(2e0*np.pi)
            nu[k], xi[k] = [nu_tpsa.get(), nu_tpsa.get(ind_1(delta_))]
    else:
        nu = np.nan
        xi = np.nan
        print("\ncompute_nu_xi: unstable ", stable)

    return stable, nu, xi


def find_closest_nu(nu, w):
    min = 1e30
    for k in range(w.size):
        nu_k = acos2(w[k].imag, w[k].real) / (2e0 * np.pi)
        diff = np.abs(nu_k - nu)
        if diff < min:
            [ind, min] = [k, diff]
    return ind


def sort_eigen_vec(n_dof, nu, w):
    """
    Args:
        nu:    eigenvalues / computed tunes
        n_dof: degrees of freedom
        w:     complex eigen values (from eigen vector computation
               procedure)

    Todo:
        Check that vectors are long enough ...
    """
    order = np.zeros(2 * n_dof, int)
    for k in range(n_dof):
        order[2 * k] = find_closest_nu(nu[k], w)
        order[2 * k + 1] = find_closest_nu(1e0 - nu[k], w)
    return order


def compute_A(n_dof, eta, u):
    """

    Original version thanks to Johan
    Diagonalize the Poincaré map:
      M = B * A * R * A^-1 * B^-1

    Todo:
       what's is special about this inverse
    """
    sign = np.sign

    n = 2 * n_dof

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
            sgn_vec * (u1[:, 2 * i].real + sgn_im * u1[:, 2 * i].imag * 1j)
            / scl,
            sgn_vec
            * (u1[:, 2 * i + 1].real + sgn_im * u1[:, 2 * i + 1].imag * 1j)
            / scl,
        ]

    for i in range(n_dof):
        [A[:n, 2 * i], A[:n, 2 * i + 1]] = \
            [u1[:, 2 * i].real, u1[:, 2 * i].imag]

    if n_dof == 2:
        # For coasting beam translate to momentum dependent fix point.

        B = np.identity(6)
        B[x_, delta_], B[px_, delta_] = eta[x_], eta[px_]
        B[ct_, x_], B[ct_, px_] = eta[px_], -eta[x_]
        logger.info("\ncompute_A\nB:\n" + mat2txt(B))

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


def compute_M_diag(
        n_dof: int, M: np.ndarray
) -> [bool, np.ndarray, np.ndarray, np.ndarray]:
    """

    Args:
        M:     transfer matrix
        n_dof: degrees of freedom

    Return:

    See xxx reference
    """
    stable = check_if_stable(n_dof, M)
    if stable:
        n = 2 * n_dof

        nu_symp = compute_nu_symp(n_dof, M)
        logger.debug("computed tunes (for symplectic matrix): %s", nu_symp)

        # Diagonalise M.
        
        [w, u] = np.linalg.eig(M[:n, :n])

        # nu_eig = acos2(w.imag, w.real) / (2e0 * np.pi)
        nu_eig = np.zeros(n)
        for k in range(n):
            if abs(w[k].real) > 1e0:
                print("\ncompute_M_diag: |arg| for acos2 > 1 ", w[k])
                return False, np.nan, np.nan, np.nan

            nu_eig[k] = acos2(w[k].imag, w[k].real) / (2e0 * np.pi)

        logger.debug(
            "\nu:\n"
            + mat2txt(u)
            + "\nnu_symp:\n"
            + vec2txt(nu_symp)
            + "\nnu_eig:\n"
            + vec2txt(nu_eig)
            + "\nlambda:\n"
            + vec2txt(w)
        )

        order = sort_eigen_vec(n_dof, nu_symp, w)

        w_ord = np.zeros(n, complex)
        u_ord = np.zeros((n, n), complex)
        nu_eig_ord = np.zeros(n, float)
        for k in range(n):
            w_ord[k] = w[order[k]]
            u_ord[:, k] = u[:, order[k]]
            nu_eig_ord[k] = acos2(w_ord[k].imag, w_ord[k].real) / (2e0 * np.pi)

        logger.debug(
            "\norder:\n"
            + vec2txt(order)
            + "\nnu_eig_ord:\n"
            + vec2txt(nu_eig_ord)
            + "\nlambda_ord:\n"
            + vec2txt(w_ord)
            + "\nu_ord^-1.M.u_ord:\n"
            + mat2txt(chop_array(np.linalg.inv(u_ord) @ M[:n, :n] @ u_ord,
                                 1e-15))
        )

        eta = compute_dispersion(M)

        logger.debug("\neta:\n" + vec2txt(eta))

        [A, A_inv, u1] = compute_A(n_dof, eta, u_ord)
        A, _ = compute_A_CS(2, A)

        logger.debug(
            "\nu1:\n"
            + mat2txt(u1)
            + "\nu1^-1.M.u1:\n"
            + mat2txt(chop_array(np.linalg.inv(u1) @ M[:n, :n] @ u1, 1e-13))
            + "\nu1^T.omega.u1:\n"
            + mat2txt(chop_array(u1.T @ omega_block_matrix(n_dof) @ u1, 1e-13))
            + "\nA:\n"
            + mat2txt(chop_array(A, 1e-13))
            + "\nA^T.omega.A:\n"
            + mat2txt(chop_array(A[:n, :n].T @ omega_block_matrix(n_dof)
                                 @ A[:n, :n], 1e-13))
        )

        R = A_inv @ M @ A

        nu = np.zeros(3, float)
        alpha_rad = np.zeros(3, float)
        for k in range(n_dof):
            nu[k] = \
                np.arctan2(R[2 * k][2 * k + 1], R[2 * k][2 * k]) / (2e0 * np.pi)
            if (nu[k] < 0e0) and (k < 2):
                nu[k] += 1e0
            if n_dof == 3:
                alpha_rad[k] = \
                    np.log(np.sqrt(np.absolute(w_ord[2 * k])
                                   * np.absolute(w_ord[2 * k + 1])))

        A_CS, dnu_cs = compute_A_CS(n_dof, A)
        logger.debug(
            "\n\nA_CS:\n"
            + mat2txt(chop_array(A_CS, 1e-10))
            + "\n\nR:\n"
            + mat2txt(chop_array(R, 1e-10))
            + "\n\nnu        = {:18.16f} {:18.16f} {:18.16f}".
            format(nu[X_], nu[Y_], nu[Z_])
        )
        logger.debug(f"dnu_cs  {dnu_cs}")
        if n_dof == 3:
            logger.debug(
                "alpha_rad = {:13.6e} {:13.6e} {:13.6e}".format(
                    alpha_rad[X_], alpha_rad[Y_], alpha_rad[Z_]
                )
            )
    else:
        A         = np.nan
        A_inv     = np.nan
        alpha_rad = np.nan

    return stable, A, A_inv, alpha_rad


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
    # not accidentally zero
    if dnu < -eps:
        dnu += 1.0

    return alpha, beta, dnu


def _extract_tps(elem: tslib.Observer) -> tslib.ss_vect_tps:
    """Extract tps data from the elments' observer
    """
    ob = elem.get_observer()
    if not ob.has_truncated_power_series_a():
        msg = f"Observer {ob} has no  TPSA"
        logger.error(msg)
        raise AssertionError(msg)
    else:
        logger.debug("Data available for eleeent index %5d name %s",
                     ob.get_observed_index(), ob.get_observed_name())
    return ob.get_truncated_power_series_a()


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


def tps2twiss(tpsa: gtpsa.ss_vect_tpsa) -> (np.ndarray, np.ndarray):
    """Extract twiss information from truncated power series

    Warning:
        It is thue users responsibility to ensure that the matrix
        is a matrix rotationg towards to Floquet space
    """
    Aj = tpsa.jacobian()
    r = transform_matrix_extract_twiss(Aj)
    return r


def compute_map_and_diag(
        n_dof: int,
        acc: tslib.Accelerator,
        calc_config: tslib.ConfigType,
        *,
        desc: gtpsa.desc,
        tpsa_order: int = 1):
    """propagate once around ring. use this map to find phase space origin

    Todo:
         Rename phase space origin fix point
    """
    t_map = compute_map(acc, calc_config, desc=desc, tpsa_order=tpsa_order)
    M = np.array(t_map.jacobian())
    logger.info("\ncompute_map_and_diag\nM:\n" + mat2txt(M))

    stable, A, A_inv, alpha_rad = compute_M_diag(n_dof, M)
    if stable:
        Atest = gtpsa.ss_vect_tpsa(desc, 1)
        Atest.set_zero()
        Atest.set_jacobian(A)
        acc.propagate(calc_config, Atest)
    else:
        print("\ncompute_map_and_diag: unstable")
    return stable, t_map, A


def compute_Twiss_along_lattice(
    n_dof: int,
    acc: tslib.Accelerator,
    calc_config: tslib.ConfigType = None,
    *,
    A: gtpsa.ss_vect_tpsa = None,
    desc: gtpsa.desc = None,
    tpsa_order: int = 1
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
    if calc_config is None:
        calc_config = tslib.ConfigType()

    if A is None:
        _, A = \
            compute_map_and_diag(n_dof, acc, calc_config, desc=desc,
                                 tpsa_order=tpsa_order)
    logger.info("\ncompute_Twiss_along_lattice\nA:\n" + mat2txt(A))

    # Not really required ... but used for convenience
    observers = instrument_with_standard_observers(acc)

    # Propagate through the accelerator
    A_map = gtpsa.ss_vect_tpsa(desc, 1)
    A_map.set_zero()
    A_map.set_jacobian(A)
    logger.debug("\ncompute_Twiss_along_lattice\nA:\n" + prt2txt(A_map))

    for k in range(len(acc)):
        acc.propagate(calc_config, A_map, k, 1)
        # Zero the phase advance so that the fraction tune change is not
        # exceeding two pi
        Aj = A_map.jacobian()
        rjac, _ = compute_A_CS(2, Aj)
        A_map.set_jacobian(rjac)

    logger.debug("\ncompute_Twiss_along_lattice A:\n%s" + prt2txt(A_map))

    indices = [elem.index for elem in acc]
    tps_tmp = [_extract_tps(elem) for elem in acc]
    data = [tps2twiss(t) for t in tps_tmp]
    # print(type(tps_tmp), type(tps_tmp[0]))
    tps_tmp = np.array(tps_tmp, dtype=object)
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
    res = \
        info.merge(dict(twiss=twiss_parameters, dispersion=dispersion, tps=tps))
    return res


__all__ = [
    "compute_map", "compute_nu_symp", "check_if_stable_1D",
    "check_if_stable_2D", "check_if_stable_3D", "check_if_stable_3D",
    "compute_nu_xi", "compute_map_and_diag", "compute_Twiss_along_lattice",
    "jac2twiss", "compute_M_diag"
]
