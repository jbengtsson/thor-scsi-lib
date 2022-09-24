"""
"""
import thor_scsi.lib as tslib
from .accelerator import instrument_with_radiators
from .closed_orbit import compute_closed_orbit
from . import linalg
from . import linear_optics as lo
from .output import vec2txt, mat2txt
from .phase_advance import compute_nus_for_symplectic_matrix
from .phase_space_vector import ss_vect_tps2ps_jac, array2ss_vect_tps
import numpy as np
from dataclasses import dataclass

import logging
logger = logging.getLogger("thor-scsi")

X_ = 0
Y_ = 1
Z_ = 2

@dataclass
class RadiationResult:
    #: the different taus
    relaxation_constants: np.ndarray
    emittance: np.ndarray
    #: for consistency checks
    fractional_tunes: np.ndarray


def acos2(sin, cos):
    # Calculate the normalised phase advance from the trace = 2*2*pi* nu of
    # the Poincaré map; i.e., assuming mid-plane symmetry.
    # The sin part is used to determine the quadrant.
    mu = np.arccos(cos)
    if sin < 0e0:
        mu = 2e0*np.pi - mu
    return mu


def calculate_nu(M):
    tr = M.trace()
    # Check if stable.
    if tr < 2e0:
        calculate_nu(tr/2e0,  M[0][1])/(2e0*np.pi)
        return nu
    else:
        print("\ncalculate_nu: unstable\n")
        return float('nan')


def calculate_nus(n_dof, M):
    nus = np.zeros(n_dof, float)
    for k in range(n_dof):
        nus[k] = calculate_nu(M[2*k:2*k+2, 2*k:2*k+2])/(2e0*np.pi)
        if n_dof == 3:
            nus[2] = 1e0 - nus[2]
    return nus


def calculate_nu_symp(M):
    # Calculate normalised phase advance from a symplectic periodic matrix.
    dof = 2
    n = 2*dof
    I = np.identity(n)
    tr = np.zeros(3, float)
    for k in range(3):
        tr[k] = np.trace(M[2*k:2*k+2, 2*k:2*k+2])
    M4b4 = M[0:n, 0:n]
    [p1, pm1] = [np.linalg.det(M4b4-I), np.linalg.det(M4b4+I)]
    [po2, q] = [(p1-pm1)/16e0, (p1+pm1)/8e0 - 1e0]
    if tr[X_] > tr[Y_]:
        sgn = 1
    else:
        sgn = -1
    [x, y] = [-po2+sgn*np.sqrt(po2**2-q), -po2-sgn*np.sqrt(po2**2-q)]
    nu = \
        [acos2(M[0][1], x)/(2e0*np.pi), acos2(M[2][3], y)/(2e0*np.pi),
         1e0-acos2(M[4][5], tr[Z_]/2e0)/(2e0*np.pi)]
    return np.array(nu)


def find_closest_nu(nu, w):
    min = 1e30
    for k in range(w.size):
        nu_k = acos2(w[k].imag, w[k].real)/(2e0*np.pi)
        diff =  np.abs(nu_k-nu)
        if diff < min:
            [ind, min] = [k, diff]
    return ind


def sort_eigen_vec(nu, w):
    order = []
    for k in range(nu.size):
        order.append(find_closest_nu(nu[k], w))
    for k in range(nu.size):
        order.append(find_closest_nu(1e0-nu[k], w))
    return np.array(order)


def calculate_radiation(
    acc: tslib.Accelerator,
    *,
    energy,
    calc_config: tslib.ConfigType = None,
    install_radiators: bool = True,
    dof=3
):
    """

    Todo:
        Rename function
        Inspect if radiators are installed?
        check in calc config if radation is requested
        in this case install the radiators

    1. calculate fix point and Poincarè Map M with damped system
       (i.e. radiation on and cavity on (without dispersion in a second case)
    2. diagonalise M = A $\Gamma$ A$^{-1}$
    3. eigenvalues:

        - complex part: tunes,
        - real part: damping times  (refer equation)

       use eigen values of symplectic matrix to identify the planes

    4. propagate A, thin kick will create diffusion coeffs (don't forget to
       zero them before calculation starts (sum it up afterwards

    """

    if install_radiators:
        logger.debug("Installing radiators")
        # keep variable as long as you need to do calculations
        rad_del = instrument_with_radiators(acc, energy=energy)

    if calc_config is None:
        raise AssertionError
        calc_config = tslib.ConfigType()
        calc_config.radiation = True
        # is this used anywhere?
        calc_config.emittance = False
        calc_config.Cavity_on = True

    calc_config.Energy = energy
    logger.debug(
        f"calc_config radiation { calc_config.radiation} emmittance {calc_config.emittance} Cavity on {calc_config.Cavity_on}"
    )

    # Compute the fixed point ... radation is on
    r = compute_closed_orbit(acc, calc_config, delta=0e0)
    print("r.one_turn_map")
    print(mat2txt(r.one_turn_map))

    nu = calculate_nu_symp(r.one_turn_map)
    print("\nnu = [{:7.5f}, {:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_], nu[Z_]))

    # diagonalise M
    n = 2 * dof
    M_tp = np.transpose(r.one_turn_map[:n, :n])

    # Finds eigen values and eigen vectors ... same functionallity as in
    # find_phase_space_origin / find_phase_space_fix_point
    w, v = np.linalg.eig(M_tp)
    print("\nlambda:")
    for w_k in w:
        nu_k = acos2(w_k.imag, w_k.real)/(2e0*np.pi)
        print(" %7.5f" % (nu_k), end="")
    print()

    order = sort_eigen_vec(nu, w)
    print("\norder:\n", order)
    w_ord= np.array([w[order[0]], w[order[3]], w[order[1]], w[order[4]],
                     w[order[2]], w[order[5]]])
    v_ord= np.array([v[order[0]], v[order[3]], v[order[1]], v[order[4]],
                     v[order[2]], v[order[5]]])
    print("\nw:")
    print(vec2txt(w))
    print("\nw_ord:")
    print(vec2txt(w_ord))
    print("\nv:")
    print(mat2txt(v))
    print("\nv_ord:")
    print(mat2txt(v_ord))

    eta = np.zeros(6, np.float)
    A_inv, v1 = linalg.compute_A_inv_prev(dof, eta, v_ord)

    A = np.linalg.inv(A_inv)
    print("\nA:")
    print(mat2txt(A))
    return

    Atmp = np.zeros([7, 7], dtype=np.float)
    Atmp[:6, :6] = A
    print("Atmp ")
    print(mat2txt(Atmp))
    # Atmp[[2, 3, 4, 5], :] = Atmp[[4, 5, 2, 3], :]
    print("Atmp resuffled")
    print(mat2txt(Atmp))
    Ap = array2ss_vect_tps(Atmp)
    return

    print("Ap before calc")
    print(Ap)

    # Diffusion coefficients
    acc.propagate(calc_config, Ap)
    print("Ap after calc")
    print(Ap)

    r = RadiationResult(relaxation_constants=w.real, fractional_tunes=w.imag)
