"""Read lattice file and calculate radiation
"""
import logging
import xarray as xr

#logging.basicConfig(level="DEBUG")

from thor_scsi.factory import accelerator_from_config
from thor_scsi.lib import (
    ConfigType,
    ss_vect_tps,
    ss_vect_double,
    RadiationDelegate,
    RadiationDelegateKick,
    ObservedState
)
from thor_scsi.lib import phase_space_index_internal as phase_space_ind

logging.basicConfig(level="DEBUG")
from thor_scsi.factory import accelerator_from_config
from thor_scsi.utils.accelerator import instrument_with_radiators
from thor_scsi.utils.radiate import calculate_radiation
import os
import numpy as np

import thor_scsi.lib as tslib

from thor_scsi.utils import linalg
from thor_scsi.utils.closed_orbit import compute_closed_orbit
from thor_scsi.utils import linear_optics as lo

from thor_scsi.utils.output import vec2txt, mat2txt


X_ = 0
Y_ = 1
Z_ = 2


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


def sort_eigen_vec(dof, nu, w):
    order = []
    for k in range(dof):
        order.append(find_closest_nu(nu[k], w))
    for k in range(dof):
        order.append(find_closest_nu(1e0-nu[k], w))
    return np.array(order)


t_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
#t_file = os.path.join(t_dir, "b3_tst.lat")
t_file = os.path.join(t_dir, "b3_sf_40Grad_JB.lat")

acc = accelerator_from_config(t_file)
print(" ".join([elem.name for elem in acc]))
print("Length", np.sum([elem.getLength() for elem in acc]))

# b2 = acc.find("b2", 0)

# energy = 2.5e9

# cav.setVoltage(cav.getVoltage() * 1./2.)
# cav.setVoltage(0)
cav = acc.find("cav", 0)
print("acc cavity", repr(cav))
txt=\
    f"""Cavity info
frequency         {cav.getFrequency()/1e6} MHz",
voltage           {cav.getVoltage()/1e6} MV
harmonic number   {cav.getHarmonicNumber()}
    """
print(txt)

# cav.setVoltage(cav.getVoltage() * 1./2.)
# cav.setVoltage(0)
cav = acc.find("cav", 0)
print("acc cavity", repr(cav))
txt=\
    f"""Cavity info
frequency         {cav.getFrequency()/1e6} MHz",
voltage           {cav.getVoltage()/1e6} MV
harmonic number   {cav.getHarmonicNumber()}
    """
print(txt)

radiate = not True
calc_config = tslib.ConfigType()
calc_config.radiation = radiate
# is this used anywhere?
calc_config.emittance = False
calc_config.Cavity_on = not True

print(
    "calc_config",
    calc_config.radiation,
    calc_config.emittance,
    calc_config.Cavity_on,
)

calc_config.Energy = 2.5e9

dof = 2
n = 2 * dof

r = compute_closed_orbit(acc, calc_config, delta=0e0)
M = r.one_turn_map[:6, :6]
print("\nM:")
print(mat2txt(M))

# Diagonalise M.
M_tp = M.T[:n, :n]
[w, v] = np.linalg.eig(M_tp)
print("\nlambda:")
for w_k in w:
    nu_k = acos2(w_k.imag, w_k.real)/(2e0*np.pi)
    print(" %7.5f" % (nu_k), end="")
print()

nu = calculate_nu_symp(r.one_turn_map)
print("\nnu = [{:7.5f}, {:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_], nu[Z_]))

nu[Z_] = 1e0 - nu[Z_]
order = sort_eigen_vec(dof, nu, w)
print("\norder:\n", order)
print("\nw:")
print(vec2txt(w))
print("\nv:")
print(mat2txt(v))
w_ord = np.zeros(n, complex)
v_ord = np.zeros((n, n), complex)
for k in range(dof):
    [w_ord[2*k], w_ord[2*k+1]] = np.array([w[order[k]], w[order[k+dof]]])
    [v_ord[:, 2*k], v_ord[:, 2*k+1]] = [v[:, order[k]], v[:, order[k+dof]]]
print("\nw:")
print(vec2txt(w_ord))
print("\nv:")
print(mat2txt(v_ord))

eta = lo.compute_dispersion(M)
print("\neta:\n", eta)

if True:
    [A, v1] = linalg.compute_A_inv_prev(dof, eta, v_ord)
else:
    # Busted: called by compute_A_inv_with_dispersion.
    # [A, v1] = linalg.compute_A_inv(v_ord, n_dof=dof)
    [A_inv, v1] = lo.compute_A_inv_with_dispersion(v_ord, eta, n_dof=dof)

print("\nA:")
print(mat2txt(A))
[A, dnu] = lo.compute_A_CS(dof, A)
print("\ndnu:", dnu)
print("\nA_CS:")
print(mat2txt(A))
A_inv = np.linalg.inv(A)
print("\nR:")
print(mat2txt(A_inv @ M @ A))
exit()

r = calculate_radiation(
    acc, energy=2.5e9, calc_config=calc_config, install_radiators=True
)

exit()

use_tpsa = True
if not use_tpsa:
    ps = ss_vect_double()
    ps.set_zero()
    ps[phase_space_ind.x_] = 1e-3
else:
    ps = ss_vect_tps()
    ps.set_identity()


# First step:
#
# use closed orbit
# 1. calculate fix point and Poincarè Map M with damped system (i.e. radiation on
#    and cavity on (without dispersion in a second case)
# 2. diagonalise M = A $\Gamma$ A$^{-1}$
# 3. eigenvalues:
#        - complex part: tunes,
#        - real part: damping times  (refer equation)
#    use eigen values of symplectic matrix to identify the planes
# 4. propagate A, thin kick will create diffusion coeffs (don't forget to zero
#    them before calculation starts (sum it up afterwards


print(ps)
acc.propagate(calc_config, ps, 0, 2000)
print(ps)


if use_tpsa:
    # Inspect curly_H in
    for a_del in rad_del:
        name = a_del.getDelegatorName()
        idx = a_del.getDelegatorIndex()
        curly_H_x = a_del.getCurlydHx()
        txt = f"{name:10s} {idx:4d} curly_H_x {curly_H_x:5f}"
        print(txt)

    I = np.array([a_del.getSynchrotronIntegralsIncrements() for a_del in rad_del_kick])

    for a_del in rad_del_kick:
        name = a_del.getDelegatorName()
        idx = a_del.getDelegatorIndex()
        curly_H_x = a_del.getCurlydHx()
        dI = a_del.getSynchrotronIntegralsIncrements()
        D_rad = a_del.getDiffusionCoefficientsIncrements()

        txt = f"{name:10s} {idx:4d} curly_H_x {curly_H_x: 10.6e}"
        txt += "    dI " + ",".join(["{: 10.6e}".format(v) for v in dI])
        txt += "   "
        txt += "    D_rad" + ",".join(["{: 10.6e}".format(v) for v in D_rad])
        txt += "   "
        print(txt)
