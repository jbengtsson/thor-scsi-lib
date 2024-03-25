"""Use Case:
     Parametric scans/evaluations for a unit cell.
"""


import logging

# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="WARNING")
logger = logging.getLogger("thor_scsi")

from dataclasses import dataclass
import os

import numpy as np

from thor_scsi.utils import lattice_properties as lp


def set_phi_rb(lat_prop, phi_rb):
    # Optimum reverse bend angle is:
    #   phi_rb = -0.37
    phi_b  = 3.0
    b_name = ["b0", "b1", "b2", "b3", "b4", "b5"]
    b_scl = np.array(
        [1.094181/phi_b, 0.151199/phi_b, 0.151101/phi_b, 0.101861/phi_b,
         0.001569/phi_b, 0.000089/phi_b]
    )

    dphi = 3.0 - 2.0*phi_rb

    lat_prop.set_phi_fam("qf", phi_rb, True)
    for k in range(len(b_scl)):
        lat_prop.set_phi_fam(b_name[k], b_scl[k]*dphi, True)


def set_phi_var_bend_rad_rb(dphi):
    phi_b  = 3.0
    n_b    = 6

    phi_b = 0e0
    for k in range(1, n_b):
      phi_b += get_set.set_phi_fam(lat_prop._lattice, "b0", 0)
 
    get_set.set_phi_fam(lat_prop._lattice, "b0", dphi, True)
    get_set.set_phi_fam(lat_prop._lattice, "b1", dphi, True)
    get_set.set_phi_fam(lat_prop._lattice, "b2", dphi, True)
    get_set.set_phi_fam(lat_prop._lattice, "b3", dphi, True)
    get_set.set_phi_fam(lat_prop._lattice, "b4", dphi, True)
    get_set.set_phi_fam(lat_prop._lattice, "b5", dphi, True)
    get_set.set_phi_fam(lat_prop._lattice, "qf", dphi, True)


@dataclass
class MatchingState:
    n_iter: float = 0
    #: :math:`\chiˆ2`
    chi_2: float = np.nan


def set_phi_cell(name, phi):
    print("set_phi_cell: {name:8s} {phi:5.3f}")


# Global additional variables needed for optimisation function.
n_iter     = 0
chi_2_min  = 1e30
n_iter_min = 0
prms_min   = []
rbend      = ""

def opt_var_bend_radius(
        no, lat, param_list, C, bounds, phi_des, eps_x_des, nu_des,
                     beta_straight_des, weights, phi):
    """Use Case: optimise unit cell.
    """

    def compute_chi_2(
            no: int,
            lat: tslib.Accelerator,
            model_state: tslib.ConfigType,
            eps_x_des: float,
            nu_des, prms
    ) -> Tuple[float, float]:
        """Computes weightsed sum square
        """

        global n_iter
        global chi_2_min
        global n_iter_min
        global prms_min

        M_map = lo.compute_map(lat, model_state, desc=desc, tpsa_order=no)
        stable = lo.check_if_stable_3D(M_map.jacobian())

        if stable[0] and stable[1] and stable[2]:
            
            ind = np.zeros(nv, int)
            ind[delta_] = 1
            alpha_c_1 = M_map.ct.get(ind)/C
            ind[delta_] = 2
            alpha_c_2 = M_map.ct.get(ind)/C
            _, nu, xi = lo.compute_nu_xi(desc, no, M_map)

            M, A, data = \
                compute_periodic_solution(
                    lat, model_state, named_index, desc, False)

            beta_straight = np.zeros(2, dtype=float)
            beta_straight[X_] = data.twiss.sel(plane="x", par="beta").values[0]
            beta_straight[Y_] = data.twiss.sel(plane="y", par="beta").values[0]
            eta_straight_x    = \
                data.dispersion.sel(phase_coordinate="x").values[0]

            stable, M, cod, A, U_0, J, tau, eps, D_rad = \
                rad.compute_radiation(lat, model_state, 2.5e9, 1e-15, desc=desc)

            if stable:
                n = len(weights)
                dchi_2 = np.zeros(n, dtype=float)
                dchi_2[0] = \
                    weights["eps_x"] * np.sum((eps[X_] - eps_x_des) ** 2)
                dchi_2[1] = weights["nu_cell"] * np.sum((nu_cell - nu_des) ** 2)
                dchi_2[2] = \
                    weights["beta_straight"] \
                    * np.sum((beta_straight - beta_straight_des) ** 2)
                dchi_2[3] = weights["eta_straight_x"] * eta_straight_x ** 2
                dchi_2[4] = weights["xi"] * np.sum(xi ** 2)
                dchi_2[5] = weights["alpha_c_1"] / alpha_c_1 ** 2
                dchi_2[6] = weights["alpha_c_2"] * alpha_c_2 ** 2
                chi_2 = 0e0
                for k in range(n):
                    chi_2 += dchi_2[k]
            else:
                chi_2 = 1e30
        else:
            chi_2 = 1e30

        if chi_2 < chi_2_min:
            chi_2_min = chi_2
            n_iter_min = n_iter
            prms_min = prms
            print(f"\n{n_iter:4d} chi_2 = {chi_2:9.3e}\n\n  prms   =\n", prms)
            if stable:
                print(f"\n  dchi_2 =\n", dchi_2)
            print(f"\n  phi            = {compute_phi(lat):8.5f}")
            [b1, b2, b3] = \
                [param_list[0][0], param_list[1][0], param_list[2][0]]
            print(f"\n  {b1:5s}          = {get_phi_elem(lat, b1, 0):8.5f}")
            print(f"  {b2:5s}          = {get_phi_elem(lat, b2, 0):8.5f}")
            print(f"  {b3:5s}          = {get_phi_elem(lat, b3, 0):8.5f}")
            [b1, b2, b3, b4] = \
                [param_list[3][0], param_list[4][0], param_list[5][0],
                 param_list[6][0]]
            print(f"\n  {b1:5s}          = {get_phi_elem(lat, b1, 0):8.5f}")
            print(f"  {b2:5s}          = {get_phi_elem(lat, b2, 0):8.5f}")
            print(f"  {b3:5s}          = {get_phi_elem(lat, b3, 0):8.5f}")
            print(f"  {b4:5s}          = {get_phi_elem(lat, b4, 0):8.5f}")
            [b1, b2] = [rbend, param_list[8][0]]
            print(f"\n  {b1:5s}          = {get_phi_elem(lat, b1, 0):8.5f}")
            print(f"  {b2:5s}          = {get_phi_elem(lat, b2, 0):8.5f}")
            if stable:
                print(f"\n  eps_x          = {1e12*eps[X_]:5.3f}")
                print("  nu_cell        = [{:7.5f}, {:7.5f}]".
                      format(nu_cell[X_], nu_cell[Y_]))
                print("  beta_straight  = [{:7.5f}, {:7.5f}]".
                      format(beta_straight[X_], beta_straight[Y_]))
                print(f"  eta_straight_x = {eta_straight_x:10.3e}")
                print(f"  xi             = [{xi[X_]:5.3f}, {xi[Y_]:5.3f}]")
                print(f"  alpha_c_1      = {alpha_c_1:9.3e}")
                print(f"  alpha_c_2      = {alpha_c_2:9.3e}")
            else:
                print("\n Unstable.")
        return chi_2

    def f_unit_cell(prms):
        global n_iter
        global chi_2_min
        global prms_min
        global rbend

        n_iter += 1
        set_prm(lat, param_list, prms)
        dphi = \
            (phi
             - 8*(prms[0] + prms[1] + prms[2])
             - 2*(prms[3] + prms[4] + prms[5] + prms[6])
             - 2*prms[8])/8e0
        set_phi_fam(lat, rbend, dphi)
        chi_2 = compute_chi_2(no, lat, model_state, eps_x_des, nu_des, prms)
        return chi_2

    def prt_result(no, f_unit_cell, param_list, prms0, minimum):
        global n_iter
        global chi_2_min
        global n_iter_min
        global prms_min

        M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=no)
        stable = lo.check_if_stable(3, M.jacobian())
        if not stable:
            print("\nCell unstable:")
            print("  prms_min =", prms_min)
            print("  Setting prs_min for plot.")
            set_prm(lat, param_list, prms_min)
        # Compute new Twiss parameters along lattice.
        M, A, data = \
            compute_periodic_solution(
                lat, model_state, named_index, desc, False)
        plot_Twiss(data, "after.png", "Linear Optics - After")

        print("\nInitial parameter values:")
        n = n_iter_min
        n_iter = 0
        chi_2_min = 1e30
        f_unit_cell(prms0)
        print("\nFinal parameter values:")
        n_iter = n - 1
        chi_2_min = 1e30
        f_unit_cell(minimum["x"])
        print("\n Minimum:\n", minimum)

        print_Twiss(lat, data)

    global rbend

    max_iter = 1000
    f_tol    = 1e-4
    x_tol    = 1e-4

    print("\nopt_var_bend_radius:\n")
    # Initialise parameters.
    prms1 = prms0 = get_prm(lat, param_list)
    rbend = param_list[7][0]
    print("\nrbend = ", rbend)
    print("\nprms = ", prms1)

    # Methods:
    #   Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA,
    #   SLSQP, trust-constr, dogleg, truct-ncg, trust-exact, trust-krylov.

    # Powell ftol
    # CG     gtol
    minimum = optimize.minimize(
        f_unit_cell,
        prms1,
        method="CG",
        # bounds=bounds,
        # options={"xtol": x_tol, "maxiter": max_iter},
        options={"gtol": f_tol, "maxiter": max_iter},
    )

    prt_result(no, f_unit_cell, param_list, prms0, minimum)


# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 2
# Number of parameters.
nv_prm = 0
# Parameters max order.
no_prm = 0

cod_eps = 1e-15
E_0     = 3.0e9

home_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_4U")
file_name = os.path.join(home_dir, "max_4u.lat")

lat_prop = \
    lp.lattice_properties_class(nv, no, nv_prm, no_prm, file_name, E_0, cod_eps)

# Compute Twiss parameters along lattice.
stable = lat_prop.comp_per_sol()
print("\nCircumference [m]      = {:7.5f}".format(lat_prop.compute_circ()))
print("Total bend angle [deg] = {:7.5f}".
      format(lat_prop.compute_phi(lat_prop._lattice)))
lat_prop.prt_M()
if not stable:
    assert False
Twiss = lat_prop.get_Twiss(len(lat_prop._lattice)-1)
lat_prop.prt_Twiss_param(Twiss)

# Compute radiation properties.
stable = lat_prop.compute_radiation()
lat_prop.prt_rad()
lat_prop.prt_M_rad()

types = lat_prop.get_types()
lat_prop.prt_Twiss(types)

if False:
    lat_prop.plt_Twiss(types, False)

if not False:
    phi, eps_x, J_x, J_z, alpha_c = \
        lat_prop.unit_cell_rev_bend(15, -0.97, set_phi_rb)

    lat_prop.plt_scan_phi_rb(
        "plt_scan_phi_rb.png", phi, eps_x, J_x, J_z, alpha_c, True)
