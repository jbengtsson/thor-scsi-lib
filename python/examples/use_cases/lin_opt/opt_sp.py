
#!/usr/bin/env python3
"""
Optimize a higher-order achromat using one superperiod.
- Live interactive chi² plot during optimization
- Final dashboard figure (chi² + constraints)
- Self-contained single Python file
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import ClassVar
from scipy import optimize as opt
import logging

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

from thor_scsi.utils import lattice_properties as lp
from thor_scsi.utils import nonlin_dyn as nld_class
from thor_scsi.utils import linear_optics as lo, index_class as ind, prm_class as pc

from thor_scsi.utils import index_class as ind

ind = ind.index_class()


# -----------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------

PRM_RANGE = {
    "phi":       [0.0, 2.0],
    "phi_rbend": [-0.3, 0.3],
    "b_2":       [-10.0, 10.0],
    "b_2_bend":  [-1.5, -1.5]
}

DESIGN_VAL = {
    "eps_x_des"    : 50e-12,
    "phi_1_des"    : 1.3,
    "phi_rb_1_des" : -0.2,
    "b_2_des"      : 2.0,
    "nu_uc_des"    : np.array([0.4, 0.1]),
    "nu_sp_des"    : np.array([58.20/20.0, 17.28/20.0]),
    "beta_des"     : [5.0, 3.0],
    "dnu_des"      : [0.75, 0.25]
}

def get_weights():
    return {
        "eps_x"       : 1e17,
        "dphi"        : 1e-1,
        "phi_1"       : 0e-2,
        "phi_rb"      : 0e-3,
        "b_2"         : 0e-3,
        "alpha^(1)_c" : 1e-14,
        "alpha^(2)_c" : 1e1,
        "U_0"         : 1e-15,
        "etap_x_uc"   : 1e2,
        "alpha_uc"    : 1e-1,
        "nu_uc_x"     : 1e-2,
        "nu_uc_y"     : 1e-2,
        "eta_x"       : 1e2,
        "nu_sp_x"     : 1e0,
        "nu_sp_y"     : 1e0,
        "beta_x"      : 0e-6,
        "beta_y"      : 0e-6,
        "dnu_x"       : 1e-3,
        "dnu_y"       : 1e-3,
        "xi"          : 1e-7,
        "eta^(2)_x"   : 1e-6
    }

COD_EPS = 1e-10
E_0     = 3.0e9
A_MAX   = np.array([6e-3, 3e-3])
DELTA_MAX = 3e-2
BETA_INJ  = np.array([3.0, 3.0])

# -----------------------------------------------------------------------
# LATTICE SETUP
# -----------------------------------------------------------------------

def setup_lattice(lat_file: str):
    lat_prop = lp.lattice_properties_class(lat_file, E_0, COD_EPS, 2)
    b_3_list = ["s1_h3", "s2_h3", "s3_h3", "s4_h3"]
    nld = nld_class.nonlin_dyn_class(lat_prop, A_MAX, BETA_INJ, DELTA_MAX, b_3_list)
    nld.zero_mult(3)
    nld.zero_mult(4)
    if not lat_prop.comp_per_sol():
        raise RuntimeError("Lattice unstable in periodic solution")
    if not lat_prop.compute_radiation():
        raise RuntimeError("Lattice unstable in radiation computation")
    return lat_prop, nld

def find_uc_and_sp(lat_prop):
    uc_list = np.array([lat_prop._lattice.find("d2_h2_sl_d0a", 0).index, lat_prop._lattice.find("d2_h2_sl_d0a", 2).index])
    sp_list = np.array([lat_prop._lattice.find("lsborder", 0).index, lat_prop._lattice.find("lsborder", 1).index])
    return uc_list, sp_list

def get_bends(lat_prop):
    bend_names = ["d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3"]
    rbend_names = ["r1_h2", "r2_h2"]
    bend_list = [pc.bend_class(lat_prop, [name]) for name in bend_names]
    return bend_list, rbend_names

def get_prms(bend_list, rbend_list, lat_prop):
    prm = [("q1_h2", "b_2", PRM_RANGE["b_2"]), ("q2_h2", "b_2", PRM_RANGE["b_2"])]
    for b in bend_list:
        prm.append((b, "b_2_bend", PRM_RANGE["b_2_bend"]))
        prm.append((b, "phi_bend", PRM_RANGE["phi"]))
    for r in rbend_list:
        prm.append((r, "b_2", PRM_RANGE["b_2"]))
        prm.append((r, "phi", PRM_RANGE["phi_rbend"]))
    prm_list = pc.prm_class(lat_prop, prm)
    dprm_list = np.full(len(prm), 1e-4)
    return prm_list, dprm_list

# -----------------------------------------------------------------------
# OPTIMIZER
# -----------------------------------------------------------------------

class OptSP:
    def __init__(self, prm_class):
        self.lat_prop = prm_class.lattice[0]
        self.nld = prm_class.lattice[1]
        self.uc_list = prm_class.s_loc[0]
        self.sp_list = prm_class.s_loc[1]
        self.des_val_list = prm_class.design_vals
        self.weights = prm_class.weights
        self.bend_list = prm_class.dipoles[0]
        self.rbend_list = prm_class.dipoles[1]
        self.prm_list = prm_class.params[0]
        self.dprm_list = prm_class.params[1]
        self.history = []
        self.constr_history = []
        self.chi_2_min = 1e30
        self.n_iter = -1
        plt.ion()
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], label="chi²")
        self.ax.set_yscale("log")
        self.ax.set_xlabel("Iteration")
        self.ax.set_ylabel("chi²")
        self.ax.set_title("Live Convergence")
        self.ax.grid(True, which='both', ls='--')
        self.fig.show()
        self.fig.canvas.draw()

    def live_update(self):
        self.line.set_data(np.arange(len(self.history)), self.history)
        self.ax.relim()
        self.ax.autoscale_view()
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def compute_lat_prop(self):
        return self.lat_prop.comp_per_sol()

    def compute_constr(self):
        self.constr = {}
        eps_x = self.lat_prop._eps[ind.X]
        self.constr["eps_x"] = (eps_x - self.des_val_list["eps_x_des"])**2
        return self.constr

    def compute_chi_2(self):
        self.compute_constr()
        return sum(self.weights[k] * self.constr.get(k, 0) for k in self.weights)

    def f_sp(self, prm):
        self.n_iter += 1
        self.prm_list.set_prm(prm)
        if not self.compute_lat_prop():
            self.history.append(1e30)
            self.constr_history.append({})
            self.live_update()
            return 1e30
        chi_2 = self.compute_chi_2()
        self.history.append(chi_2)
        self.constr_history.append(self.constr.copy())
        self.live_update()
        if chi_2 < self.chi_2_min:
            logging.info(f"Iteration {self.n_iter}, chi² = {chi_2:.3e}")
            self.chi_2_min = chi_2
        return chi_2

    def run(self):
        prm, bounds = self.prm_list.get_prm()
        return opt.minimize(self.f_sp, prm, method="CG", bounds=bounds, options={"maxiter": 10000, "gtol": 1e-7, "eps": self.dprm_list})

# -----------------------------------------------------------------------
# DASHBOARD
# -----------------------------------------------------------------------

def plot_dashboard(history, constr_history, outname="dashboard.png"):
    keys = list(constr_history[0].keys())
    data = {k: [ch.get(k,1e-30) for ch in constr_history] for k in keys}
    fig, axs = plt.subplots(1,2,figsize=(12,5))
    axs[0].plot(history, marker='o')
    axs[0].set_yscale('log')
    axs[0].set_xlabel("Iteration")
    axs[0].set_ylabel("chi²")
    axs[0].set_title("Convergence of chi²")
    axs[0].grid(True, which='both', ls='--')
    for k,v in data.items():
        axs[1].plot(v, label=k)
    axs[1].set_yscale('log')
    axs[1].set_xlabel("Iteration")
    axs[1].set_ylabel("Constraint value")
    axs[1].set_title("Constraints")
    axs[1].grid(True, which='both', ls='--')
    axs[1].legend(fontsize='small')
    plt.suptitle("Superperiod Optimization Dashboard")
    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig(outname,dpi=150)
    logging.info(f"Dashboard saved to {outname}")
    plt.show()

# -----------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lattice", help=".lat file name (without extension)")
    parser.add_argument("--dashboard", default="dashboard.png", help="Dashboard output filename")
    args = parser.parse_args()
    home_dir = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "MAX_IV")
    lat_file = os.path.join(home_dir, args.lattice + ".lat")
    lat_prop, nld = setup_lattice(lat_file)
    uc_list, sp_list = find_uc_and_sp(lat_prop)
    bend_list, rbend_list = get_bends(lat_prop)
    prm_list, dprm_list = get_prms(bend_list, rbend_list, lat_prop)
    weight_list = get_weights()
    @dataclass
    class PrmClass:
        lattice:     ClassVar[list] = [lat_prop, nld]
        s_loc:       ClassVar[list] = [uc_list, sp_list]
        design_vals: ClassVar[dict] = DESIGN_VAL
        weights:     ClassVar[list] = weight_list
        dipoles:     ClassVar[list] = [bend_list, rbend_list]
        params:      ClassVar[list] = [prm_list, dprm_list]
    opt = OptSP(PrmClass)
    result = opt.run()
    print(result)
    plot_dashboard(opt.history, opt.constr_history, outname=args.dashboard)

if __name__ == "__main__":
    main()
