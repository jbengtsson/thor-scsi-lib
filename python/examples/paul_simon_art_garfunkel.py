# Art - Philosophy - Science: Ancient Greece
#
# 𝐿𝑖𝑏𝑟𝑎𝑟𝑦 𝑜𝑓 𝐴𝑙𝑒𝑥𝑎𝑛𝑑𝑟𝑖𝑎
# https://www.bbc.co.uk/programmes/b00j0q53
#
# D. Knuth 𝑇ℎ𝑒 𝐴𝑟𝑡 𝑜𝑓 𝐶𝑜𝑚𝑝𝑢𝑡𝑒𝑟 𝑃𝑟𝑜𝑔𝑟𝑎𝑚𝑚𝑖𝑛𝑔
#
# 𝐼 𝑡ℎ𝑜𝑢𝑔ℎ𝑡 𝑡ℎ𝑎𝑡 𝐼 𝑤𝑎𝑠 𝑎 𝑝𝑒𝑟𝑓𝑒𝑐𝑡𝑖𝑜𝑛𝑖𝑠𝑡 𝑢𝑛𝑡𝑖𝑙 𝐼 𝑚𝑒𝑡 𝐾𝑛𝑢𝑡ℎ.
# F.C. Graham
# https://mathweb.ucsd.edu/~fan/paint/math.html
#
# M. Frangouli 𝐻𝑒𝑟𝑒'𝑠 𝑡𝑜 𝑡ℎ𝑒 𝐻𝑒𝑟𝑜𝑒𝑠
# https://youtu.be/RHOj622NMwo
#
# Vangelis 𝐶ℎ𝑎𝑟𝑖𝑜𝑡𝑠 𝑜𝑓 𝐹𝑖𝑟𝑒
# https://youtu.be/CSav51fVlKU
# https://youtu.be/CwzjlmBLfrQ
#
# 𝑇ℎ𝑒 𝑅𝑒𝑎𝑙 𝐶ℎ𝑎𝑟𝑖𝑜𝑡𝑠 𝑜𝑓 𝐹𝑖𝑟𝑒
# https://youtu.be/zox9Q_Ljtiw?t=771
# https://youtu.be/K6QSsOwnHUQ
# https://youtu.be/aY2eNtg6X2k
# https://youtu.be/Ilceqb9YxFc
#
# Statue 𝑆𝑝𝑖𝑟𝑖𝑡 𝑜𝑓 𝑡ℎ𝑒 𝑀𝑎𝑟𝑎𝑡ℎ𝑜𝑛, Hopkinton, MA & sister town Marathon, Greece.
#
# Kathrine Switzer Boston Marathon, 1968.
# https://kathrineswitzer.com/1967-boston-marathon-the-real-story
#
# Stylianos Kyriakides Boston Marathon, 1946.
#
# Simon & Garfunkel 𝑇ℎ𝑒 𝐵𝑜𝑥𝑒𝑟
# https://youtu.be/l3LFML_pxlY


import logging
# Levels: DEBUG, INFO, WARNING, ERROR, and CRITICAL.
logging.basicConfig(level="ERROR")

import gtpsa

import os
import numpy as np
import matplotlib.pyplot as plt

import thor_scsi.lib as tslib

from thor_scsi.factory import accelerator_from_config

from thor_scsi.utils import linear_optics as lo, courant_snyder as cs, \
    radiate as rad

from thor_scsi.utils.phase_space_vector import map2numpy
from thor_scsi.utils.output import mat2txt, vec2txt
from thor_scsi.utils.twiss_output import twiss_ds_to_df, df_to_tsv

# Turn on interactive mode (default is off).
plt.ion()

X_, Y_, Z_ = [
    tslib.spatial_index.X,
    tslib.spatial_index.Y,
    tslib.spatial_index.Z
]

def plt_twiss(ds):
    fig, (gr_1, gr_2) = plt.subplots(2)

    gr_1.set_title("Linear Optics")
    gr_1.set_xlabel("s [m]")
    gr_1.set_ylabel(r"$\beta_{x,y}$ [m]")
    gr_1.plot(ds.s, ds.twiss.sel(plane="x", par="beta"), label=r"$\beta_x$")
    gr_1.plot(ds.s, ds.twiss.sel(plane="y", par="beta"), label=r"$\beta_y$")
    gr_1.legend()

    gr_2.set_xlabel("s [m]")
    gr_2.set_ylabel(r"$\eta_x$ [m]")
    gr_2.plot(ds.s, ds.dispersion.sel(phase_coordinate="x"), label=r"$\eta_x$")
    fig.tight_layout()


def plt_curly_H(ds):
    fig, gr = plt.subplots(1)
    gr.set_title("curly_H")
    gr.set_xlabel(r"$\eta_x$ [m]")
    gr.set_ylabel(r"$\eta_x'$ [rad]")
    gr.plot(ds.dispersion.sel(phase_coordinate="x"),
            ds.dispersion.sel(phase_coordinate="px"))
    fig.tight_layout()


def prt_fam(lat, fam_name):
    print("\nprt_fam:")
    for q in lat.elementsWithName(fam_name):
        print("  {:4s} {:3d} {:6.3f} {:4.2f}".
              format(q.name, q.index, q.getMultipoles().getMultipole(2).real,
                     q.getLength()))


def set_db_2_fam(lat, fam_name, db_2):
    for q in lat.elements_with_name(fam_name):
        b_2 = q.get_multipoles().get_multipole(2).real
        q.get_multipoles().set_multipole(2, b_2+db_2)


def set_db_2L_fam(lat, fam_name, db_2L):
    q = lat.find(fam_name, 0)
    L = q.get_length()
    set_db_2_fam(lat, q.name, db_2L/L)


def compute_dnu_db_2L(lat, model_state, fam_name, db_2L):
    n_dof = 2
    nu = np.zeros([2, n_dof])
    for k in range(-1, 2, 2):
        set_db_2L_fam(lat, fam_name, k*db_2L)
        #M = map2numpy(compute_map(lat, model_state))[:6, :6]
        a_map = lo.compute_map(lat, model_state, desc=desc)
        M = a_map.jacobian()
        set_db_2L_fam(lat, fam_name, -k*db_2L)
        nu[(k+1)//2] = lo.compute_nu_symp(n_dof, M)

    dnu_db_2 = (nu[1]-nu[0])/(2e0*db_2L)
    return dnu_db_2


def tweak_nu(fam_names, dnu_x, dnu_y):
    n = len(fam_names)
    dnu = np.array([dnu_x, dnu_y])
    dnu_db_2L = np.zeros((2, n), dtype='float')
    for k in range(n):
        dnu_db_2L[:, k] = \
            compute_dnu_db_2L(lat, model_state, fam_names[k],   1e-4)
    print("\ndnu_db_2L:\n", mat2txt(dnu_db_2L))
    u, s, v_t = np.linalg.svd(dnu_db_2L, full_matrices=False)
    dnu_db_2L_inv = (v_t.T @ np.diag(s**-1) @ u.T)
    db_2L = dnu_db_2L_inv @ dnu
    print("\ndb_2L = ", vec2txt(db_2L))
    for k in range(n):
        set_db_2L_fam(lat, fam_names[k], db_2L[k])


def print_Twiss(str, Twiss):
    """

    todo:
        rename str e.g. header? prefix?

    """
    # eta, alpha, beta = Twiss[0], Twiss[1], Twiss[2]
    # that way I also check that Twiss has exactly three parameters
    eta, alpha, beta = Twiss
    print(str, end="")
    print(f"  eta    = [{eta[X_]:9.3e}, {eta[Y_]:9.3e}]")
    print(f"  alpha  = [{alpha[X_]:9.3e}, {alpha[Y_]:9.3e}]")
    print(f"  beta   = [{beta[X_]:5.3f}, {beta[Y_]:5.3f}]")


def compute_periodic_solution(lat, model_state, named_index, desc):
    """
    Todo:
        model_state: rename to calculation_configuration or calc_config
    """
    # Compute the periodic solution for a super period.
    # Degrees of freedom - RF cavity is off; i.e., coasting beam.
    n_dof = 2
    model_state.radiation = False
    model_state.Cavity_on = False

    stable, M, A = lo.compute_map_and_diag(n_dof, lat, model_state, desc=desc)
    res = cs.compute_Twiss_A(A)
    Twiss = res[:3]
    print_Twiss("Twiss:\n", Twiss)
    A_map = gtpsa.ss_vect_tpsa(desc, no)
    A_map.set_jacobian(A)
    ds = \
        lo.compute_Twiss_along_lattice(
            n_dof, lat, model_state, A=A_map, desc=desc, mapping=named_index)

    return M, A, ds


t_dir = os.path.join(
    os.environ["HOME"], "Nextcloud", "thor_scsi", "JB", "BESSY-III")
t_file = os.path.join(t_dir, "alsu-7ba-20180503c.lat")
# t_file = os.path.join(t_dir, "b3_sf_40Grad_JB.lat")

lat = accelerator_from_config(t_file)
model_state = tslib.ConfigType()

cav = lat.find("cav", 0)
cav.set_harmonic_number(327)
cav.set_phase(180.0)
print("\nCavity", repr(cav))
print(f"""\nCavity info:
  f [MHz] {1e-6*cav.get_frequency()}
  V [MV]  {1e-6*cav.get_voltage()}
  h       {cav.get_harmonic_number()}
  phi     {cav.get_phase()}
""", end="")

n_dof = 2

# Number of phase-space coordinates.
nv = 7
# Variables max order.
no = 2

E = 2.0e9

model_state.radiation = False
model_state.Cavity_on = False

named_index = gtpsa.IndexMapping(dict(x=0, px=1, y=2, py=3, delta=4, ct=5))

desc = gtpsa.desc(nv, no)

# Compute linear map & diagonalise.
M, A, data = \
    compute_periodic_solution(lat, model_state, named_index, desc)

# Compute linear chromaticity from trace of 2nd order map.
M = lo.compute_map(lat, model_state, desc=desc, tpsa_order=2)
stable, nu, xi = lo.compute_nu_xi(desc, no, M)

print("\nM:\n", mat2txt(M.jacobian()[:6, :6]), "\n")
print("  nu = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))
print("  xi = [{:7.5f}, {:7.5f}]".format(xi[X_], xi[Y_]))

stable, A, A_inv, _ = lo.compute_M_diag(n_dof, M.jacobian())

ds = \
    lo.compute_Twiss_along_lattice(
        n_dof, lat, model_state, desc=desc, mapping=named_index)
df = twiss_ds_to_df(ds)
# print("\nTwiss - ds:\n", ds)
# print("\nTwiss - df:\n", df)

with open("twiss.tsf", "wt") as fp:
    fp.write(df_to_tsv(df))
df.to_json("twiss.json")

if True:
    plt_twiss(ds)
    plt_curly_H(ds)

rad.compute_radiation(lat, model_state, E, 1e-15, desc=desc)

model_state.radiation = False
model_state.Cavity_on = False

# Tweak the tune using all the quadrupole families.
b_2 = ['q1', 'q2', 'q3', 'q4', 'mb3h', 'qfh', 'mqfh']
tweak_nu(b_2, 0.00501, -0.00670)
a_map = lo.compute_map(lat, model_state, desc=desc)
M = a_map.jacobian()[:6, :6]
nu = lo.compute_nu_symp(n_dof, M)
print("\nM:\n" + mat2txt(M))
print("\nnu = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))

plt.ioff()
plt.show()
