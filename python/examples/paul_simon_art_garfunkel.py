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
from thor_scsi.utils.linear_optics import (
    compute_map,
    compute_nu_symp,
    compute_M_diag,
    compute_Twiss_along_lattice
)
from thor_scsi.utils.radiate import compute_radiation
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


def prt_fam(acc, fam_name):
    print("\nprt_fam:")
    for q in acc.elementsWithName(fam_name):
        print("  {:4s} {:3d} {:6.3f} {:4.2f}".
              format(q.name, q.index, q.getMultipoles().getMultipole(2).real,
                     q.getLength()))


def set_db_2_fam(acc, fam_name, db_2):
    for q in acc.elements_with_name(fam_name):
        b_2 = q.get_multipoles().get_multipole(2).real
        q.get_multipoles().set_multipole(2, b_2+db_2)


def set_db_2L_fam(acc, fam_name, db_2L):
    q = acc.find(fam_name, 0)
    L = q.get_length()
    set_db_2_fam(acc, q.name, db_2L/L)


def compute_dnu_db_2L(acc, calc_config, fam_name, db_2L):
    dof = 2
    nu = np.zeros([2, dof])
    for k in range(-1, 2, 2):
        set_db_2L_fam(acc, fam_name, k*db_2L)
        #M = map2numpy(compute_map(acc, calc_config))[:6, :6]
        a_map = compute_map(acc, calc_config, desc=desc)
        M = a_map.jacobian()
        set_db_2L_fam(acc, fam_name, -k*db_2L)
        nu[(k+1)//2] = compute_nu_symp(dof, M)

    dnu_db_2 = (nu[1]-nu[0])/(2e0*db_2L)
    return dnu_db_2


def tweak_nu(fam_names, dnu_x, dnu_y):
    n = len(fam_names)
    dnu = np.array([dnu_x, dnu_y])
    dnu_db_2L = np.zeros((2, n), dtype='float')
    for k in range(n):
        dnu_db_2L[:, k] = \
            compute_dnu_db_2L(acc, calc_config, fam_names[k],   1e-4)
    print("\ndnu_db_2L:\n", mat2txt(dnu_db_2L))
    u, s, v_t = np.linalg.svd(dnu_db_2L, full_matrices=False)
    dnu_db_2L_inv = (v_t.T @ np.diag(s**-1) @ u.T)
    db_2L = dnu_db_2L_inv @ dnu
    print("\ndb_2L = ", vec2txt(db_2L))
    for k in range(n):
        set_db_2L_fam(acc, fam_names[k], db_2L[k])


t_dir = \
    os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi", "JB",
                 "BESSY-III")
t_file = os.path.join(t_dir, "alsu-7ba-20180503c.lat")
# t_file = os.path.join(t_dir, "b3_sf_40Grad_JB.lat")

acc = accelerator_from_config(t_file)
calc_config = tslib.ConfigType()

cav = acc.find("cav", 0)
cav.set_harmonic_number(327)
cav.set_phase(180.0)
print("\nCavity", repr(cav))
print(f"""\nCavity info:
  f [MHz] {1e-6*cav.get_frequency()}
  V [MV]  {1e-6*cav.get_voltage()}
  h       {cav.get_harmonic_number()}
  phi     {cav.get_phase()}
""", end="")

dof = 2

E = 2.0e9
calc_config.radiation = False
calc_config.Cavity_on = False


desc = gtpsa.desc(6, 2)

M = compute_map(acc, calc_config, desc=desc)
jac = M.jacobian()
nu = compute_nu_symp(dof, jac)
print("\nM:\n" + mat2txt(jac))
print("\nnu = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))

stable, A_mat, _, _ = compute_M_diag(dof, jac)
A = gtpsa.ss_vect_tpsa(desc, 1)
A.set_jacobian(A_mat)

ds = compute_Twiss_along_lattice(2, acc, calc_config, A=A, desc=desc)
df = twiss_ds_to_df(ds)
# print("\nTwiss - ds:\n", ds)
# print("\nTwiss - df:\n", df)

with open("twiss.tsf", "wt") as fp:
    fp.write(df_to_tsv(df))
df.to_json("twiss.json")

if True:
    plt_twiss(ds)
    plt_curly_H(ds)


compute_radiation(acc, calc_config, E, 1e-15, desc=desc)

calc_config.radiation = False
calc_config.Cavity_on = False

# Tweak the tune using all the quadrupole families.
b_2 = ['q1', 'q2', 'q3', 'q4', 'mb3h', 'qfh', 'mqfh']
tweak_nu(b_2, 0.00501, -0.00670)
a_map = compute_map(acc, calc_config, desc=desc)
M = a_map.jacobian()
nu = compute_nu_symp(dof, M)
print("\nM:\n" + mat2txt(M))
print("\nnu = [{:7.5f}, {:7.5f}]".format(nu[X_], nu[Y_]))

plt.ioff()
plt.show()
nu
