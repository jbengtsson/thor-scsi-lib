"""Test of Unit Case prototype code:

    * Read lattice file.
    * Compute twiss parameters.
    * Compute closed orbit.
"""

import os
import numpy as np

from thor_scsi.lib import (
    ss_vect_tps,
    ConfigType
)

from examples import twiss_jb as tw

from thor_scsi.factory import accelerator_from_config
#from thor_scsi.utils import linear_optics
#from thor_scsi.utils.extract_info import accelerator_info


[X_, Y_, Z_]                    = [0, 1, 2]
[x_, px_, y_, py_, ct_, delta_] = [0, 1, 2, 3, 5, 4]
[Dip, Quad, Sext]               = [1, 2, 3]


def test_stuff(n, conf):
    map = ss_vect_tps()

    map.set_identity()
    acc.propagate(conf, map,  n, len(acc))
    acc.propagate(conf, map,  0, n)
    M = tw.get_mat(map)
    # Reduce matrix from 7x7 to 6x6.
    M = M[0:6, 0:6]
    tw.prt_np_mat("\nPoincaré Map:\n", "14.6e", M)

    m51 = M[x_][x_]*M[px_][delta_] - M[px_][x_]*M[x_][delta_]
    m52 = M[x_][px_]*M[px_][delta_] - M[px_][px_]*M[x_][delta_]
    print(f"\nm_15 = {m51:21.14e} m_16 = {m52:21.14e}")
    print(f"m_15 = {M[ct_][x_]:21.14e} m_16 = {M[ct_][px_]:21.14e}")

    A = tw.compute_ring_twiss(M)
    tw.prt_np_mat("\nA:\n", "14.6e", A)

    n_dof = 2
    n     = 2*n_dof

    [nus, stable] = tw.compute_nus(n_dof, M)
    print(f"\ncompute_nus:\n  nu    = [{nus[X_]:5.3f}, {nus[Y_]:5.3f}]")

    [nus, stable] = tw.compute_nus_symp_mat(n_dof, M)
    print(f"\ncompute_nus_symp_mat:\n  nu    = [{nus[X_]:5.3f}, {nus[Y_]:5.3f}]")

    [eta, alpha, beta, nu, stable] = tw.compute_twiss_M(M)
    tw.prt_twiss("\ncompute_twiss_M:\n", eta, alpha, beta, nu)

    [eta, alpha, beta, nu] = tw.compute_twiss_A(A)
    tw.prt_twiss("\ncompute_twiss_A:\n", eta, alpha, beta, nu)

    [eta, alpha, beta, nu] = tw.compute_twiss_A_A_tp(A)
    tw.prt_twiss("\ncompute_twiss_A_A_tp:\n", eta, alpha, beta, nu)

    # Cross check.
    tw.prt_np_mat("\nA^-1*M*A:\n", "14.6e",
                  np.linalg.multi_dot([np.linalg.inv(A), M, A]))


t_dir  = os.path.join(os.environ["HOME"], "Nextcloud", "thor_scsi")
t_file = os.path.join(t_dir, "b3_tst.lat")

acc  = accelerator_from_config(t_file)
conf = ConfigType()

if False:
    test_stuff(15, conf)

cod, found, s_loc = tw.compute_closed_orbit(acc, conf, 0e0, 10, 1e-10)
elem = acc.find("uq1", 1)
#print(repr(elem.getMultipoles()))
muls = elem.getMultipoles()
muls.setMultipole(Dip, 1e-3 - 1e-3j)
print("\n", repr(elem))
cod, found, s_loc = tw.compute_closed_orbit(acc, conf, 0e0, 10, 1e-10)
exit()

t_map = tw.compute_map(acc, conf)
M = tw.get_mat(t_map)

# Reduce matrix from 7x7 to 6x6.
M = M[0:6, 0:6]
tw.prt_np_mat("\nPoincaré Map:\n", "14.6e", M)

[eta, alpha, beta, nu, stable] = tw.compute_twiss_M(M)
tw.prt_twiss("\ncompute_twiss_M:\n", eta, alpha, beta, nu)

A = tw.compute_ring_twiss(M)
tw.prt_np_mat("\nA:\n", "14.6e", A)

# Cross check.
tw.prt_np_mat("\nA^-1*M*A:\n", "14.6e",
              np.linalg.multi_dot([np.linalg.inv(A), M, A]))

tw.compute_twiss_lat("linlat.out", acc, conf, tw.get_map(A))

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
