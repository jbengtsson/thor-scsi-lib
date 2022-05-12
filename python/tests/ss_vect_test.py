import pytest
import thor_scsi.lib as tlib
from thor_scsi.utils import phase_space_vector
import numpy as np


def test20_tps_access():
    tps = tlib.tps()
    assert tps.peek([0] * 6) == 0
    return

    # These lines let the test crash currently
    idx = [0] * 6
    idx[3] = 5
    assert tps.peek(idx) == 0

    tps.pook(idx, 42)
    assert tps.peek(idx) == 42
    print(tps)
    assert 1 == 0


def test25_sst_vect_access():
    t_map = tlib.ss_vect_tps()
    t_map.set_identity()

    tps = t_map[0]
    assert tps[0] == 1
    assert tps[1] == 0

    tps = t_map[1]
    assert tps[0] == 0
    assert tps[1] == 1
    assert tps[2] == 0

    print(tps.cst())
    assert tps.cst() == 0


def test30_element_access():
    """Test accessing element directly

    Todo:
        Clean up of the interface
    """
    t_map = tlib.ss_vect_tps()
    t_map.set_identity()

    print(t_map)
    # assert t_map[0, 1] == 1
    # assert t_map[1, 2] == 1
    # assert t_map[2, 3] == 1
    # assert t_map[3, 4] == 1
    # assert t_map[4, 5] == 1
    # assert t_map[5, 6] == 1


def test40_arma_support():
    """Test that numpy arrays can be transparantly handed by arma
    """
    npm = np.eye(6)
    mat = tlib.Matrix(npm)
    mat.swap_rows(4, 5)
    mat.swap_cols(5, 4)

    test = mat - np.eye(6)
    assert test.sum() == 0


def test40_ps_swap():
    """Test that identity matrix is an identity matrix again
    """
    npm = np.eye(6)
    mat = phase_space_vector.jac_intern_to_std(npm)
    test = mat - np.eye(6)
    assert test.sum() == 0


def test41_omega_matrix():
    """Omega matrix created correctly
    """
    a = phase_space_vector.omega_block_matrix(3)
    print("omega?", a)
    assert a[0, 1] == 1
    a[0, 1] = 0
    assert a[1, 0] == -1
    a[1, 0] = 0

    assert a[2, 3] == 1
    a[2, 3] = 0
    assert a[3, 2] == -1
    a[3, 2] = 0
    assert a[4, 5] == 1
    a[4, 5] = 0
    assert a[5, 4] == -1
    a[5, 4] = 0

    print("zero?", a)
    assert np.absolute(a).sum() == 0


def test42_ps_offset_identity():
    """Test that identity matrix is an identity matrix again
    """
    npm = np.eye(6)
    mat = phase_space_vector.jac_intern_to_std(npm)
    assert np.absolute(mat - np.eye(6)).sum() == 0


def test43_ps_offset_marked():
    """Test that identity matrix is an identity matrix again
    """
    npm = np.zeros([6, 6])
    marker = np.eye(4, 4, 1) * 5 + np.eye(4, 4, -1) * -7
    npm[:4, :4] = marker
    npm[4, 5] = 2
    npm[5, 4] = 3
    print(npm)
    mat = phase_space_vector.jac_intern_to_std(npm)
    print(mat)
    assert mat[5, 4] == 2
    mat[5, 4] = 0
    assert mat[4, 5] == 3
    mat[4, 5] = 0
    mat[:4, :4] -= marker
    assert np.absolute(mat).sum() == 0


def test50_ss_vect_arma():
    """Test access of ss_vect as armadillo matrix
    """
    return

    t_map = tlib.ss_vect_tps()
    t_map.set_identity()

    a_mat = tlib.Matrix(np.eye(6))
    print(type(a_mat), dir(a_mat))
    mat = tlib.ss_vect_tps_to_mat(t_map)
    mat = np.array(mat, copy=False)[:6, :6]
    print(type(mat), dir(mat))
    mat.swap_rows(5, 6)
    print(mat)
    mat.swap_cols(5, 6)
    print(mat)
    test = mat - np.eye(6)
    assert np.absolute(test).sum() == 0
    # assert 0


def mat2map(M):
    """Check that functions are equivalent
    """
    [Id, t_map] = [tlib.ss_vect_tps(), tlib.ss_vect_tps()]
    Id.set_identity()
    t_map.set_zero()
    for j in range(len(M)):
        for k in range(len(M[0])):
            t_map[j] += M[j][k] * Id[k]
    return t_map


def test60_sst_vect_arma():
    """
    """
    idx = np.arange(4)
    test = np.zeros([7, 7])
    test[idx + 2, idx] = np.arange(1, 4 + 1) * 3
    test[1, 3] = 3
    print(test)
    t = mat2map(test[:6, :6])
    t2 = tlib.mat_to_ss_vect_tps(tlib.Matrix(test))

    for i in range(len(t)):
        for j in range(6):
            assert test[i, j] == t[i][j]
            assert test[i, j] == t2[i][j]


def test61_sst_vect_arma():
    test = np.random.uniform(size=[7, 7])
    t = mat2map(test[:6, :6])
    t2 = tlib.mat_to_ss_vect_tps(tlib.Matrix(test))

    for i in range(6):
        for j in range(6):
            assert test[i, j] == t[i][j]
            assert test[i, j] == t2[i][j]
