import pytest

from thor_scsi.utils.phase_space_vector import omega_block_matrix
import numpy as np


def test_omega_matrix():
    mat = omega_block_matrix(3)
    nr, nc = mat.shape

    mat_check = mat.copy()
    # Check for ones and zero out ... prepare for final test
    for i in range(3):
        index = i * 2 + 1

        upper = mat[index - 1, index]
        assert upper == pytest.approx(1, 1e-12)
        mat[index - 1, index] = 0
        del upper

        lower = mat[index, index - 1]
        assert lower == pytest.approx(-1, 1e-12)
        mat[index, index - 1] = 0
        del lower

    print(mat, mat_check)
    assert np.sum(np.absolute(mat.ravel())) == pytest.approx(0, abs=1e-12)
