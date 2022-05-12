import pytest

from thor_scsi.utils import linalg


def test_omega_matrix():
    mat = linalg.omega_block_matrix(3)
    nr, nc = mat.shape


def test_invert_matix():

    # linalg.compute_A_inv(v)
