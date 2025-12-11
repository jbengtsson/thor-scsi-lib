import pytest

from thor_scsi.utils import linalg


def test_omega_matrix():
    mat = linalg.omega_block_matrix(3)
    nr, nc = mat.shape
    assert nr == 6 and nc == 6


def test_invert_matix():
    # Placeholder: ensure module exposes omega_block_matrix and PlaneDistance
    assert hasattr(linalg, "omega_block_matrix")
    assert hasattr(linalg, "PlaneDistance")
