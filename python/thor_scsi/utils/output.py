"""Support for formatting np.ndarray's in a Tracy compatible fashion
"""
import numpy as np


def sign(x: float) -> int:
    """Compute
    """

    return np.sign(x).astype(np.int)

    if x > 0e0:
        return 1
    elif x < 0e0:
        return -1
    else:
        return 0


def chop_vec(vec, eps):
    for k in range(vec.size):
        if np.abs(vec[k]) < eps:
            vec[k] = 0e0
    return vec


def chop_mat(mat, eps):
    for k in range(mat[:, 0].size):
        chop_vec(mat[k, :], eps)
    return mat


def chop_cmplx_vec(vec, eps):
    for k in range(vec.size):
        [x, y] = [vec[k].real, vec[k].imag]
        if np.abs(x) < eps:
            x = 0e0
        if np.abs(y) < eps:
            y = 0e0
        vec[k] = complex(x, y)
    return vec


def chop_cmplx_mat(mat, eps):
    for k in range(mat[:, 0].size):
        chop_cmplx_vec(mat[k, :], eps)
    return mat


def complex2txt(val: complex) -> str:
    """Convert a complex number to text

    Todo:
        Check if provided by standard package
    """
    if val.imag < 0e0:
        sign_str = "-"
    else:
        sign_str = "+"

    return f"{val.real:12.3e} {sign_str} {np.abs(val.imag):9.3e}i"


def vec2txt(vec: np.ndarray) -> str:
    """Convert a (numpy) vector to precisely formatted text
    """
    a = np.atleast_1d(vec)
    flag = np.iscomplexobj(a)
    if flag:
        return " ".join([complex2txt(v) for v in a])
    return " ".join(["{:15.6e}".format(v) for v in a])


def mat2txt(mat: np.ndarray) -> str:
    """Convert a (numpy) matrix to precisely formatted text
    """
    a = np.atleast_2d(mat)
    return "\n".join([vec2txt(v) for v in a])


__all__ = ["vec2txt", "mat2txt", "complex2txt"]
