"""Support for formatting np.ndarray's in a Tracy compatible fashion
"""

import io
from contextlib import redirect_stdout

import numpy as np


def prt2txt(input):
    with io.StringIO() as buf, redirect_stdout(buf):
        print(input, end='')
        output = buf.getvalue()
    return output


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


def chop_array(arr: np.ndarray, eps: float, *, copy: bool = True) -> np.ndarray:
    """
    """
    if copy:
        narr = arr.copy()
    narr[np.absolute(narr) < eps] = 0e0
    return narr


def chop_vec(vec, eps):
    return chop_array(vec, eps)


def chop_mat(mat, eps):
    return chop_array(mat, eps)


def chop_cmplx_vec(vec, eps, *, copy=True):
    return chop_array(vec, eps, copy=copy)


def chop_cmplx_mat(mat, eps, *, copy=True):
    return chop_array(mat, eps, copy=copy)


def cmplx2txt(val: complex) -> str:
    """Convert a complex number to text

    Todo:
        Check if provided by standard package
    """
    if val.imag < 0e0:
        sign_str = "-"
    else:
        sign_str = "+"

    return f"{val.real:13.6e} {sign_str} {np.abs(val.imag):9.3e}i"


def vec2txt(vec: np.ndarray) -> str:
    """Convert a (numpy) vector to precisely formatted text
    """
    a = np.atleast_1d(vec)
    flag = np.iscomplexobj(a)
    if flag:
        return " ".join([cmplx2txt(v) for v in a])
    return " ".join([" {:13.6e}".format(v) for v in a])
    # return " ".join([" {:23.16e}".format(v) for v in a])


def mat2txt(mat: np.ndarray, *, name=None) -> str:
    """Convert a (numpy) matrix to precisely formatted text
    """
    a = np.atleast_2d(mat)
    txt = "\n".join([vec2txt(v) for v in a])
    if name is not None:
        txt = name + txt
    return txt


__all__ = ["prt2txt", "vec2txt", "mat2txt", "complex2txt"]
