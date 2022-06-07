import pytest
import thor_scsi.lib as tslib
import numpy as np


def test00_instantiate_multipoles():
    n_max = 20
    mul = tslib.TwoDimensionalMultipoles(n_max)


def test10_multiply_multipoles():
    mul = tslib.TwoDimensionalMultipoles()
    mul.setMultipole(1, 2)
    mul.setMultipole(2, 4)

    assert mul.getMultipole(1).real == pytest.approx(2)
    assert mul.getMultipole(2).real == pytest.approx(4)
    print(mul)

    mul *= 2
    assert mul.getMultipole(1).real == pytest.approx(4)
    assert mul.getMultipole(2).real == pytest.approx(8)


def test11_multiply_multipoles_scalar_int():
    mul = tslib.TwoDimensionalMultipoles()
    mul.setMultipole(1, 2)
    mul.setMultipole(2, 4)

    assert mul.getMultipole(1).real == pytest.approx(2)
    assert mul.getMultipole(2).real == pytest.approx(4)
    print(mul)

    n_mul = mul * 2
    # Former multipoles still the same
    assert mul.getMultipole(1).real == pytest.approx(2)
    assert mul.getMultipole(2).real == pytest.approx(4)

    assert n_mul.getMultipole(1).real == pytest.approx(4)
    assert n_mul.getMultipole(2).real == pytest.approx(8)


def test11_multiply_multipoles_scalar_float():
    mul = tslib.TwoDimensionalMultipoles()
    mul.setMultipole(1, 2)
    mul.setMultipole(2, 4)

    assert mul.getMultipole(1).real == pytest.approx(2)
    assert mul.getMultipole(2).real == pytest.approx(4)
    print(mul)

    n_mul = mul * float(2.0)
    # Former multipoles still the same
    assert mul.getMultipole(1).real == pytest.approx(2)
    assert mul.getMultipole(2).real == pytest.approx(4)

    assert n_mul.getMultipole(1).real == pytest.approx(4)
    assert n_mul.getMultipole(2).real == pytest.approx(8)

def test12_multiply_multipoles_scalar_float():
    """float on the other side ...

    Commutative law needs to be explicitly implemented
    """
    mul = tslib.TwoDimensionalMultipoles()
    mul.setMultipole(1, 2)
    mul.setMultipole(2, 4)

    assert mul.getMultipole(1).real == pytest.approx(2)
    assert mul.getMultipole(2).real == pytest.approx(4)
    print(mul)

    n_mul = float(2.0) * mul
    # Former multipoles still the same
    assert mul.getMultipole(1).real == pytest.approx(2)
    assert mul.getMultipole(2).real == pytest.approx(4)

    assert n_mul.getMultipole(1).real == pytest.approx(4)
    assert n_mul.getMultipole(2).real == pytest.approx(8)


def test20_multiply_multipoles_scalar_vector():
    h_max = 4
    mul = tslib.TwoDimensionalMultipoles(h_max)
    refs = np.array([2, 3, 5, 7])
    mul.setMultipole(1, refs[0])
    mul.setMultipole(2, refs[1])
    mul.setMultipole(3, refs[2])
    mul.setMultipole(4, refs[3])

    assert mul.getMultipole(1).real == pytest.approx(refs[0])
    assert mul.getMultipole(2).real == pytest.approx(refs[1])
    assert mul.getMultipole(3).real == pytest.approx(refs[2])
    assert mul.getMultipole(4).real == pytest.approx(refs[3])

    factors = np.array([-11 / 13.0, 17.0 / 19.0, 23 / 29.0, 31.0 / 37.0])
    n_refs = refs * factors
    mul *= factors

    assert mul.getMultipole(1).real == pytest.approx(n_refs[0])
    assert mul.getMultipole(2).real == pytest.approx(n_refs[1])
    assert mul.getMultipole(3).real == pytest.approx(n_refs[2])
    assert mul.getMultipole(4).real == pytest.approx(n_refs[3])


def test21_multiply_multipoles_scalar_vector():
    h_max = 4
    mul = tslib.TwoDimensionalMultipoles(h_max)
    refs = np.array([2, 3, 5, 7])
    mul.setMultipole(1, refs[0])
    mul.setMultipole(2, refs[1])
    mul.setMultipole(3, refs[2])
    mul.setMultipole(4, refs[3])

    assert mul.getMultipole(1).real == pytest.approx(refs[0])
    assert mul.getMultipole(2).real == pytest.approx(refs[1])
    assert mul.getMultipole(3).real == pytest.approx(refs[2])
    assert mul.getMultipole(4).real == pytest.approx(refs[3])

    factors = np.array([-11 / 13.0, 17.0 / 19.0, 23 / 29.0, 31.0 / 37.0])
    n_refs = refs * factors
    n_mul = mul * factors

    assert mul.getMultipole(1).real == pytest.approx(refs[0])
    assert mul.getMultipole(2).real == pytest.approx(refs[1])
    assert mul.getMultipole(3).real == pytest.approx(refs[2])
    assert mul.getMultipole(4).real == pytest.approx(refs[3])

    assert n_mul.getMultipole(1).real == pytest.approx(n_refs[0])
    assert n_mul.getMultipole(2).real == pytest.approx(n_refs[1])
    assert n_mul.getMultipole(3).real == pytest.approx(n_refs[2])
    assert n_mul.getMultipole(4).real == pytest.approx(n_refs[3])


def test31_multipoles_element_access():
    """

    Todo:
       resolve the access issue...
       Find
    """
    h_max = 4
    mul = tslib.TwoDimensionalMultipoles(h_max)
    print(repr(mul))
    refs = np.array([2, 3, 5, 7])
    ## tmp = np.array(mul, copy=False)
    ## tmp[:] = refs
    ## print(tmp)
    ## print(repr(mul))
    ##
    ## assert mul.getMultipole(1).real == pytest.approx(refs[0])
    ## assert mul.getMultipole(2).real == pytest.approx(refs[1])
    ## assert mul.getMultipole(3).real == pytest.approx(refs[2])
    ## assert mul.getMultipole(4).real == pytest.approx(refs[3])
