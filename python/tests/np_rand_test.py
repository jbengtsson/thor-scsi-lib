import pytest
from functools import partial

from numpy.random import Generator, PCG64

gen = Generator(PCG64(seed=355113))


def test00_reproducible():
    """Test that the generator returns the same seed every time

    assuming that the same seed is used
    """
    test = gen.integers(low=-2**30, high=2**30)
    assert test == 352929011

def test01_reproducible_second():
    """Check that second run gives the same result
    """
    test = gen.integers(low=-2**30, high=2**30)
    assert test == -13297206


def test10_reproducible_uniform():
    """Check that uniform returns the same value

    """
    test = gen.uniform(low=2, high=3)
    assert test >= 2
    assert test <= 3
    assert pytest.approx(test, 2.4252)
