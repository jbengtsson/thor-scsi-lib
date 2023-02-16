import pytest
import thor_scsi.lib as tslib
import numpy as np
import gtpsa

def test00_instantiate_multipoles():
    n_max = 20
    mul = tslib.TwoDimensionalMultipoles(n_max)


def test10_multiply_multipoles():
    mul = tslib.TwoDimensionalMultipoles(0e0)
    mul.set_multipole(1, 2)
    mul.set_multipole(2, 4)

    assert mul.get_multipole(1).real == pytest.approx(2)
    assert mul.get_multipole(2).real == pytest.approx(4)
    print(mul)

    mul *= 2
    assert mul.get_multipole(1).real == pytest.approx(4)
    assert mul.get_multipole(2).real == pytest.approx(8)


def test11_multiply_multipoles_scalar_int():
    mul = tslib.TwoDimensionalMultipoles(0e0)
    mul.set_multipole(1, 2)
    mul.set_multipole(2, 4)

    assert mul.get_multipole(1).real == pytest.approx(2)
    assert mul.get_multipole(2).real == pytest.approx(4)
    print(mul)

    n_mul = mul * 2
    # Former multipoles still the same
    assert mul.get_multipole(1).real == pytest.approx(2)
    assert mul.get_multipole(2).real == pytest.approx(4)

    assert n_mul.get_multipole(1).real == pytest.approx(4)
    assert n_mul.get_multipole(2).real == pytest.approx(8)


def test11_multiply_multipoles_scalar_float():
    mul = tslib.TwoDimensionalMultipoles(0e0)
    mul.set_multipole(1, 2)
    mul.set_multipole(2, 4)

    assert mul.get_multipole(1).real == pytest.approx(2)
    assert mul.get_multipole(2).real == pytest.approx(4)
    print(mul)

    n_mul = mul * float(2.0)
    # Former multipoles still the same
    assert mul.get_multipole(1).real == pytest.approx(2)
    assert mul.get_multipole(2).real == pytest.approx(4)

    assert n_mul.get_multipole(1).real == pytest.approx(4)
    assert n_mul.get_multipole(2).real == pytest.approx(8)

def test12_multiply_multipoles_scalar_float():
    """float on the other side ...

    Commutative law needs to be explicitly implemented
    """
    mul = tslib.TwoDimensionalMultipoles(0e0)
    mul.set_multipole(1, 2)
    mul.set_multipole(2, 4)

    assert mul.get_multipole(1).real == pytest.approx(2)
    assert mul.get_multipole(2).real == pytest.approx(4)
    print(mul)

    n_mul = complex(2.0) * mul
    # Former multipoles still the same
    assert mul.get_multipole(1).real == pytest.approx(2)
    assert mul.get_multipole(2).real == pytest.approx(4)

    assert n_mul.get_multipole(1).real == pytest.approx(4)
    assert n_mul.get_multipole(2).real == pytest.approx(8)


def test20_multiply_multipoles_scalar_vector():
    h_max = 4
    mul = tslib.TwoDimensionalMultipoles(0e0, h_max)
    refs = np.array([2, 3, 5, 7])
    mul.set_multipole(1, refs[0])
    mul.set_multipole(2, refs[1])
    mul.set_multipole(3, refs[2])
    mul.set_multipole(4, refs[3])

    assert mul.get_multipole(1).real == pytest.approx(refs[0])
    assert mul.get_multipole(2).real == pytest.approx(refs[1])
    assert mul.get_multipole(3).real == pytest.approx(refs[2])
    assert mul.get_multipole(4).real == pytest.approx(refs[3])

    factors = np.array([-11 / 13.0, 17.0 / 19.0, 23 / 29.0, 31.0 / 37.0])
    n_refs = refs * factors
    mul *= factors

    assert mul.get_multipole(1).real == pytest.approx(n_refs[0])
    assert mul.get_multipole(2).real == pytest.approx(n_refs[1])
    assert mul.get_multipole(3).real == pytest.approx(n_refs[2])
    assert mul.get_multipole(4).real == pytest.approx(n_refs[3])


def test21_multiply_multipoles_scalar_vector():
    h_max = 4
    mul = tslib.TwoDimensionalMultipoles(0e0, h_max)
    refs = np.array([2, 3, 5, 7])
    mul.set_multipole(1, refs[0])
    mul.set_multipole(2, refs[1])
    mul.set_multipole(3, refs[2])
    mul.set_multipole(4, refs[3])

    assert mul.get_multipole(1).real == pytest.approx(refs[0])
    assert mul.get_multipole(2).real == pytest.approx(refs[1])
    assert mul.get_multipole(3).real == pytest.approx(refs[2])
    assert mul.get_multipole(4).real == pytest.approx(refs[3])

    factors = np.array([-11 / 13.0, 17.0 / 19.0, 23 / 29.0, 31.0 / 37.0])
    n_refs = refs * factors
    n_mul = mul * factors

    assert mul.get_multipole(1).real == pytest.approx(refs[0])
    assert mul.get_multipole(2).real == pytest.approx(refs[1])
    assert mul.get_multipole(3).real == pytest.approx(refs[2])
    assert mul.get_multipole(4).real == pytest.approx(refs[3])

    assert n_mul.get_multipole(1).real == pytest.approx(n_refs[0])
    assert n_mul.get_multipole(2).real == pytest.approx(n_refs[1])
    assert n_mul.get_multipole(3).real == pytest.approx(n_refs[2])
    assert n_mul.get_multipole(4).real == pytest.approx(n_refs[3])


def test31_multipoles_element_access():
    """

    Todo:
       resolve the access issue...
       Find
    """
    h_max = 4
    mul = tslib.TwoDimensionalMultipoles(0e0, h_max)
    print(repr(mul))
    refs = np.array([2, 3, 5, 7])
    ## tmp = np.array(mul, copy=False)
    ## tmp[:] = refs
    ## print(tmp)
    ## print(repr(mul))
    ##
    ## assert mul.get_multipole(1).real == pytest.approx(refs[0])
    ## assert mul.get_multipole(2).real == pytest.approx(refs[1])
    ## assert mul.get_multipole(3).real == pytest.approx(refs[2])
    ## assert mul.get_multipole(4).real == pytest.approx(refs[3])

def test40_variance_multipoles():
    import gtpsa
    import gtpsa._gtpsa_variant_test
    import sys
    mo = 10
    desc = gtpsa.desc(8, mo)


    h_max = 4
    mul = tslib.TwoDimensionalMultipolesTpsa(0e0, h_max)
    # Make the dipole and sextupole field dependent parameters
    df = gtpsa.ctpsa(desc, mo)
    df.set(0, 1+0.001j)
    # Dependence on variable .. convention after the state space
    df.setv(1, [0j] * 6 + [1+0j] + [0j])
    mul.set_multipole(1, gtpsa._gtpsa_variant_test.CTpsaOrComplex(df))

    qf = gtpsa.ctpsa(desc, mo)
    qf.set(0, 1 + 2j)
    # Dependence on variable .. convention after the state space
    qf.setv(1, [0j] * 7 + [1+0j])
    mul.set_multipole(3, gtpsa._gtpsa_variant_test.CTpsaOrComplex(qf))
    # print("multipoles", repr(mul))

    x = gtpsa.tpsa(desc, mo)
    x.set(0, 1e-3)
    x.setv(1, [1] + [0] * 5)
    y = gtpsa.tpsa(desc, mo)
    y.set(0, 2e-3)
    y.setv(1, [0, 0, 1] + [0] * 3 )
    z = gtpsa.ctpsa(x, y)
    zv = gtpsa._gtpsa_variant_test.CTpsaOrComplex(z)
    z.print("z", 0e0)
    (z**2).print("(z**2)", 1e-12)
    (zv*zv).to_object().print("(zv**2)", 1e-12)
    (zv*zv).to_object().print("(zv**2)", 1e-12)
    sys.stdout.flush()
    print("z cst", z.get())
    print("zv", zv)
    r = mul.field(zv).to_object()
    print("order", r.order)
    print(r.print("result"))
    assert(1==0)
