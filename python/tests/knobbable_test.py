import pytest
import gtpsa
import gtpsa.utils_df
from thor_scsi.utils import knobs
import thor_scsi.lib as tslib
from thor_scsi.pyflame import Config


def test_knob_quadrupole():
    length = 0.2
    K = 1.2
    N = 4
    name = "qd"
    config = Config()
    C = Config()
    C.setAny("L", length)
    C.setAny("name", name)
    C.setAny("K", K)
    C.setAny("N", N)
    quad = tslib.Quadrupole(C)
    nquad = knobs.convert_magnet_to_knobbable(quad)
    assert nquad.name == quad.name
    assert nquad.get_length() == pytest.approx(quad.get_length(), rel=1e-12)
    # test of to_object should be avoided
    K_check = nquad.get_main_multipole_strength().to_object()
    assert K_check == pytest.approx(K, rel=1e-12)


def test_knob_horizontal_steerer():
    length = 0.0
    K = 0.1
    N = 1
    name = "h_st"
    config = Config()
    C = Config()
    C.setAny("L", length)
    C.setAny("name", name)
    C.setAny("K", K)
    C.setAny("N", N)
    h_st = tslib.HorizontalSteerer(C)
    nh_st = knobs.convert_magnet_to_knobbable(h_st)
    assert h_st.name == nh_st.name
    assert h_st.get_length() == pytest.approx(length, rel=1e-12)

    named_index_d = dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6)
    named_index = gtpsa.IndexMapping(named_index_d)
    desc = gtpsa.desc(6, 2, 1, 2)
    knobs.make_magnet_knobbable(nh_st, po=1, desc=desc, named_index=named_index)


def test_knob_vertical_steerer():
    length = 0.0
    K = 0.2
    N = 1
    name = "v_st"
    config = Config()
    C = Config()
    C.setAny("L", length)
    C.setAny("name", name)
    C.setAny("K", K)
    C.setAny("N", N)
    v_st = tslib.VerticalSteerer(C)
    nv_st = knobs.convert_magnet_to_knobbable(v_st)
    assert nv_st.name == v_st.name
    assert nv_st.get_length() == pytest.approx(length, rel=1e-12)


def test_quadrupole_dependence_on_offset():
    """check that the parameter dependence on dx works as expected

    Todo:
        Inspect why number of integration steps does not increase
        prediction accuracy
    """
    length = 0.25
    K = 1
    N = 1000
    name = "a_quad"
    config = Config()
    C = Config()
    C.setAny("L", length)
    C.setAny("name", name)
    C.setAny("K", K)
    C.setAny("N", N)
    quad = tslib.QuadrupoleTpsa(C)
    quad.set_number_of_integration_steps(N)
    assert quad.get_number_of_integration_steps() == N

    named_index_d = dict(x=0, px=1, y=2, py=3, delta=4, ct=5, K=6, dx=7, dy=8)
    named_index = gtpsa.IndexMapping(named_index_d)
    desc = gtpsa.desc(6, 2, 3, 2)
    knobs.make_magnet_knobbable(
        quad, po=2, offset=True, desc=desc, named_index=named_index
    )

    calc_config = tslib.ConfigType()
    ps = gtpsa.ss_vect_tpsa(desc, 2, 6, named_index)
    ps.set_identity()
    quad.propagate(calc_config, ps)
    # offset of 1 (unit =m)
    expected_px_vs_dx =  K * (1) * length
    expected_x_vs_dx = expected_px_vs_dx * 1/2 * length
    print("\nTpsa estimate: x  ", ps.x.get(dx=1),  "expected", expected_x_vs_dx)
    print("Tpsa estimate: px ", ps.px.get(dx=1), "expected", expected_px_vs_dx)

    # assert ps.x.get(dx=1) == pytest.approx(expected_dx_dependence, rel=1e-12)

    # Estimate effect based on propagation
    # consistency check
    quad = tslib.Quadrupole(C)
    offset = 1e-3
    quad.set_dx(offset)
    ps = gtpsa.ss_vect_double(0e0, index_mapping=named_index)
    ps.set_zero()
    quad.propagate(calc_config, ps)
    print(ps.x, ps.px)
    assert ps.x == pytest.approx(expected_x_vs_dx * offset, rel=1e-2)
    assert ps.px == pytest.approx(expected_px_vs_dx * offset, rel=1e-1)
