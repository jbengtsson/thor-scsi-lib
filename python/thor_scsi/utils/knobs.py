from .. import lib as tslib
import gtpsa

corresponding_types = {
    tslib.Quadrupole: tslib.QuadrupoleTpsa,
    tslib.HorizontalSteerer: tslib.HorizontalSteererTpsa,
    tslib.VerticalSteerer: tslib.VerticalSteererTpsa,
}


def convert_magnet_to_knobbable(a_magnet: tslib.Mpole) -> tslib.MpoleTpsa:
    config = a_magnet.config()
    corresponding_type = corresponding_types[type(a_magnet)]
    return corresponding_type(config)


def make_magnet_offset_knobbable(
    magnet: tslib.FieldKick,
    *,
    po: int,
    desc: gtpsa.desc,
    named_index: gtpsa.IndexMapping
):
    """

    Args:
        magnet: magnet to make knobbable
        po: maximum order of the parameter
        desc:
        named_index:

    Returns:

    """
    dxv = magnet.get_dx().to_object()
    dyv = magnet.get_dy().to_object()
    dx = gtpsa.tpsa(desc, po, mapping=named_index)
    dy = gtpsa.tpsa(desc, po, mapping=named_index)
    dx.set_knob(dxv, "dx")
    dy.set_knob(dyv, "dy")
    dx.name = magnet.name + "_dx"
    dy.name = magnet.name + "_dy"
    magnet.set_dx(gtpsa.TpsaOrDouble(dx))
    magnet.set_dy(gtpsa.TpsaOrDouble(dy))


def make_magnet_offset_unknobbable(magnet: tslib.FieldKick):
    dx = float(magnet.get_dx().to_object().get())
    dy = float(magnet.get_dy().to_object().get())
    magnet.set_dx(gtpsa.TpsaOrDouble(dx))
    magnet.set_dy(gtpsa.TpsaOrDouble(dy))


def make_magnet_strength_knobbable(
    magnet: tslib.Mpole,
    *,
    po: int,
    multipole_number: int = None,
    desc: gtpsa.desc,
    named_index: gtpsa.IndexMapping
):
    if multipole_number is None:
        multipole_number = magnet.get_main_multipole_number()
    k_orig = magnet.get_main_multipole_strength().to_object()
    k = gtpsa.ctpsa(desc, po, mapping=named_index)
    k.name = magnet.name + "_K"
    k.set_knob(k_orig, "K")
    magnet.get_multipoles().set_multipole(multipole_number, gtpsa.CTpsaOrComplex(k))


def make_magnet_strength_unknobbable(magnet: tslib.Mpole, multipole_number: int = None):
    if multipole_number is None:
        multipole_number = magnet.get_main_multipole_number()
    mul = complex(
        magnet.get_multipoles().get_multipole(multipole_number).to_object().get()
    )
    magnet.get_multipoles().set_multipole(multipole_number, gtpsa.CTpsaOrComplex(mul))


def make_magnet_knobbable(
    magnet, *, po: int, desc, named_index, multipole_number: int = None, offset=False
):
    make_magnet_strength_knobbable(
        magnet,
        po=po,
        multipole_number=multipole_number,
        desc=desc,
        named_index=named_index,
    )
    if offset:
        make_magnet_offset_knobbable(magnet, po=po, desc=desc, named_index=named_index)
    return magnet


def make_magnet_unknobbable(magnet, multipole_number: int = None, offset: bool = False):
    """Replace knobbed variables with standard types"""
    if offset:
        make_magnet_offset_unknobbable(magnet)
    make_magnet_strength_unknobbable(magnet, multipole_number=multipole_number)


def convert_if_steerer(elem):
    if isinstance(elem, (tslib.HorizontalSteerer, tslib.VerticalSteerer)):
        n_elem = convert_magnet_to_knobbable(elem)
    else:
        n_elem = elem
    return n_elem


def convert_if_quadrupole(elem):
    if isinstance(elem, tslib.Quadrupole):
        n_elem = convert_magnet_to_knobbable(elem)
    else:
        n_elem = elem
    return n_elem
