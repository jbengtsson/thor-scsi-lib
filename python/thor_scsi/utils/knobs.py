from .. import lib as ts
import gtpsa

corresponding_types = {
    ts.Bending:           ts.BendingTpsa,
    ts.Sextupole:         ts.SextupoleTpsa,
    ts.Quadrupole:        ts.QuadrupoleTpsa,
    ts.HorizontalSteerer: ts.HorizontalSteererTpsa,
    ts.VerticalSteerer:   ts.VerticalSteererTpsa,
}


def convert_magnet_to_knobbable(a_magnet: ts.Mpole) -> ts.MpoleTpsa:
    config = a_magnet.config()
    corresponding_type = corresponding_types[type(a_magnet)]
    return corresponding_type(config)


def make_magnet_offset_knobbable(
        magnet: ts.FieldKick,
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
    # dx.set_variable(dxv, "dx")
    # dy.set_variable(dyv, "dy")
    dx.name = magnet.name + "_dx"
    dy.name = magnet.name + "_dy"
    magnet.set_dx(gtpsa.TpsaOrDouble(dx))
    magnet.set_dy(gtpsa.TpsaOrDouble(dy))


def make_magnet_offset_unknobbable(magnet: ts.FieldKick):
    dx = float(magnet.get_dx().to_object().get())
    dy = float(magnet.get_dy().to_object().get())
    magnet.set_dx(gtpsa.TpsaOrDouble(dx))
    magnet.set_dy(gtpsa.TpsaOrDouble(dy))


def make_magnet_strength_knobbable(
        magnet: ts.Mpole,
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
    # k.set_variable(k_orig, "K")
    magnet.get_multipoles().set_multipole(
        multipole_number, gtpsa.CTpsaOrComplex(k))


def make_magnet_strength_unknobbable(
        magnet: ts.Mpole, multipole_number: int = None):
    if multipole_number is None:
        multipole_number = magnet.get_main_multipole_number()
        mul = complex(
            magnet.get_multipoles().get_multipole(multipole_number).
            to_object().get()
        )
        magnet.get_multipoles().set_multipole(
            multipole_number, gtpsa.CTpsaOrComplex(mul))


def make_magnet_knobbable(
        magnet, *, po: int, desc, named_index, multipole_number: int = None,
        offset=False
):
    make_magnet_strength_knobbable(
        magnet,
        po=po,
        multipole_number=multipole_number,
        desc=desc,
        named_index=named_index,
    )
    if offset:
        make_magnet_offset_knobbable(
            magnet, po=po, desc=desc, named_index=named_index)
    return magnet


def make_magnet_unknobbable(
        magnet, multipole_number: int = None, offset: bool = False):
    """Replace knobbed variables with standard types"""
    if offset:
        make_magnet_offset_unknobbable(magnet)
    make_magnet_strength_unknobbable(magnet, multipole_number=multipole_number)


def convert_if_steerer(elem):
    if isinstance(elem, (ts.HorizontalSteerer, ts.VerticalSteerer)):
        n_elem = convert_magnet_to_knobbable(elem)
    else:
        n_elem = elem
    return n_elem


def convert_if_quadrupole(elem):
    if isinstance(elem, ts.Quadrupole):
        n_elem = convert_magnet_to_knobbable(elem)
    else:
        n_elem = elem
    return n_elem


__all__ = [
    "convert_magnet_to_knobbable", "", "make_magnet_offset_knobbable",
    "make_magnet_offset_unknobbable", "make_magnet_strength_knobbable",
    "make_magnet_strength_unknobbable", "make_magnet_knobbable",
    "make_magnet_unknobbable", "convert_if_steerer", "convert_if_quadrupole"
]
