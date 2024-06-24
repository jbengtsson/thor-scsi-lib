# neeed to find out where to place it
# standard_canonocail_variable_dimensions?
import gtpsa
canonical_variable_default = dict(x=0, px=1, y=2, py=3, delta=4, ct=5)
canonical_variables_mapping_default = gtpsa.IndexMapping(canonical_variable_default)


__all__ = [
    "accelerator",
    "closed_orbit",
    "linear_optics",
    "twiss_output",
    "radiate",
    "twiss_output",
    "canonical_variables_mapping_default",
    "canonical_variable_default"
]
