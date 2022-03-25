#include <pybind11/stl.h>
// #include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <tps/tpsa_lin.h>
#include <tps/enums.h>
#include <string>
#include <sstream>

#include "thor_scsi.h"

namespace py = pybind11;

void py_thor_scsi_init_tps(py::module &m)
{

	py::enum_<spatial_ind>(m, "spatial_ind")
		.value("X_", spatial_ind::X_)
		.value("Y_", spatial_ind::Y_)
		.value("Z_", spatial_ind::Z_);

	py::enum_<phase_space_ind>(m, "phase_space_ind")
		.value("x_",     phase_space_ind::x_)
		.value("px_",    phase_space_ind::px_)
		.value("y_",     phase_space_ind::y_)
		.value("py_",    phase_space_ind::py_)
		.value("delta_", phase_space_ind::delta_)
		.value("ct_",    phase_space_ind::ct_);


	// Polymorphic number class.
	// Classes.

	py::class_<tps>(m, "tps")
		.def(py::init<>())
		.def(py::init<const double>())
		.def(py::init<const double, const int>())
		.def("__str__", &tps::pstr)
		.def("__repr__", &tps::repr);


	// declare_field<double>(m, "_double");
	// declare_field<tps>(m, "_tps");

	py::class_<ss_vect<tps>>(m, "ss_vect_tps")
		.def("set_identity", &ss_vect<tps>::set_identity)
		.def("set_zero",     &ss_vect<tps>::set_zero)
		.def("__repr__",     &ss_vect<tps>::repr)
		.def("__str__",      &ss_vect<tps>::pstr)
		.def(py::init<>());

	py::class_<ss_vect<double>>(m, "ss_vect_double")
		//.def("set_identity", &ss_vect<double>::set_identity)
		.def("set_zero", &ss_vect<double>::set_zero)
		.def("__repr__",  &ss_vect<double>::repr)
		.def("__str__",  &ss_vect<double>::pstr)
		.def(py::init<>());

	m.def("lists_to_ss_vect_tps", &stlmattomap_save);
	m.def("ss_vect_tps_to_lists", &maptostlmat);


}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
