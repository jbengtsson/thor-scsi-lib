#include <pybind11/stl.h>
// #include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <tps/enums.h>
#include <string>
#include <sstream>

#include "thor_scsi.h"

namespace py = pybind11;

void py_thor_scsi_init_tps(py::module_ &m)
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


	declare_field<double>(m, "_double");
	declare_field<tps>(m, "_tps");

	//py::class_<ss_vect<tps>>(m, "ss_vect_tps")
//		.def(py::init<>());

	//py::class_<ss_vect<double>>(m, "ss_vect_double")
	//	.def(py::init<>());


}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
