#include <pybind11/stl.h>
// #include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include "thor_scsi.h"
#include <thor_scsi/core/machine.h>
#include <thor_scsi/std_machine/std_machine.h>

//namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace py = pybind11;

void py_thor_scsi_init_machine(py::module &m)
{

	// m.def("register_elements", &register_elements);
	// needs to be done only once
	register_elements();

	py::class_<tsc::Machine, std::shared_ptr<tsc::Machine>>(m, "Accelerator")
		.def("find", &tsc::Machine::find)
		.def("propagate", py::overload_cast<tsc::ConfigType&, tsc::ss_vect_dbl&, size_t, int>(&tsc::Machine::propagate), "xxx")
		.def("propagate", py::overload_cast<tsc::ConfigType&, tsc::ss_vect_tps&, size_t, int>   (&tsc::Machine::propagate),  "xxx")

		.def("__len__", &tsc::Machine::size)
		.def("__getitem__", py::overload_cast<size_t>(&tsc::Machine::at))
		.def(py::init<const Config &>());



}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
