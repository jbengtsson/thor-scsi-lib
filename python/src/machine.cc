#include <pybind11/stl.h>
// #include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include "thor_scsi.h"
#include <thor_scsi/core/machine.h>
#include <thor_scsi/std_machine/std_machine.h>

//namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace py = pybind11;

void py_thor_scsi_init_machine(py::module_ &m)
{

	// m.def("register_elements", &register_elements);
	// needs to be done only once
	register_elements();

	py::class_<tsc::Machine>(m, "Accelerator")
		.def("find", &tsc::Machine::find)
		//.def("__len__", &tsc::Machine::size)
		//.def("__getattr__", py::overload_cast<size_t>(&tsc::Machine::at))
		.def(py::init<const Config &>());


}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
