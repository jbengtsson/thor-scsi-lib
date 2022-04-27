#include <pybind11/stl.h>
// #include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include "thor_scsi.h"
#include <thor_scsi/core/machine.h>
#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/std_machine/accelerator.h>

//namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace ts = thor_scsi;
namespace py = pybind11;


void py_thor_scsi_init_accelerator(py::module &m)
{

	// m.def("register_elements", &register_elements);
	// needs to be done only once
	ts::register_elements();

	py::class_<tsc::Machine::LogRecord> log_record(m, "LogRecord");
	py::class_<tsc::Machine::Logger, std::shared_ptr<tsc::Machine::Logger>> logger(m, "Logger");

	const char prop_doc[] = "propagate phase space through elements";
	const int imax = std::numeric_limits<int>::max();

	py::class_<ts::Accelerator, std::shared_ptr<ts::Accelerator>>(m, "Accelerator")
		.def("find",                 &ts::Accelerator::find)
		// unique_ptr not working ... check for memory management ...
		// should now be fixed with shared_ptrs
		.def("elementsWithName",     &ts::Accelerator::elementsWithName)
		// unique_ptr not working ... check for memory management ...
		.def("elementsWithNameType", &ts::Accelerator::elementsWithNameType)
		.def("setLogger",            &ts::Accelerator::set_logger)
		.def("setTrace",             &ts::Accelerator::set_trace, "register the stream to log to")
		.def("__len__",              &ts::Accelerator::size)
		.def("__getitem__", py::overload_cast<size_t>(&ts::Accelerator::at))
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_dbl&, size_t, int>(&ts::Accelerator::propagate), prop_doc,
		     py::arg("calc_config"), py::arg("ps"), py::arg("start") = 0, py::arg("max_elements") = imax)
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_tps&, size_t, int>(&ts::Accelerator::propagate), prop_doc,
		     py::arg("calc_config"), py::arg("ps"), py::arg("start") = 0, py::arg("max_elements") = imax)
		.def(py::init<const Config &>());



}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
