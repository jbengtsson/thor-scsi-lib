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

static const size_t n_turns=1;
static const int imax = std::numeric_limits<int>::max();

enum MachineLogLevel {
	error = THOR_SCSI_ERROR,
	warn = THOR_SCSI_WARN,
// compatible with python
	warning = THOR_SCSI_WARNING,
	info = THOR_SCSI_INFO,
	debug = THOR_SCSI_DEBUG,
// Super verbose logging, a la the rf cavity element,
	fine = THOR_SCSI_FINE
};

void py_thor_scsi_init_accelerator(py::module &m)
{

	// m.def("register_elements", &register_elements);
	// needs to be done only once
	ts::register_elements();

	py::class_<tsc::Machine::LogRecord> log_record(m, "LogRecord");
	py::class_<tsc::Machine::Logger, std::shared_ptr<tsc::Machine::Logger>> logger(m, "Logger");

	const char prop_doc[] = "propagate phase space through elements";

    py::enum_<MachineLogLevel>(m, "accelerator_log_level")
	    .value("error",   MachineLogLevel::error)
	    .value("warn",    MachineLogLevel::warn)
	    .value("warning", MachineLogLevel::warning)
	    .value("info",    MachineLogLevel::info)
	    .value("debug",   MachineLogLevel::debug)
	    .value("fine",    MachineLogLevel::fine);


	const char acc_init_conf_doc[] = \
"initalise accelerator with a configuration file\n\
\n\
Args:\n\
   conf: configuration for the accelerator (aka parsed lattise file)\n\
   add_marker_at_start: add a marker at the start of the lattice\n\
\n\
Warning:\n\
  consider if the marker is not better added manually to the lattice\n\
  file\n";

	const char acc_init_list_doc[] = \
"initalise accelerator with a list of elements\n\
\n\
Args:\n\
   elements:            list of elements\n\
   add_marker_at_start: add a marker at the start of the lattice\n\
\n\
Warning:\n\
  consider if the marker is not better added manually to the lattice\n\
  file";

	py::class_<ts::Accelerator, std::shared_ptr<ts::Accelerator>>(m, "Accelerator")
		.def("find",                 &ts::Accelerator::find)
		// unique_ptr not working ... check for memory management ...
		// should now be fixed with shared_ptrs
		.def("elements_with_name",     &ts::Accelerator::elementsWithName)
		// unique_ptr not working ... check for memory management ...
		.def("elements_with_name_type", &ts::Accelerator::elementsWithNameType)
		.def("set_logger",            &ts::Accelerator::set_logger)
		.def("set_log_level",          [](ts::Accelerator& acc, int level){
						       acc.set_log_level(level);
					       }, "set the log level, have a look to the accelerator log level enum")
#if 0
		.def("set_trace",             [](ts::Accelerator& acc, const bool flag){
						      // Not used as
						      py::scoped_ostream_redirect stream(
							      std::cout,                               // std::ostream&
							      py::module_::import("sys").attr("stdout") // Python output
							      );
						      acc.set_trace(std::cout);
					      },  "activate trace output to std oute"
			)
#endif
		//.def("__copy__",             &ts::Accelerator::clone, "make a copy of the accelerator")
		.def("__len__",              &ts::Accelerator::size)
		.def("__getitem__", py::overload_cast<size_t>(&ts::Accelerator::at))
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_dbl&, size_t, int, size_t, bool>(&ts::Accelerator::propagate), prop_doc,
		     py::arg("calc_config"), py::arg("ps"), py::arg("start") = 0, py::arg("max_elements") = imax, py::arg("n_turns") = n_turns,
		     py::arg("tracy_compatible") = false)
             /*
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_tps&, size_t, int, size_t, bool>(&ts::Accelerator::propagate), prop_doc,
		     py::arg("calc_config"), py::arg("ps"), py::arg("start") = 0, py::arg("max_elements") = imax, py::arg("n_turns") = n_turns,
		     py::arg("tracy_compatible") = false)
              */
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_tpsa&, size_t, int, size_t, bool>(&ts::Accelerator::propagate), prop_doc,
		     py::arg("calc_config"), py::arg("ps"), py::arg("start") = 0, py::arg("max_elements") = imax, py::arg("n_turns") = n_turns,
		     py::arg("tracy_compatible") = false)
		.def(py::init<const Config &, bool>(), acc_init_list_doc,
		     py::arg("config object"), py::arg("add_marker_at_start") = false)
		.def(py::init<const std::vector<std::shared_ptr<thor_scsi::core::ElemType>>&, bool>(), acc_init_list_doc,
		     py::arg("elements"), py::arg("add_marker_at_start") = false)
		;



}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
