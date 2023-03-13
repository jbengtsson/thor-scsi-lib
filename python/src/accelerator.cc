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


static const char acc_init_conf_doc[] = \
"initalise accelerator with a configuration file\n\
\n\
Args:\n\
   conf: configuration for the accelerator (aka parsed lattise file)\n\
   add_marker_at_start: add a marker at the start of the lattice\n\
\n\
Warning:\n\
  consider if the marker is not better added manually to the lattice\n\
  file\n";

static const char acc_init_list_doc[] = \
"initalise accelerator with a list of elements\n\
\n\
Args:\n\
   elements:            list of elements\n\
   add_marker_at_start: add a marker at the start of the lattice\n\
\n\
Warning:\n\
  consider if the marker is not better added manually to the lattice\n\
  file";

static const char prop_doc[] = "propagate phase space through elements";

template<typename Types, typename Class>
void add_methods_accelerator(py::class_<Class> t_acc)
{
	// using double_type = typename Types::double_type;

	t_acc
		.def("find",                 &Class::find)
		// unique_ptr not working ... check for memory management ...
		// should now be fixed with shared_ptrs
		.def("elements_with_name",     &Class::elementsWithName)
		// unique_ptr not working ... check for memory management ...
		.def("elements_with_name_type", &Class::elementsWithNameType)
		.def("set_logger",            &Class::set_logger)
		.def("set_log_level",          [](Class& t_acc, int level){
			t_acc.set_log_level(level);
		}, "set the log level, have a look to the accelerator log level enum")
#if 0
		.def("set_trace",             [](Class& t_acc, const bool flag){
			// Not used as
			py::scoped_ostream_redirect stream(
				std::cout,                               // std::ostream&
				py::module_::import("sys").attr("stdout") // Python output
				);
			t_acc.set_trace(std::cout);
		},  "activate trace output to std oute"
			)
#endif
		//.def("__copy__",             &Class::clone, "make a copy of the accelerator")
		.def("__len__",              &Class::size)
		.def("__getitem__", py::overload_cast<size_t>(&Class::at))
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_dbl&, size_t, int, size_t, bool>(&Class::propagate), prop_doc,
		     py::arg("calc_config"), py::arg("ps"), py::arg("start") = 0, py::arg("max_elements") = imax, py::arg("n_turns") = n_turns,
		     py::arg("tracy_compatible") = false)
		/*
		  .def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_tps&, size_t, int, size_t, bool>(&Class::propagate), prop_doc,
		     py::arg("calc_config"), py::arg("ps"), py::arg("start") = 0, py::arg("max_elements") = imax, py::arg("n_turns") = n_turns,
		     py::arg("tracy_compatible") = false)
		*/
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_tpsa&, size_t, int, size_t, bool>(&Class::propagate), prop_doc,
		     py::arg("calc_config"), py::arg("ps"), py::arg("start") = 0, py::arg("max_elements") = imax, py::arg("n_turns") = n_turns,
		     py::arg("tracy_compatible") = false)
		.def(py::init<const Config &, bool>(), acc_init_list_doc,
		     py::arg("config object"), py::arg("add_marker_at_start") = false)

		.def(py::init<std::vector<std::shared_ptr<thor_scsi::core::CellVoid>>&, bool>(), acc_init_list_doc,
		     py::arg("elements"), py::arg("add_marker_at_start") = false)

		;
}

void py_thor_scsi_init_accelerator(py::module &m)
{

	// m.def("register_elements", &register_elements);
	// needs to be done only once
	ts::register_elements();

	py::class_<tsc::Machine::LogRecord> log_record(m, "LogRecord");
	py::class_<tsc::Machine::Logger, std::shared_ptr<tsc::Machine::Logger>> logger(m, "Logger");


	py::enum_<MachineLogLevel>(m, "accelerator_log_level")
		.value("error",   MachineLogLevel::error)
		.value("warn",    MachineLogLevel::warn)
		.value("warning", MachineLogLevel::warning)
		.value("info",    MachineLogLevel::info)
		.value("debug",   MachineLogLevel::debug)
		.value("fine",    MachineLogLevel::fine);



	py::class_<ts::Accelerator, std::shared_ptr<ts::Accelerator>> acc(m, "Accelerator");
	add_methods_accelerator<tsc::StandardDoubleType, ts::Accelerator>(acc);

	py::class_<ts::AcceleratorTpsa, std::shared_ptr<ts::AcceleratorTpsa>> accK(m, "AcceleratorTpsa");
	  add_methods_accelerator<tsc::TpsaVariantType, ts::AcceleratorTpsa>(accK);


}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
