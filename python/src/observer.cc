#include <pybind11/pybind11.h>
#include <thor_scsi/elements/standard_observer.h>

namespace py = pybind11;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

class PyObserver: public tsc::Observer {
public:
	using tsc::Observer::Observer;

	void view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, view, elem, ps, os, cnt);
	}
	void view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, view, elem, ps, os, cnt);
	}
	void view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState os, const int cnt) override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, view, elem, ps, os, cnt);
	}

	void show(std::ostream& strm, int level) const override {
		PYBIND11_OVERRIDE_PURE(void, tsc::Observer, strm, level);
	}
};

void py_thor_scsi_init_observers(py::module &m)
{

	py::enum_<tsc::ObservedState>(m, "ObservedState")
		.value("event",     tsc::ObservedState::event)
		.value("start",     tsc::ObservedState::start)
		.value("end",       tsc::ObservedState::end)
		.value("undefined", tsc::ObservedState::undefined)
		.value("failure",   tsc::ObservedState::failure);


	py::class_<tsc::Observer,  PyObserver, std::shared_ptr<tsc::Observer>> observer(m, "Observer");
	observer.def(py::init<>());

	py::class_<tse::StandardObserver, tsc::Observer, std::shared_ptr<tse::StandardObserver>> std_observer(m, "StandardObserver");
	std_observer
		.def("__str__",                      &tse::StandardObserver::pstr)
		.def("__repr__",                     &tse::StandardObserver::repr)
		.def("get_observed_name",            &tse::StandardObserver::getObservedName)
		.def("get_observed_index",           &tse::StandardObserver::getObservedIndex)
		.def("has_phase_space",              &tse::StandardObserver::hasPhaseSpace)
		.def("get_phase_space",              &tse::StandardObserver::getPhaseSpace)
		.def("has_truncated_power_series",   &tse::StandardObserver::hasTruncatedPowerSeries)
		.def("get_truncated_power_series",   &tse::StandardObserver::getTruncatedPowerSeries)
		.def("has_truncated_power_series_a", &tse::StandardObserver::hasTruncatedPowerSeriesA)
		.def("get_truncated_power_series_a", &tse::StandardObserver::getTruncatedPowerSeriesA)
		.def("reset",                        &tse::StandardObserver::reset)
		.def(py::init<>());

}
