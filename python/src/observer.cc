#include <pybind11/pybind11.h>
#include <thor_scsi/elements/standard_observer.h>
#include <gtpsa/python/objects_with_named_index.h>

namespace py = pybind11;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;
namespace gpy = gtpsa::python;

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

namespace thor_scsi::python {
    class StandardObserverWithIndex : public tse::StandardObserver {
	std::shared_ptr<gpy::IndexMapping> m_mapping = nullptr;
	using base = tse::StandardObserver;

    public:
	StandardObserverWithIndex(std::shared_ptr<gpy::IndexMapping> mapping = gpy::default_index_mapping_ptr)
	    : m_mapping(mapping)
	    {}

	auto getMapping(void) const {
	    return std::const_pointer_cast<gpy::IndexMapping>(this->m_mapping);
	}

	auto getTruncatedPowerSeriesA(void) {
	    // can I avoid this copy?
	    gtpsa::ss_vect<gtpsa::tpsa> *vec = base::getTruncatedPowerSeriesA().get();
	    return gpy::StateSpaceWithNamedIndex<gtpsa::tpsa>(*vec, this->getMapping());
	}

	auto getPhaseSpace(void) {
	    // can I avoid this copy?
	    gtpsa::ss_vect<double> vec = base::getPhaseSpace();
	    return gpy::StateSpaceWithNamedIndex<double>(vec, this->getMapping());
	}
    };
}
namespace tpy = thor_scsi::python;

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

	py::class_<tse::StandardObserver, tsc::Observer, std::shared_ptr<tse::StandardObserver>> std_observer_intern(m, "_StandardObserver");
	std_observer_intern
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


	py::class_<tpy::StandardObserverWithIndex, std::shared_ptr<tpy::StandardObserverWithIndex>> std_observer(m, "StandardObserver", std_observer_intern);
	std_observer
	    .def("get_phase_space",              &tpy::StandardObserverWithIndex::getPhaseSpace)
	    .def("get_truncated_power_series_a", &tpy::StandardObserverWithIndex::getTruncatedPowerSeriesA)
	    .def(py::init<std::shared_ptr<gpy::IndexMapping>>())
	    ;

}
