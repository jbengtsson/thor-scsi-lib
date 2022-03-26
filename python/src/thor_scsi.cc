#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
// #include <pybind11/operators.h>
#include <vector>
#include <tuple>

#include <thor_scsi/version.h>
#include <tps/enums.h>
// #include <thor_scsi/core/defines.h>
#include <thor_scsi/elements/constants.h>
#include <thor_scsi/core/cells.h>
// #include <thor_scsi/core/elements_enums.h>
#include "thor_scsi.h"
#include <string>
#include <sstream>

namespace py = pybind11;
namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;


// Polymorphic number template class wrapper.



/*
 * get_eps_x: provide the memory and return it to python
 *            TODO: check for memory leaks!
 */


PYBIND11_MODULE(lib, scsi) {
    scsi.doc() = "Self-Consistent Symplectic Integrator (SCSI)";

    py_thor_scsi_init_tps(scsi);
    py_thor_scsi_init_elements(scsi);
    py_thor_scsi_init_accelerator(scsi);
    py_thor_scsi_init_config_type(scsi);
    // py_thor_scsi_init_lattice(scsi);


    // Beamline class.

    // Constants.
    //scsi.attr("HOMmax") = HOMmax;
    scsi.attr("c0")     = tse::speed_of_light;
    // scsi.attr("q_e")    = tsc::q_e;
    // scsi.attr("m_e")    = tsc::m_e;
    // scsi.attr("mu_0")   = tsc::mu_0;
    // scsi.attr("eps_0")  = tsc::eps_0;
    // scsi.attr("r_e")    = tsc::r_e;
    // scsi.attr("h_bar")  = tsc::h_bar;
    /*
    */
    // Enums.

    #if 0
    py::enum_<tse::PartsKind>(scsi, "PartsKind")
      .value("drift",      tse::PartsKind::drift)
      .value("Wigl",       tse::PartsKind::Wigl)
      .value("Mpole",      tse::PartsKind::Mpole)
      .value("Cavity",     tse::PartsKind::Cavity)
      .value("marker",     tse::PartsKind::marker)
      .value("undef",      tse::PartsKind::undef)
      .value("Insertion",  tse::PartsKind::Insertion)
      .value("FieldMap",   tse::PartsKind::FieldMap)
      .value("Spreader",   tse::PartsKind::Spreader)
      .value("Recombiner", tse::PartsKind::Recombiner)
      .value("Solenoid",   tse::PartsKind::Solenoid)
      .value("Map",        tse::PartsKind::Map)
      .export_values();
    #endif

    #if 0
    py::enum_<tse::MpoleKind>(scsi, "MpoleKind")
      .value("All",   tse::MpoleKind::All)
      .value("Dip",   tse::MpoleKind::Dip)
      .value("Quad",  tse::MpoleKind::Quad)
      .value("Sext",  tse::MpoleKind::Sext)
      .value("Oct",   tse::MpoleKind::Oct)
      .value("Dec",   tse::MpoleKind::Dec)
      .value("Dodec", tse::MpoleKind::Dodec);
    #endif

#if 0
    py::enum_<tse::PlaneKind>(scsi, "Horizontal")
      .value("Horizontal", tse::PlaneKind::Horizontal)
      .value("Vertical",   tse::PlaneKind::Vertical);
#endif

    // Classes.

#if 0
    py::class_<tse::CellType>(scsi, "CellType")
      // int
      .def_readwrite("Fnum",       &tse::CellType::Fnum)
      .def_readwrite("Knum",       &tse::CellType::Knum)
      // double
      .def_readwrite("S",          &tse::CellType::S)
      .def_readwrite("curly_dH_x", &tse::CellType::curly_dH_x)
      // std::vector<double>:
      .def_readwrite("dI",         &tse::CellType::dI)
      .def_readwrite("dS",         &tse::CellType::dS)
      .def_readwrite("dT",         &tse::CellType::dT)
      .def_readwrite("Eta",        &tse::CellType::Eta)
      .def_readwrite("Etap",       &tse::CellType::Etap)
      .def_readwrite("Alpha",      &tse::CellType::Alpha)
      .def_readwrite("Beta",       &tse::CellType::Beta)
      .def_readwrite("Nu",         &tse::CellType::Nu)
      .def_readwrite("BeamPos",    &tse::CellType::BeamPos)
      // std::vector< std::vector<double> >
      .def_readwrite("maxampl",    &tse::CellType::maxampl)
      .def_readwrite("A",          &tse::CellType::A)
      .def_readwrite("sigma",      &tse::CellType::sigma)
      // tse::CellType
      .def(py::init<>());
#endif

}
