#include <pybind11/stl.h>
#include "thor_scsi.h"
#include <pybind11/pybind11.h>
#include <thor_scsi/core/config.h>

namespace tsc = thor_scsi::core;
namespace py = pybind11;

void py_thor_scsi_init_config_type(py::module &m)
{
  py::class_<tsc::ConfigType, std::shared_ptr<tsc::ConfigType>>(m, "ConfigType")
      // bool
      .def_readwrite("trace",          &tsc::ConfigType::trace)
      .def_readwrite("reverse_elem",   &tsc::ConfigType::reverse_elem)
      .def_readwrite("stable",         &tsc::ConfigType::stable)
      .def_readwrite("ErrFlag",        &tsc::ConfigType::ErrFlag)
      .def_readwrite("Cavity_on",      &tsc::ConfigType::Cavity_on)
      .def_readwrite("radiation",      &tsc::ConfigType::radiation)
      .def_readwrite("emittance",      &tsc::ConfigType::emittance)
      .def_readwrite("quad_fringe",    &tsc::ConfigType::quad_fringe)
      .def_readwrite("H_exact",        &tsc::ConfigType::H_exact)
      .def_readwrite("Cart_Bend",      &tsc::ConfigType::Cart_Bend)
      .def_readwrite("dip_edge_fudge", &tsc::ConfigType::dip_edge_fudge)
      .def_readwrite("pathlength",     &tsc::ConfigType::pathlength)
      .def_readwrite("Aperture_on",    &tsc::ConfigType::Aperture_on)
      .def_readwrite("EPU",            &tsc::ConfigType::EPU)
      .def_readwrite("mat_meth",       &tsc::ConfigType::mat_meth)
      .def_readwrite("IBS",            &tsc::ConfigType::IBS)
      .def_readwrite("tuneflag",       &tsc::ConfigType::tuneflag)
      .def_readwrite("chromflag",      &tsc::ConfigType::chromflag)
      .def_readwrite("codflag",        &tsc::ConfigType::codflag)
      .def_readwrite("mapflag",        &tsc::ConfigType::mapflag)
      .def_readwrite("passflag",       &tsc::ConfigType::passflag)
      .def_readwrite("overflag",       &tsc::ConfigType::overflag)
      .def_readwrite("chambre",        &tsc::ConfigType::chambre)
      // long int
      // int
      // anything below this line ... expect to be removed
      // reason: these data are rather summary than configuration
      // double
      .def_readwrite("dPcommon",       &tsc::ConfigType::dPcommon)
      .def_readwrite("dPparticle",     &tsc::ConfigType::dPparticle)
      .def_readwrite("Energy",         &tsc::ConfigType::Energy)
      .def_readwrite("dE",             &tsc::ConfigType::dE)
      .def(py::init<>());
}
