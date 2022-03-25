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
      .def_readwrite("Cell_nLoc",      &tsc::ConfigType::Cell_nLoc)
      .def_readwrite("Elem_nFam",      &tsc::ConfigType::Elem_nFam)
      .def_readwrite("CODimax",        &tsc::ConfigType::CODimax)
      // int
      .def_readwrite("bpm",            &tsc::ConfigType::bpm)
      .def_readwrite("hcorr",          &tsc::ConfigType::hcorr)
      .def_readwrite("vcorr",          &tsc::ConfigType::vcorr)
      .def_readwrite("qt",             &tsc::ConfigType::qt)
      .def_readwrite("gs",             &tsc::ConfigType::gs)
      .def_readwrite("ge",             &tsc::ConfigType::ge)
      .def_readwrite("RingType",       &tsc::ConfigType::RingType)
      .def_readwrite("lossplane",      &tsc::ConfigType::lossplane)
      // double
      .def_readwrite("dPcommon",       &tsc::ConfigType::dPcommon)
      .def_readwrite("dPparticle",     &tsc::ConfigType::dPparticle)
      .def_readwrite("delta_RF",       &tsc::ConfigType::delta_RF)
      .def_readwrite("Omega",          &tsc::ConfigType::Omega)
      .def_readwrite("U0",             &tsc::ConfigType::U0)
      .def_readwrite("Alphac",         &tsc::ConfigType::Alphac)
      .def_readwrite("Energy",         &tsc::ConfigType::Energy)
      .def_readwrite("dE",             &tsc::ConfigType::dE)
      .def_readwrite("CODeps",         &tsc::ConfigType::CODeps)
      .def_readwrite("Qb",             &tsc::ConfigType::Qb)
      .def_readwrite("alpha_z",        &tsc::ConfigType::alpha_z)
      .def_readwrite("beta_z",         &tsc::ConfigType::beta_z)
      .def_readwrite("beta0",          &tsc::ConfigType::beta0)
      .def_readwrite("gamma0",         &tsc::ConfigType::gamma0)
      // std::vector<double>:
      .def_readwrite("TotalTune",      &tsc::ConfigType::TotalTune)
      .def_readwrite("Chrom",          &tsc::ConfigType::Chrom)
      .def_readwrite("alpha_rad",      &tsc::ConfigType::alpha_rad)
      .def_readwrite("D_rad",          &tsc::ConfigType::D_rad)
      .def_readwrite("J",              &tsc::ConfigType::J)
      .def_readwrite("tau",            &tsc::ConfigType::tau)
      .def_readwrite("D_IBS",          &tsc::ConfigType::D_IBS)
      .def_readwrite("eps",            &tsc::ConfigType::eps)
      .def_readwrite("epsp",           &tsc::ConfigType::epsp)
      .def_readwrite("CODvect",        &tsc::ConfigType::CODvect)
      .def_readwrite("wr",             &tsc::ConfigType::wr)
      .def_readwrite("wi",             &tsc::ConfigType::wi)
      // std::vector< std::vector<double> >
      .def_readwrite("OneTurnMat",     &tsc::ConfigType::OneTurnMat)
      .def_readwrite("Ascr",           &tsc::ConfigType::Ascr)
      .def_readwrite("Ascrinv",        &tsc::ConfigType::Ascrinv)
      .def_readwrite("Vr",             &tsc::ConfigType::Vr)
      .def_readwrite("Vi",             &tsc::ConfigType::Vi)
      .def(py::init<>());
}
