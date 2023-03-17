#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/elements/cavity.h>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/elements/constants.h>
#include <gtpsa/utils.hpp>

#include <cmath>
#include <ostream>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;


tse::CavityType::CavityType(const Config &config) :  LocalGalilean(config)
{
	this->setFrequency(config.get<double>("Frequency"));
	this->setPhase(config.get<double>("Phase", 0e0));
	this->setVoltage(config.get<double>("Voltage"));
	this->setHarmonicNumber(config.get<double>("HarmonicNumber"));
}


template<typename T>
void tse::CavityType::_localPropagate(tsc::ConfigType &conf, gtpsa::ss_vect<T> &ps)
{

	const double L = this->PL, c0 = speed_of_light;
	const bool debug = false;

	drift_propagate(conf, L/2e0, ps);

	THOR_SCSI_LOG(DEBUG)
		<< "\n  cavity on    = " << conf.Cavity_on
		<< "\n  f_RF         = " << this->Pfreq
		<< "\n  V_RF         = " << this->Pvolt
		<< "\n  radiation on = " << conf.radiation;

	if (conf.Cavity_on && this->Pvolt != 0e0) {
		auto energy = conf.Energy;
		if(!std::isfinite(energy)){
			throw std::runtime_error(
				"Energy is NaN and cavity calculation requested");
		}
		// if(np.finite(energy))
		const T delta = - this->Pvolt / (energy)
			*sin(2e0*M_PI*this->Pfreq/c0*ps[ct_]+this->phi);

		THOR_SCSI_LOG(DEBUG)
			<< "\n  delta:\n" << delta
			<< "\n  phi_RF = " << this->phi
			<< "\n  E      = " << conf.Energy
			<< "\n  c_0    = " << c0
			<< "\n  M_PI   = " << M_PI
			<< "\n  ct:\n" << ps[ct_];

		ps[delta_] += delta;

#ifdef THOR_SCSI_USE_RADIATION
		if (conf.radiation) conf.dE -= gtpsa::cst(delta);
#endif
		if (conf.pathlength) ps[ct_] -= this->Ph/this->Pfreq*c0;
	}
	drift_propagate(conf, L/2e0, ps);
}

// template void tse::CavityType::_localPropagate(tsc::ConfigType &conf, ss_vect<double>             &ps);
// template void tse::CavityType::_localPropagate(tsc::ConfigType &conf, ss_vect<tps>                &ps);
template void tse::CavityType::_localPropagate(tsc::ConfigType &conf, gtpsa::ss_vect<double>      &ps);
// template void tse::CavityType::_localPropagate(tsc::ConfigType &conf, gtpsa::ss_vect<tps>         &ps);
template void tse::CavityType::_localPropagate(tsc::ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps);

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
