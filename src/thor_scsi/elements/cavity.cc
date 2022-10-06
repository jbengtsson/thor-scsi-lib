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
	this->setVoltage(config.get<double>("Voltage"));
	this->setHarmonicNumber(config.get<double>("harnum"));
}


template<typename T>
void tse::CavityType::_localPropagate(tsc::ConfigType &conf, gtpsa::ss_vect<T> &ps)
{

	const double L = this->PL, c0 = speed_of_light;
	const bool debug = false;

	drift_propagate(conf, L/2e0, ps);

	if(debug){
		std::cout << "cavity on " << conf.Cavity_on
			  <<  "   frequency " << this->Pfreq
			  <<  "   voltage " << this->Pvolt
			  <<  "radiation on " << conf.radiation

			  << std::endl;
	}

	if (conf.Cavity_on && this->Pvolt != 0e0) {
		auto energy = conf.Energy;
		if(!std::isfinite(energy)){
			throw std::runtime_error("Energy is not NaN and cavity calculation requested");
		}
		// if(np.finite(energy))
		const T delta = - this->Pvolt / (energy)
			*sin(2e0*M_PI*this->Pfreq/c0*ps[ct_]+this->phi);

		if(debug){
			std::cout << "Cavity computed delta " << delta
				  << " cavity phase " << this->phi
				  << " energy "<< conf.Energy
				  << " speed of light " << c0
				  << " M_PI" << M_PI
				  << " ct  " << ps[ct_]
				  << std::endl;
		}

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
template void tse::CavityType::_localPropagate(tsc::ConfigType &conf, gtpsa::ss_vect<tps>         &ps);
template void tse::CavityType::_localPropagate(tsc::ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps);

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
