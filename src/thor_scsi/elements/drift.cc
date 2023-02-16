#include <thor_scsi/core/config.h>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/elements/drift.h>
// #include <tps/ss_vect.h>
// #include <tps/tps.h>
#include <tps/tps_type.h>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

template<class C>
template<typename T>
void tse::DriftTypeWithKnob<C>::_propagate(const tsc::ConfigType &conf, gtpsa::ss_vect<T> &ps)
{
	// auxilliar function implementing a "non lattice element" drift
	drift_propagate(conf, this->PL, ps);

	if (conf.emittance && !conf.Cavity_on){
		// Needs A^-1.
		// curly_dH_x = is_tps<tps>::get_curly_H(ps);
		;
	}
}

template void tse::DriftType::_propagate(const tsc::ConfigType &conf, gtpsa::ss_vect<double>      &ps);
template void tse::DriftType::_propagate(const tsc::ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps);
// template void tse::DriftType::_propagate(const tsc::ConfigType &conf, gtpsa::ss_vect<tps>         &ps);

template void tse::DriftTypeTpsa::_propagate(const tsc::ConfigType &conf, gtpsa::ss_vect<double>      &ps);
template void tse::DriftTypeTpsa::_propagate(const tsc::ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps);
// template void tse::DriftTypeTpsa::_propagate(const tsc::ConfigType &conf, gtpsa::ss_vect<tps>         &ps);


/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
