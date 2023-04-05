#include <thor_scsi/custom/aircoil_interpolation.h>

namespace tsc = thor_scsi::core;
namespace tsu = thor_scsi::custom;

template<class C>
void tsu::AirCoilMagneticFieldKnobbed<C>::show(std::ostream& strm, int level) const
{
    strm << "AirCoilMagneticFieldKnobbed({ ";
    bool first = true;
    for(const auto& f: this->m_filaments){
	if(!first){
	    strm << ", ";
	}
	strm << "{";
	f.show(strm, level);
	strm <<"}";
	first = false;
    }
    strm << "})";
}


template void tsu::AirCoilMagneticFieldKnobbed<tsc::StandardDoubleType>::show(std::ostream& strm, int level) const;
template void tsu::AirCoilMagneticFieldKnobbed<tsc::TpsaVariantType>::show(std::ostream& strm, int level) const;
