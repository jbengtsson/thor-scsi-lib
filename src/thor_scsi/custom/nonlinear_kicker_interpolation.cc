#include <thor_scsi/custom/nonlinear_kicker_interpolation.h>

namespace tsc = thor_scsi::core;
namespace tsu = thor_scsi::custom;


const std::vector<tsu::aircoil_filament_t> tsu::construct_aircoil_filaments(
    const std::vector<tsu::aircoil_filament_t>& filaments_one_quater)
{

    const double n_filaments = filaments_one_quater.size() * 4;
    std::vector<tsu::aircoil_filament_t> filaments;
    filaments.reserve(n_filaments);

    for(const auto& f: filaments_one_quater){
	filaments.push_back({  f.x,  f.y,  f.current});
	filaments.push_back({ -f.x,  f.y,  f.current});
	filaments.push_back({  f.x, -f.y,  f.current});
	filaments.push_back({ -f.x, -f.y,  f.current});

    }
    return filaments;
}


template<class C>
void tsu::NonlinearKickerKnobbed<C>::show(std::ostream& strm, int level) const
{
    strm << "NonlinearKickerKnobbed({ ";
    tsu::AirCoilMagneticFieldKnobbed<C>::show(strm, level);
    strm << "})";
}

template void tsu::NonlinearKickerKnobbed<tsc::StandardDoubleType>::show(std::ostream& strm, int level) const;
template void tsu::NonlinearKickerKnobbed<tsc::TpsaVariantType>::show(std::ostream& strm, int level) const;
